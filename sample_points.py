
Map.centerObject(Sample, 10);
Map.addLayer(Sample, {}, 'Study Area');
Map.addLayer(Sample, {}, 'Sample Points');

// Cloud masking
var cloudMask = function(img) {
  var q = img.select('QA_PIXEL');
  var cloud = 1 << 3;
  return img.updateMask(q.bitwiseAnd(cloud).eq(0))
           .multiply(0.0000275)
           .add(-0.2);
};

// Landsat composite (2020–2023)
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(Sample)
  .filterDate('2020-06-01', '2023-09-30')
  .filterMetadata('CLOUD_COVER', 'less_than', 40)
  .map(cloudMask);

var composite = l8.median().clip(Sample);

// Add salinity indices
var addIndices = function(image) {
  var bands = {
    B: image.select('SR_B2'),
    G: image.select('SR_B3'),
    R: image.select('SR_B4'),
    NIR: image.select('SR_B5'),
    SWIR1: image.select('SR_B6'),
    SWIR2: image.select('SR_B7')
  };

  var indices = ee.Image.cat([
    bands.NIR.divide(bands.SWIR1).rename('SI1'),
    bands.G.multiply(bands.R).sqrt().rename('SI2'),
    bands.G.pow(2).multiply(bands.R.pow(2)).multiply(bands.NIR.pow(2)).sqrt().rename('SI3'),
    bands.G.pow(2).multiply(bands.R.pow(2)).sqrt().rename('SI4'),
    bands.SWIR1.divide(bands.SWIR2).rename('SI5'),
    bands.NIR.subtract(bands.SWIR1).divide(bands.NIR.add(bands.SWIR1)).rename('SI6'),
    bands.SWIR1.subtract(bands.SWIR2).divide(bands.SWIR1.add(bands.SWIR2)).rename('SI7'),
    bands.B.divide(bands.R).rename('SI8'),
    bands.B.subtract(bands.R).divide(bands.B.add(bands.R)).rename('SI9'),
    bands.R.subtract(bands.NIR).divide(bands.R.add(bands.NIR)).rename('NDSI'),
    (bands.NIR.multiply(bands.R).subtract(bands.G.multiply(bands.B)))
      .divide(bands.NIR.multiply(bands.R).add(bands.G.multiply(bands.B))).sqrt().rename('CRSI'),
    bands.NIR.subtract(bands.R).divide(bands.NIR.add(bands.R)).rename('NDVI')
  ]);

  return image.addBands(indices);
};

composite = addIndices(composite);

// Create grid tiles for tiling the study area
var grid = Sample.geometry().coveringGrid(ee.Projection('EPSG:4326').atScale(7000));
Map.addLayer(grid, {}, 'Grid Tiles');

// Loop over first few tiles (adjust "i < N" as needed)
var tileList = grid.toList(grid.size());
print(grid.size())
for (var i = 0; i < 110; i++) { // Try 4 tiles first. Increase if needed.
  var tile = ee.Feature(tileList.get(i)).geometry();

  // Clip composite and filter samples to tile
  var compositeTile = composite.clip(tile);
  var samplesTile = Sample.filterBounds(tile);

  var sampledTile = compositeTile.sampleRegions({
    collection: samplesTile,
    properties: ['salinity'], // Replace with your salinity column name
    scale: 30
  });

  // Export sample table
  Export.table.toDrive({
    collection: sampledTile,
    description: 'Salinity_Samples_Tile_' + i,
    folder: 'DL',
    fileFormat: 'CSV',
    selectors: ['SI1', 'SI2', 'SI3', 'SI4', 'SI5', 'SI6', 'SI7', 'SI8', 'SI9', 'NDSI', 'CRSI', 'NDVI', 'salinity']
  });

  // Export clipped image
  Export.image.toDrive({
    image: compositeTile.toFloat(),
    description: 'Salinity_Composite_Tile_' + i,
    folder: 'DL',
    region: tile,
    scale: 60,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
}
