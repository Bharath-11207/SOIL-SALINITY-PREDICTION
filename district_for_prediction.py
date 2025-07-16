// Center map and display boundaries
Map.centerObject(studyArea, 10);
Map.addLayer(studyArea, {}, 'Study Area');

// Function to mask clouds
var cloudMask = function(img) {
  var q = img.select('QA_PIXEL');
  var cloud = 1 << 3;
  return img.updateMask(q.bitwiseAnd(cloud).eq(0))
           .multiply(0.0000275)
           .add(-0.2);
};

// Landsat-8 Composite (2020-2023)
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(studyArea)
  .filterDate('2023-01-01', '2023-01-10')
  .filterMetadata('CLOUD_COVER', 'less_than', 40)
  .map(cloudMask);

var composite = l8.median().clip(studyArea);

// Function to add salinity indices
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

  return image.addBands(indices)
              .select(['SI1','SI2','SI3','SI4','SI5','SI6','SI7','SI8','SI9','NDSI','CRSI','NDVI']);
};

var indexImage = addIndices(composite);

// Visual check
Map.addLayer(indexImage.select('NDVI'), {min: -1, max: 1, palette: ['red', 'yellow', 'green']}, 'NDVI');
Map.addLayer(indexImage.select('SI1'), {min: 0, max: 5}, 'SI1');

// Export full composite with all indices
Export.image.toDrive({
  image: indexImage.toFloat(),
  description: 'Srikakulam_12Band_Composite_232',
  folder: 'DL',
  region: studyArea,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
