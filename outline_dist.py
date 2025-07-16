var collection = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA")
  .filterBounds(studyArea)
  .filterDate('2014-06-01', '2021-09-30')
  .filter(ee.Filter.calendarRange(6, 9, 'month'))
  .filter(ee.Filter.lt('CLOUD_COVER', 40));

// Cloud mask function
function maskL8(image) {
  var qa = image.select('QA_PIXEL');
  var cloudBitMask = 1 << 3;
  var cloudShadowBitMask = 1 << 4;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
              .and(qa.bitwiseAnd(cloudShadowBitMask).eq(0));
  return image.updateMask(mask);
}

// Function to calculate salinity indices
function calculateSalinityIndices(img) {
  var bands = {
    B: img.select('B2'), G: img.select('B3'), R: img.select('B4'),
    NIR: img.select('B5'), SWIR1: img.select('B6'), SWIR2: img.select('B7')
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

  return img.addBands(indices);
}

// Process image collection
var corrected = collection.map(maskL8);
var salinityImage = corrected.map(calculateSalinityIndices).median().clip(studyArea);

// Sample the image at your labeled points
var sampledPoints = salinityImage.sampleRegions({
  collection: salinitySamples,
  properties: ['salinity_class'],  // Must match your shapefile attribute name
  scale: 30
});

// Export the sample points with salinity indices to CSV
Export.table.toDrive({
  collection: sampledPoints,
  description: 'Salinity_Sample_Export',
  fileFormat: 'CSV',
  folder: 'Soil_Salinity'
});