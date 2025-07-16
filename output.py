

var image = ee.Image("projects/ee-harshav1575/assets/salinity_prediction_map_1");

Map.centerObject(image, 9); 

var visParams = {
  min: 0,      
  max: 1,       
  palette: ['red', 'green'] 
};

var clipped = image.clip(region);
Map.centerObject(region, 9);
Map.addLayer(clipped, visParams, 'Clipped Salinity Map');
