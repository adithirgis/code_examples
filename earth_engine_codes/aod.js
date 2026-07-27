var dataset = ee.Image('USGS/SRTMGL1_003');
var elevation = dataset.select('elevation');

// Clip to the output image to the Bengaluru Boundaries.
var clipped = elevation.clip(geometrybox);

Map.addLayer(clipped, { min: 0, max: 5000}, "clipped");
Map.centerObject(clipped);

Export.image.toDrive({
  image: clipped, 
  description: 'Ban_SRTM_Elevation',  
  scale: 30, 
  region: geometrybox,
  maxPixels: 1E10
})