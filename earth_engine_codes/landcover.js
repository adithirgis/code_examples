var admin1 = ee.FeatureCollection("FAO/GAUL_SIMPLIFIED_500m/2015/level1");

var admin1 = ee.FeatureCollection("FAO/GAUL/2015/level2");
var delhi = admin1.filter(ee.Filter.eq('ADM2_NAME','Delhi'));

var geometry_delhi = delhi.geometry()
var dataset = ee.ImageCollection("ESA/WorldCover/v100").first();

// Clip to the output image to the Delhi Boundaries.
var clipped = dataset.clip(geometrybox);

var visualization = {
  bands: ['Map'],
};

Map.centerObject(clipped);
Map.addLayer(clipped, visualization, "Landcover");

Export.image.toDrive({
  image: clipped, 
  description: 'Delhi_LULC_2020', 
  scale: 10, 
  region: geometrybox,
  maxPixels: 1E10
})