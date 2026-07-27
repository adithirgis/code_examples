var admin1 = ee.FeatureCollection("FAO/GAUL_SIMPLIFIED_500m/2015/level1");

var admin1 = ee.FeatureCollection("FAO/GAUL/2015/level2");
var bangalore1 = admin1.filter(ee.Filter.eq('ADM2_NAME','Bangalore Urban'));

var geometry_ban = bangalore1.geometry()
var dataset = ee.ImageCollection("ESA/WorldCover/v100").first();

// Clip to the output image to the Bengaluru Boundaries.
var clipped = dataset.clip(geometrybox);

var visualization = {
  bands: ['Map'],
};

Map.centerObject(clipped);
Map.addLayer(clipped, visualization, "Landcover");

Export.image.toDrive({
  image: clipped, 
  description: 'Ban_LULC_2020', 
  scale: 10, 
  region: geometrybox,
  maxPixels: 1E10
})