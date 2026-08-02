var collection = ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_NO2')
.select('tropospheric_NO2_column_number_density')
.filterDate('2019-01-01', '2019-12-31');

var band_viz = {
  min: 0,
  max: 0.00009,
  palette: ['black', 'blue', 'purple', 'cyan', 'green', 'yellow', 'red']
};

// Create a geometry representing an export region.
var geometrybox = delhi;

var no2 = collection.mean().clip(geometrybox);

Map.addLayer(no2, band_viz, 'S5P N02');
Map.setCenter(77.27, 13.11, 8);

Export.image.toDrive({
  image: no2,
  description: 'Delhi_Tropospheric_NO2_Sentinel5P_2019',
  scale: 1113.2,
  region: geometrybox
});