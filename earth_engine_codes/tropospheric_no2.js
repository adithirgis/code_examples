var collection = ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_NO2')
.select('tropospheric_NO2_column_number_density')
.filterDate('2022-03-01', '2022-06-30');

var band_viz = {
  min: 0,
  max: 0.00009,
  palette: ['black', 'blue', 'purple', 'cyan', 'green', 'yellow', 'red']
};

// Create a geometry representing an export region.
var geometry1 = ee.Geometry.Rectangle([77.06675, 12.54762, 78.104963, 13.692457]);

var no2 = collection.mean().clip(geometry1);

Map.addLayer(no2, band_viz, 'S5P N02');
Map.setCenter(77.27, 13.11, 8);

Export.image.toDrive({
  image: no2,
  description: 'BNG_Tropospheric_NO2_Sentinel5P_2022June',
  scale: 1113.2,
  region: geometry1
});