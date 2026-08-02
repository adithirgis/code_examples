// set of points in Delhi (same buffers as your NDVI script)
var pointList = ee.List(table);

// generate the buffers
var bufferBy = function(size) {
  return function(feature) {
    return feature.buffer(size);
  };
};

var getCentroids = function(feature) {
  return feature.set({polyCent: feature.centroid()});
};

var pointCollectionBfs100  = table.map(bufferBy(100)).map(getCentroids);
var pointCollectionBfs300  = table.map(bufferBy(300)).map(getCentroids);
var pointCollectionBfs500  = table.map(bufferBy(500)).map(getCentroids);
var pointCollectionBfs1000 = table.map(bufferBy(1000)).map(getCentroids);
var pointCollectionBfs5000     = table.map(bufferBy(5000)).map(getCentroids);

Map.addLayer(pointCollectionBfs5000, {color: 'red'}, 'buffer_5000m');
Map.addLayer(pointCollectionBfs1000, {color: 'yellow'}, 'buffer_1000m');
Map.addLayer(pointCollectionBfs500, {color: 'white'}, 'buffer_500m');
Map.addLayer(pointCollectionBfs300, {color: 'green'}, 'buffer_300m');
Map.addLayer(pointCollectionBfs100, {color: 'black'}, 'buffer_100m');

var coll  = pointCollectionBfs5000;
var coll1 = pointCollectionBfs1000;
var coll2 = pointCollectionBfs500;
var coll3 = pointCollectionBfs300;
var coll4 = pointCollectionBfs100;

// ESA WorldCover
var esa = ee.ImageCollection("ESA/WorldCover/v100").first();

// Class values (ESA WorldCover codes) and matching names
var classValues = [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100];
var classNames = [
  'tree_cover', 'shrubland', 'grassland', 'cropland', 'builtup',
  'bare_land', 'snow_ice', 'per_water_bod', 'herb_wetland',
  'mang', 'moss_lichen'
];

// pixel area in m^2
var pixelArea = ee.Image.pixelArea();

// one band per class = area (m^2) of that class at each pixel
var classAreaBands = classValues.map(function(val, i) {
  return esa.eq(val).multiply(pixelArea).rename(classNames[i]);
});

var classAreaImage = ee.Image.cat(classAreaBands);

// Sum area per class within each buffer
var pts_lulc5000 = classAreaImage.reduceRegions({
  collection: coll,
  reducer: ee.Reducer.sum(), 
  scale: 10
});

var pts_lulc1000 = classAreaImage.reduceRegions({
  collection: coll1,
  reducer: ee.Reducer.sum(),
  scale: 10
});

var pts_lulc500 = classAreaImage.reduceRegions({
  collection: coll2,
  reducer: ee.Reducer.sum(),
  scale: 10
});

var pts_lulc300 = classAreaImage.reduceRegions({
  collection: coll3,
  reducer: ee.Reducer.sum(),
  scale: 10
});

var pts_lulc100 = classAreaImage.reduceRegions({
  collection: coll4,
  reducer: ee.Reducer.sum(),
  scale: 10
});

print(pts_lulc100); // check the class-area properties are attached

// selectors must list CODE, Location + every class values band
var selectorsList = ['CODE', 'Location'].concat(classValues);

Export.table.toDrive({
  collection: pts_lulc100,
  description: 'Delhi_LULC_100m_2019',
  selectors: selectorsList,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: pts_lulc300,
  description: 'Delhi_LULC_300m_2019',
  selectors: selectorsList,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: pts_lulc500,
  description: 'Delhi_LULC_500m_2019',
  selectors: selectorsList,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: pts_lulc1000,
  description: 'Delhi_LULC_1000m_2019',
  selectors: selectorsList,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: pts_lulc5000,
  description: 'Delhi_LULC_5000m_2019',
  selectors: selectorsList,
  fileFormat: 'CSV'
});