// set of points
var pointList= ee.List(table)

// generate the buffers
var bufferBy = function(size) {
  return function(feature) {
    return feature.buffer(size);   
  };
};

// adds the centroid to the buffer as a feature
var getCentroids = function(feature) {
  return feature.set({polyCent: feature.centroid()});
};

var pointCollectionBfs100 = table.map(bufferBy(100))
.map(getCentroids);

var pointCollectionBfs300 = table.map(bufferBy(300))
.map(getCentroids);

var pointCollectionBfs500 = table.map(bufferBy(500))
.map(getCentroids);

var pointCollectionBfs1000 = table.map(bufferBy(1000))
.map(getCentroids);

var pointCollectionBfs = table.map(bufferBy(5000))
.map(getCentroids);

Map.addLayer(pointCollectionBfs, {color: 'red'}, 'buffers5km');
Map.addLayer(pointCollectionBfs1000, {color: 'yellow'}, 'buffers1km');
Map.addLayer(pointCollectionBfs500, {color: 'white'}, 'buffers500m');
Map.addLayer(pointCollectionBfs300, {color: 'green'}, 'buffers300m');
Map.addLayer(pointCollectionBfs100, {color: 'black'}, 'buffers100m');

var coll = pointCollectionBfs ;
var coll1=pointCollectionBfs1000;
var coll2=pointCollectionBfs500;
var coll3=pointCollectionBfs300;
var coll4=pointCollectionBfs100;

var collName='pointBfs' ;

// Calculating NDVI
var s2 = ee.ImageCollection("COPERNICUS/S2");
var rgbVis = {min: 0.0, max: 3000, bands: ['B4', 'B3', 'B2']};

//Step 2: Keep only the relevant bands and filter for cloud coverage.
//Filter to metadata less than the given value ee.Filter.lt means.
var filtered = s2.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
.filter(ee.Filter.date('2022-06-01', '2022-06-30'))
//  .filter(ee.Filter.bounds(geometrybox))

print(filtered.size()); // to understand the no of data available days
//Reduces an image collection by calculating the mean of all values at each pixel across the stack of all matching bands. Bands are matched by name.
var composite = filtered.mean()
//.clip(geometrybox)
//Map.addLayer(composite, rgbVis, 'Cloudy pixel')  

// Calculate  Normalized Difference Vegetation Index (NDVI)
// 'NIR' (B8) and 'RED' (B4)

// Write a function that computes NDVI for an image and adds it as a band
function addNDVI(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('ndvi');
  return image.addBands(ndvi);
}

// Map the function over the collection
//If we want to apply some computation - such as calculating an index - to many images, you need to use map()
var withNdvi = filtered.map(addNDVI);
//print(withNdvi)

var composite = withNdvi.mean()
var ndviComposite = composite.select('ndvi')

var pts_ndvi5000 = ndviComposite.reduceRegions({
  collection: coll,
  reducer: ee.Reducer.mean(),
  scale: 10
})

var pts_ndvi1000 = ndviComposite.reduceRegions({
  collection: coll1,
  reducer: ee.Reducer.mean(),
  scale: 10
})

var pts_ndvi500 = ndviComposite.reduceRegions({
  collection: coll2,
  reducer: ee.Reducer.mean(),
  scale: 10
})

var pts_ndvi300 = ndviComposite.reduceRegions({
  collection: coll3,
  reducer: ee.Reducer.mean(),
  scale: 10
})

var pts_ndvi100 = ndviComposite.reduceRegions({
  collection: coll4,
  reducer: ee.Reducer.mean(),
  scale: 10
})

//Combine above NDVI buffer values into 1 single variable containing the data within each line set.
var ndvii = ee.FeatureCollection([
  pts_ndvi100,
  pts_ndvi300,
  pts_ndvi500,
  pts_ndvi1000,
  pts_ndvi5000
]).flatten()

print(ndvii)

Export.table.toDrive({
  collection: ndvii,
  description: 'Buffer_63sites_NDVI_June_2022',
  selectors: ['CODE','Location','mean'],
  fileFormat: 'CSV'
})