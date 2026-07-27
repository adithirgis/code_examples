var dataset_u = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY")
.select('u_component_of_wind_10m')
.filter(ee.Filter.date('2022-06-01', '2022-07-01'))

var dataset_v = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY")
.select('v_component_of_wind_10m')
.filter(ee.Filter.date('2022-06-01', '2022-07-01')) 
//  print(dataset)
print(dataset_u);
print(dataset_v); 

// Clip to the output image to the Bengaluru Boundaries.
var clipped_u = dataset_u.mean().clip(geometry); 
var clipped_v = dataset_v.mean().clip(geometry); 
  
var visualization = {
  bands: ['u_component_of_wind_10m'],
  min: -1.5,
  max: -0.9,
  palette: [
    "#000080","#0000D9","#4000FF","#8000FF","#0080FF","#00FFFF",
    "#00FF80","#80FF00","#DAFF00","#FFFF00","#FFF500","#FFDA00",
    "#FFB000","#FFA400","#FF4F00","#FF2500","#FF0A00","#FF00FF",
  ]
};

var vis = {
  bands: ['v_component_of_wind_10m'],
  min: 0.1,
  max: 0.6,
  palette: [
    "#000080","#0000D9","#4000FF","#8000FF","#0080FF","#00FFFF",
    "#00FF80","#80FF00","#DAFF00","#FFFF00","#FFF500","#FFDA00",
    "#FFB000","#FFA400","#FF4F00","#FF2500","#FF0A00","#FF00FF",
  ]
};

Map.addLayer(clipped_u, visualization, "u_component_of_wind_10m");
Map.addLayer(clipped_v, vis, "v_component_of_wind_10m");

Export.image.toDrive({
  image: clipped_u, 
  description: 'Ban_u_wind_June_2022',  
  scale: 11132 ,
  region: geometry,
  maxPixels: 1E10
})

Export.image.toDrive({
  image: clipped_v, 
  description: 'Ban_v_wind_June_2022',  
  scale: 11132 , 
  region: geometry,
  maxPixels: 1E10
})

var start_period = ee.Date('2022-06-01');
var end_period = ee.Date('2022-07-01');

var ERA5 = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY")
.filter(ee.Filter.date(start_period, end_period))

var ERA5NL=ERA5.map(function(image){ 
  return image.clip(geometry)});

var ERA5windspeed = ERA5NL.map(function(image){
  var wind_10m = image.expression(
    'sqrt(u**2 + v**2)', {
      'u': image.select('u_component_of_wind_10m'),
      'v': image.select('v_component_of_wind_10m')
    }).rename('windspeed');
  var time = image.get('system:time_start');
  return wind_10m.set('system:time_start', time) } );

var ERA5meanspeed = ERA5windspeed.map(function(image){
  var meandict = image.select('windspeed').reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geometry
  })
  return image.set(meandict); 
});

var visualization = {
  min: 0.9,
  max: 1.7,
  palette: [
    "#000080","#0000D9","#4000FF","#8000FF","#0080FF","#00FFFF",
    "#00FF80","#80FF00","#DAFF00","#FFFF00","#FFF500","#FFDA00",
    "#FFB000","#FFA400","#FF4F00","#FF2500","#FF0A00","#FF00FF",
  ]
};

Map.addLayer(ERA5meanspeed, visualization, "windspeed10m");

Export.image.toDrive({
  image: ERA5meanspeed.mean(), 
  description: 'Ban_windspeed_June_2022',  
  scale: 11132 , 
  region: geometry,
  maxPixels: 1E10
  
})