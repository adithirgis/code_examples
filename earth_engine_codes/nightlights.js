// import the nightlight image collection
var nightlight = ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG");

// set the start and end dates
var start = ee.Date.fromYMD(2019,01,1);
var end = ee.Date.fromYMD(2019,12,31);

// filter for the period of interest
var nightlights2019 = nightlight.filterDate(start,end);
// take the mean note that this operation transforms the image collection into an image
nightlights2019 = ee.Image(nightlights2019.mean());
// select the avg_rad band
nightlights2019 = nightlights2019.select("avg_rad");
// clip for the area of interest
nightlights2019 = nightlights2019.clip(geometrybox);
Map.addLayer(nightlights2019,{min:0,max:100,palette:['000000','700000','808080','FFFF00','ffffff','ffffff','ffffff']},"nightlights 2021");

Export.image.toDrive({
  image: nightlights2019, 
  description: 'Delhi_NTLI_2019',  
  scale: 463.83, 
  region: geometrybox,
  maxPixels: 1E10
})
