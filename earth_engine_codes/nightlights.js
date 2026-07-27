// import the nightlight image collection
var nightlight = ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG");

// set the start and end dates
var start = ee.Date.fromYMD(2021,12,1);
var end = ee.Date.fromYMD(2022,05,31);

// filter for the period of interest
var nightlights2021 = nightlight.filterDate(start,end);
// take the mean note that this operation transforms the image collection into an image
nightlights2021 = ee.Image(nightlights2021.mean());
// select the avg_rad band
nightlights2021 = nightlights2021.select("avg_rad");
// clip for the area of interest
nightlights2021 = nightlights2021.clip(geometrybox);
Map.addLayer(nightlights2021,{min:0,max:100,palette:['000000','700000','808080','FFFF00','ffffff','ffffff','ffffff']},"nightlights 2021");

Export.image.toDrive({
  image: nightlights2021, 
  description: 'Ban_NTLI_Dec2021_May2022',  
  scale: 463.83, 
  region: geometrybox,
  maxPixels: 1E10
})
