var dataset = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY")
.select('total_precipitation')
.filter(ee.Filter.date('2022-06-01', '2022-07-01'))
print(dataset)

var clipped = dataset.sum().clip(geometry);

Map.centerObject(clipped);

var visualization = {
  bands: ['total_precipitation'],
  min: 0.009,
  max: 0.015,
  palette: [
    "#000080","#0000D9","#4000FF","#8000FF","#0080FF","#00FFFF",
    "#00FF80","#80FF00","#DAFF00","#FFFF00","#FFF500","#FFDA00",
    "#FFB000","#FFA400","#FF4F00","#FF2500","#FF0A00","#FF00FF",
  ]
};

Map.addLayer(clipped, visualization, "total_precipitation");

Export.image.toDrive({
  image: clipped, 
  description: 'Ban_TotalPrecipitation_June_2022',
  scale: 11132 , 
  region: geometry,
  maxPixels: 1E10
})