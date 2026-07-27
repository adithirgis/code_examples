var dataset = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY")
.select('temperature_2m')
.filter(ee.Filter.date('2022-06-01', '2022-07-01'))
print(dataset)

// Clip to the output image to the Bengaluru Boundaries.
var clipped = dataset.mean().clip(geometry);

Map.centerObject(clipped);

var visualization = {
  bands: ['temperature_2m'],
  min: 296.0,
  max: 294.0,
  palette: [
    "#000080","#0000D9","#4000FF","#8000FF","#0080FF","#00FFFF",
    "#00FF80","#80FF00","#DAFF00","#FFFF00","#FFF500","#FFDA00",
    "#FFB000","#FFA400","#FF4F00","#FF2500","#FF0A00","#FF00FF",
  ]
};

Map.addLayer(clipped, visualization, "Air temperature [K] at 2m height");
// syntax to export the image to google drive
Export.image.toDrive({
  image: clipped, 
  description: 'Ban_Temperature_June_2022',  // image name
  scale: 11132 , //pixel size
  region: geometry,
  maxPixels: 1E10
})