var dataset = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY")
.select('temperature_2m')
.filter(ee.Filter.date('2019-01-01', '2019-12-31'))
print(dataset)

// Clip to the output image to the Delhi Boundaries
var clipped = dataset.mean().clip(geometrybox);

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
  description: 'Delhi_Temperature_2019',  // image name
  scale: 11132 , //pixel size
  region: geometrybox,
  maxPixels: 1E10
})