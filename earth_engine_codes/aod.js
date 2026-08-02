// load Delhi boundary
var geometrybox = delhi

// import data
var AOD = ee.ImageCollection("MODIS/006/MCD19A2_GRANULES");

var modisBands = ['Optical_Depth_055'];

// helper function to extract the QA bits
function getQABits(image, start, end, newName) {
  // Compute the bits we need to extract
  var pattern = 0;
  for (var i = start; i <= end; i++) {
    pattern += Math.pow(2, i);
  }
  // Return a single band image of the extracted QA bits, giving the band a new name
  return image.select([0], [newName])
  .bitwiseAnd(pattern)
  .rightShift(start);
}

// A function to mask out cloudy pixels.
function maskQuality(image) {
  // Select the QA band.
  var QA = image.select('AOD_QA');
  // Get the internal_cloud_algorithm_flag bit.
  var internalQuality = getQABits(QA,0, 2, 'internal_quality_flag');
  // Return an image masking out cloudy areas.
  return image.updateMask(internalQuality.eq(1));
}

// create cloud free composite
var AODmaskQ = AOD.filter(ee.Filter.or(
  ee.Filter.date('2019-01-01', '2019-12-31')
))
.map(maskQuality)
.select(modisBands)
.filterBounds(geometrybox);

print(AODmaskQ);

// create composite with quality assurance (without clouds) 
var AODwithoutmask = AOD.filter(ee.Filter.or(
  ee.Filter.date('2019-01-01', '2019-12-31')
))
.select(modisBands)
.filterBounds(geometrybox);


// vis parameters
var viz = {
  min: 0,
  max:1.2,
  bands:['Optical_Depth_055'],
  palette: ['black', 'blue', 'purple', 'cyan', 'green', 'yellow', 'red']
};
var composite1 = AODmaskQ.mean().multiply(0.001).clip(geometrybox)
var composite2 = AODwithoutmask.mean().multiply(0.001).clip(geometrybox)

print (composite1)
// add the cloud free composite
Map.addLayer(composite1,viz,'Quality mask');
// add the cloud composite
Map.addLayer(composite2,viz,'Without mask');

Export.image.toDrive({
  image:composite1, 
  description: 'Delhi_AOD_2019',
  scale: 1000,
  folder:"LUR Predictors GEE",
  region:geometrybox
  
});