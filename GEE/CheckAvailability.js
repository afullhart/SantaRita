/////////////////////////////
//Check data availability for June and Sept 2019 for roi
///////////////////////////


//Why .median()?
//During May 2019, Sentinel-2 likely captured multiple images over your ROI. 
//By calling .median(), Earth Engine looks at every overlapping pixel across all those dates and takes the median value. 
//This is a highly effective way to automatically remove transient clouds or shadows from your final visualization!


// 1. Define the Region of Interest
var roi = ee.Geometry.BBox(-111.00, 31.70, -110.75, 31.95);

Map.centerObject(roi, 11);
Map.addLayer(roi, {color: 'blue'}, 'Region of Interest (ROI)', false);

// 2. Define Datasets and Time Periods
var datasets = {
  'Landsat 8': 'LANDSAT/LC08/C02/T1_L2',
  'Sentinel-2': 'COPERNICUS/S2_SR_HARMONIZED',
  'NAIP': 'USDA/NAIP/DOQQ'
};

var periods = {
  'May 2019': ['2019-05-01', '2019-06-01'],
  'September 2019': ['2019-09-01', '2019-10-01']
};

print('--- Spatial Coverage Report ---');
print('Calculating true footprint using the Binary Mean method...');

// 3. Loop through datasets and periods to check coverage
Object.keys(datasets).forEach(function(dsName) {
  Object.keys(periods).forEach(function(periodName) {
    var id = datasets[dsName];
    var dates = periods[periodName];

    var col = ee.ImageCollection(id)
      .filterBounds(roi)
      .filterDate(dates[0], dates[1]);

    col.size().evaluate(function(size) {
      if (size === 0) {
        print(dsName + ' (' + periodName + '): 0.0% coverage ⚠️ PARTIAL (0 images found)');
        return;
      }

      var binaryCoverage = col.select(0).count().unmask(0).gt(0).rename('coverage');

      // Calculate the mean value of the 1s and 0s inside the ROI.
      var stats = binaryCoverage.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: roi,
        scale: 30, // 30m resolution for the check
        maxPixels: 1e9
      });

      // Multiply the mean fraction by 100 to get the percentage
      var coveragePct = ee.Number(stats.get('coverage')).multiply(100);

      // Pull the percentage to the client side
      coveragePct.evaluate(function(pct) {
        var finalPct = pct || 0; 
        var status = finalPct >= 99.5 ? '✅ COMPLETE' : '⚠️ PARTIAL';
        print(dsName + ' (' + periodName + '): ' + finalPct.toFixed(1) + '% coverage ' + status + ' (' + size + ' images)');
      });
    });
  });
});



///////////////////////////////////////////////////////////
//Add filtered image collection to check with the UI widget
//////////////////////////////////////////////////////////

var sent2_ic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)                          // FILTER LOCATION FIRST
  .filterDate('2019-05-01', '2019-05-31')     // FILTER DATE SECOND
  .select(['B4', 'B3', 'B2']);

var trueColorImage = sent2_ic.median().clip(roi);

var visParams = {
  bands: ['B4', 'B3', 'B2'], // Ensure Red, Green, Blue order
  min: 0,
  max: 3000                  // Stretches the brightness so it isn't too dark
};

Map.addLayer(trueColorImage, visParams, 'Sentinel-2 True Color');



///////////////////////////////////////////////////
//UI that prints system:id of images upon clicking
/////////////////////////////////////////////////

var inspectorPanel = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px',
    width: '350px'
  }
});

Map.add(inspectorPanel);
inspectorPanel.add(ui.Label('Click anywhere on the map to get the NAIP image ID.'));

Map.onClick(function(coords) {
  
  inspectorPanel.clear();
  inspectorPanel.add(ui.Label('Fetching ID...'));

  var clickPoint = ee.Geometry.Point(coords.lon, coords.lat);

  var intersectingImages = sent2_ic.filterBounds(clickPoint);

  var topImageId = intersectingImages.first().get('system:id');

  topImageId.evaluate(function(idString) {
    inspectorPanel.clear();
    
    if (idString) {
      inspectorPanel.add(ui.Label('Image ID:', {fontWeight: 'bold'}));
      inspectorPanel.add(ui.Label(idString, {fontSize: '11px', color: 'blue'}));
    } else {
      inspectorPanel.add(ui.Label('No image found at this exact pixel.'));
    }
  });
});

Map.style().set('cursor', 'crosshair');


