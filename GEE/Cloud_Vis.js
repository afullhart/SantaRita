// =========================================================================
// SETUP & ASSETS
// =========================================================================
var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm');
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');

// Extract the raw Geometry from the first feature directly
var bounds_geom = bounds_fc.first().geometry();
var v_extent = bounds_geom.bounds();

// The exact visualization parameters you requested
var rgbVis = {
  bands: ['B4', 'B3', 'B2'],
  min: 0.0,
  max: 0.3, 
  gamma: 1.4
};

// =========================================================================
// CORE LOGIC: FIND BEST WINDOW ON-THE-FLY
// =========================================================================

function getBestWindowForMonthYear(y, m) {
  y = ee.Number(y);
  m = ee.Number(m);
  
  var refDate = ee.Date.fromYMD(y, m, 1);
  var daysInMonth = refDate.advance(1, 'month').difference(refDate, 'day');
  var startDays = ee.List.sequence(1, daysInMonth.subtract(6));
  
  var windows = startDays.map(function(d) {
    var startDay = ee.Number(d);
    var startDate = ee.Date.fromYMD(y, m, startDay);
    var endDateFilter = startDate.advance(7, 'day');
    var endDateInclusive = startDate.advance(6, 'day');
    
    var s2_window = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
      .filterBounds(bounds_geom)
      .filterDate(startDate, endDateFilter);
      
    var imgCount = s2_window.size();
    var meanCldImg = s2_window.select('MSK_CLDPRB').mean();
    
    var meanCld = ee.Algorithms.If(
      imgCount.eq(0),
      100, 
      meanCldImg.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: bounds_geom,
        scale: 60,
        maxPixels: 1e9
      }).get('MSK_CLDPRB')
    );
    
    meanCld = ee.Algorithms.If(ee.Algorithms.IsEqual(meanCld, null), 100, meanCld);
    
    var windowLabel = ee.String(y.format('%d')).cat('-')
                .cat(m.format('%02d')).cat('-')
                .cat(startDay.format('%02d')).cat(' to ')
                .cat(endDateInclusive.format('%02d'));
    
    return ee.Feature(null, {
      'Start_Date': startDate.format('YYYY-MM-dd'),
      'Window_Label': windowLabel,
      'Mean_Cloud_Prob': meanCld,
      'Image_Count': imgCount
    });
  });
  
  var windowsFc = ee.FeatureCollection(windows);
  return windowsFc.sort('Mean_Cloud_Prob').first();
}


// =========================================================================
// UI APP & INTERACTIVITY
// =========================================================================

// Create a main panel to hold the UI
var panel = ui.Panel({
  style: {width: '350px', padding: '15px'}
});
ui.root.insert(0, panel);

// UI Elements
var title = ui.Label('Optimal S2 Window Explorer', {fontWeight: 'bold', fontSize: '20px'});
var desc = ui.Label('Select a Year and Month. The script will find the 7-day window with the absolute lowest cloud probability and display it.');

var yearLabel = ui.Label('Select Year:', {fontWeight: 'bold'});
var yearSlider = ui.Slider({
  min: 2018, max: 2025, value: 2019, step: 1,
  style: {stretch: 'horizontal'}
});

var monthLabel = ui.Label('Select Month:', {fontWeight: 'bold'});
var monthSlider = ui.Slider({
  min: 1, max: 12, value: 5, step: 1,
  style: {stretch: 'horizontal'}
});

var statusBox = ui.Label({
  value: 'Ready. Move sliders to calculate...',
  style: {color: 'blue', margin: '20px 0', whiteSpace: 'pre-wrap'}
});

// Add elements to the panel
panel.add(title);
panel.add(desc);
panel.add(yearLabel).add(yearSlider);
panel.add(monthLabel).add(monthSlider);
panel.add(statusBox);

// Center map on the study area
Map.centerObject(bounds_geom, 13);

// ---> CORRECTED DEFAULT MAP HERE <---
Map.setOptions('ROADMAP'); 

// The main function that runs when sliders change
function updateMap() {
  var y = yearSlider.getValue();
  var m = monthSlider.getValue();
  
  statusBox.setValue('Calculating optimal window...\n(This takes about 2-3 seconds)');
  statusBox.style().set('color', 'orange');
  
  // Ask the server for the best window feature
  var bestWindowFeature = getBestWindowForMonthYear(y, m);
  
  // .evaluate() pulls the result from the Google server back to your browser UI
  bestWindowFeature.evaluate(function(feature) {
    var props = feature.properties;
    
    if (props.Image_Count === 0) {
      statusBox.setValue('No imagery found for ' + y + '-' + m + '.\n(Try a different date)');
      statusBox.style().set('color', 'red');
      Map.layers().reset(); 
      Map.addLayer(bounds_fc.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'SRER Bounds');
      return;
    }
    
    // Update the UI Text
    var statusText = 'SUCCESS!' + 
                     '\nOptimal Dates: ' + props.Window_Label + 
                     '\nCloud Probability: ' + props.Mean_Cloud_Prob.toFixed(2) + '%' +
                     '\nImages in Window: ' + props.Image_Count;
    statusBox.setValue(statusText);
    statusBox.style().set('color', 'green');
    
    // Reconstruct the image collection for the optimal dates
    var startDate = ee.Date(props.Start_Date);
    var endDateFilter = startDate.advance(7, 'day');
    
    var s2_collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
      .filterBounds(bounds_geom)
      .filterDate(startDate, endDateFilter);
      
    // Mosaic, clip, and SCALE to match your 0.0-0.3 visualization parameters
    var s2_mosaic = s2_collection
      .mosaic()
      .clip(v_extent)
      .multiply(0.0001); 
      
    // Update the Map Layers
    Map.layers().reset(); // Clear previous layers
    
    // Layer 0: The true color composite
    Map.layers().set(0, ui.Map.Layer(s2_mosaic, rgbVis, 'Optimal S2: ' + props.Window_Label));
    
    // Layer 1: An empty red outline for your bounds
    Map.layers().set(1, ui.Map.Layer(bounds_fc.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'SRER Bounds'));
  });
}

// Attach the update function to the sliders
yearSlider.onChange(updateMap);
monthSlider.onChange(updateMap);

// Run it once on script start to load the default (May 2019)
updateMap();
