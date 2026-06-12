// =========================================================================
// SETUP & ASSETS
// =========================================================================
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_geom = bounds_fc.first().geometry();
var v_extent = bounds_geom.bounds();

// Load the PRE-CALCULATED Landsat Cloud Windows Asset
var cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass_Landsat');

// Landsat Collection 2 SR visualization parameters
var rgbVis = {
  bands: ['red', 'green', 'blue'], 
  min: 0.0,
  max: 0.3, 
  gamma: 1.2
};

// =========================================================================
// CONSOLE CHART: TIME-SERIES OF ALL OPTIMAL WINDOWS FROM ASSET
// =========================================================================

// Ensure the collection is sorted chronologically
var sortedWindows = cloud_windows.sort('system:time_start');

// Generate the Line Chart directly from the asset
var cloudChart = ui.Chart.feature.byFeature({
  features: sortedWindows,
  xProperty: 'system:time_start',
  yProperties: ['Mean_Haze_Index']
})
.setChartType('LineChart')
.setOptions({
  title: 'Landsat Haze Index Proxy (Lowest is Best)',
  hAxis: {title: 'Date', format: 'MMM yyyy'},
  vAxis: {title: 'Mean Blue Reflectance'},
  colors: ['#4575b4'],
  lineWidth: 2,
  pointSize: 3
});

// Print directly to the console
print('Time-Series: Optimal Window Clarity', cloudChart);


// =========================================================================
// UI APP & INTERACTIVITY
// =========================================================================

// Create a main panel to hold the UI
var panel = ui.Panel({
  style: {width: '350px', padding: '15px'}
});
ui.root.insert(0, panel);

// UI Elements
var title = ui.Label('Landsat 8-Day Composite Viewer', {fontWeight: 'bold', fontSize: '20px'});
var desc = ui.Label('Select a Year and Month. The script pulls the clearest 8-day composite directly from your pre-calculated asset.');

var yearLabel = ui.Label('Select Year:', {fontWeight: 'bold'});
var yearSlider = ui.Slider({
  min: 1984, max: 2025, value: 2019, step: 1, // <-- FIXED: Now begins in 1984
  style: {stretch: 'horizontal'}
});

var monthLabel = ui.Label('Select Month:', {fontWeight: 'bold'});
var monthSlider = ui.Slider({
  min: 1, max: 12, value: 5, step: 1,
  style: {stretch: 'horizontal'}
});

var statusBox = ui.Label({
  value: 'Ready. Move sliders to load imagery...',
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
Map.setOptions('HYBRID'); 

// The main function that runs when sliders change
function updateMap() {
  var y = yearSlider.getValue();
  var m = monthSlider.getValue();
  
  statusBox.setValue('Fetching composite from asset...\n(This is nearly instant!)');
  statusBox.style().set('color', 'orange');
  
  // Query the asset for this specific year and month
  var window_query = cloud_windows
    .filter(ee.Filter.eq('Year', y))
    .filter(ee.Filter.eq('Month', m));
  
  // Check if a valid window exists
  window_query.size().evaluate(function(count) {
    if (count === 0) {
      statusBox.setValue('No data found in asset for ' + y + '-' + m + '.\n(Try a different date)');
      statusBox.style().set('color', 'red');
      Map.layers().reset(); 
      Map.addLayer(bounds_fc.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'SRER Bounds');
      return;
    }
    
    var optimal_window = ee.Feature(window_query.first());
    
    optimal_window.evaluate(function(feature) {
      var props = feature.properties;
      
      // Handle the "No Data" dummy features we created in the asset script
      if (props.Window_Label.indexOf('No Data') !== -1) {
        statusBox.setValue('Asset indicates No Data available for ' + y + '-' + m);
        statusBox.style().set('color', 'red');
        Map.layers().reset(); 
        Map.addLayer(bounds_fc.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'SRER Bounds');
        return;
      }
      
      // Update the UI Text
      var statusText = 'SUCCESS!' + 
                       '\nOptimal Dates: ' + props.Window_Label + 
                       '\nHaze Index (Blue): ' + props.Mean_Haze_Index.toFixed(3);
      statusBox.setValue(statusText);
      statusBox.style().set('color', 'green');
      
      // Reconstruct the image collection for the optimal dates
      var startDate = ee.Date(props.Start_Date);
      var endDate = startDate.advance(8, 'day'); // Landsat 8-day composite window
      
      var landsat_img = ee.ImageCollection('LANDSAT/COMPOSITES/C02/T1_L2_8DAY')
        .filterBounds(bounds_geom)
        .filterDate(startDate, endDate)
        .mosaic() 
        .clip(v_extent);
        
      // Update the Map Layers
      Map.layers().reset(); 
      
      // Layer 0: The true color composite (Using the image directly, no scaling math needed!)
      Map.layers().set(0, ui.Map.Layer(landsat_img, rgbVis, 'Landsat ' + props.Window_Label));
      
      // Layer 1: An empty red outline for your bounds
      Map.layers().set(1, ui.Map.Layer(bounds_fc.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'SRER Bounds'));
    });
  });
}

// Attach the update function to the sliders
yearSlider.onChange(updateMap);
monthSlider.onChange(updateMap);

// Run it once on script start to load the default (May 2019)
updateMap();
