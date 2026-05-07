// =========================================================================
// SETUP & ASSETS
// =========================================================================
var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm');
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');

// Extract the raw Geometry from the first feature directly
var bounds_geom = bounds_fc.first().geometry();


// =========================================================================
// SPECIFIC YEAR-MONTH 7-DAY WINDOW ANALYSIS
// =========================================================================

// Define the historical period to analyze
var start_year = 2018; // Sentinel-2 SR data becomes robust globally around 2018
var end_year = 2023; 

// Build a client-side array of years so we can loop over them later to generate individual charts
var clientYears = [];
for (var y = start_year; y <= end_year; y++) {
  clientYears.push(y);
}

// Convert client-side array to an Earth Engine List for server-side mapping
var years = ee.List(clientYears);
var months = ee.List.sequence(1, 12);

// Map over every Year
var allBestWindowsList = years.map(function(y) {
  y = ee.Number(y);
  
  // Map over every Month within that Year
  var bestForYear = months.map(function(m) {
    m = ee.Number(m);
    
    // Determine the exact number of days in this specific month/year (handles leap years)
    var refDate = ee.Date.fromYMD(y, m, 1);
    var daysInMonth = refDate.advance(1, 'month').difference(refDate, 'day');
    
    // Create a rolling 7-day window. Last valid start day is (daysInMonth - 6)
    var startDays = ee.List.sequence(1, daysInMonth.subtract(6));
    
    var windows = startDays.map(function(d) {
      var startDay = ee.Number(d);
      
      // Calculate exact start and end dates for the filter
      var startDate = ee.Date.fromYMD(y, m, startDay);
      var endDateFilter = startDate.advance(7, 'day'); // filterDate is exclusive of the end date
      var endDateInclusive = startDate.advance(6, 'day'); // For the text label
      
      // Filter the Sentinel-2 record to this exact 7-day period
      var s2_window = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
        .filterBounds(bounds_geom)
        .filterDate(startDate, endDateFilter);
        
      var imgCount = s2_window.size();
        
      // Calculate the mean cloud probability
      var meanCldImg = s2_window.select('MSK_CLDPRB').mean();
      
      var meanCld = ee.Algorithms.If(
        imgCount.eq(0),
        100, // Penalize windows with absolutely no imagery
        meanCldImg.reduceRegion({
          reducer: ee.Reducer.mean(),
          geometry: bounds_geom,
          scale: 60,
          maxPixels: 1e9
        }).get('MSK_CLDPRB')
      );
      
      // Failsafe in case reduceRegion returns null
      meanCld = ee.Algorithms.If(ee.Algorithms.IsEqual(meanCld, null), 100, meanCld);
      
      // Format a clean label for charting/export
      var windowLabel = ee.String(y.format('%d')).cat('-')
                  .cat(m.format('%02d')).cat('-')
                  .cat(startDay.format('%02d')).cat(' to ')
                  .cat(endDateInclusive.format('%02d'));
      
      // We must calculate a timestamp to sort the final collection chronologically
      var timeStart = startDate.millis();
      
      return ee.Feature(null, {
        'Year': y,
        'Month': m,
        'Start_Date': startDate.format('YYYY-MM-dd'),
        'End_Date': endDateInclusive.format('YYYY-MM-dd'),
        'Window_Label': windowLabel,
        'Mean_Cloud_Prob': meanCld,
        'Image_Count': imgCount,
        'system:time_start': timeStart 
      });
    });
    
    var windowsFc = ee.FeatureCollection(windows);
    
    // Sort all windows in THIS specific month, and grab the lowest cloud probability
    return windowsFc.sort('Mean_Cloud_Prob').first();
  });
  
  return bestForYear;
});

// Flatten the List of Lists of Features into a single 1D FeatureCollection
var bestWindowsFc = ee.FeatureCollection(allBestWindowsList.flatten())
  .sort('system:time_start'); // Ensure chronological order


// =========================================================================
// OUTPUTS & EXPORTS
// =========================================================================

// Print the Master Feature Collection to the console
print('All Best 7-Day Windows (Export Data):', bestWindowsFc);

// Loop through our client-side list of years and generate ONE chart per year
clientYears.forEach(function(year) {
  
  // Filter the master collection to just this specific year
  var yearly_fc = bestWindowsFc.filter(ee.Filter.eq('Year', year));
  
  var chart = ui.Chart.feature.byFeature({
    features: yearly_fc,
    xProperty: 'Window_Label',
    yProperties: ['Mean_Cloud_Prob']
  })
  .setChartType('ColumnChart')
  .setOptions({
    title: 'Optimal Cloud Probability Windows for ' + year,
    hAxis: {
      title: '7-Day Window', 
      slantedText: true, 
      slantedTextAngle: 45
    },
    vAxis: {
      title: 'Mean Cloud Probability (%)',
      viewWindow: {min: 0, max: 100} // Locks the Y-axis so charts are easily comparable
    },
    colors: ['#1a9850'],
    legend: {position: 'none'}
  });
  
  print(chart);
});

// Export the comprehensive master table to Google Drive
Export.table.toDrive({
  collection: bestWindowsFc,
  description: 'SRER_Best_7Day_Imagery_Windows_By_Year',
  folder: 'GEE_Downloads',
  fileFormat: 'CSV'
});

