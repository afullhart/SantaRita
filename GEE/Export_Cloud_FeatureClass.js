// =========================================================================
// SETUP & ASSETS
// =========================================================================
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');

// Extract the raw Geometry from the first feature directly
var bounds_geom = bounds_fc.first().geometry().bounds();

// =========================================================================
// SPECIFIC YEAR-MONTH 8-DAY COMPOSITE ANALYSIS
// =========================================================================

// Define the historical period to analyze
var start_year = 1984; 
var end_year = 2025;   

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
    
    // Calculate exact start and end dates for the month
    var startDate = ee.Date.fromYMD(y, m, 1);
    var endDate = startDate.advance(1, 'month'); 
    
    // Filter the Landsat 8-Day Composite record to this exact month
    var landsat_month_col = ee.ImageCollection('LANDSAT/COMPOSITES/C02/T1_L2_8DAY')
      .filterBounds(bounds_geom)
      .filterDate(startDate, endDate);
      
    var imgCount = landsat_month_col.size();
    
    // Use ee.Algorithms.If to gracefully handle months with absolutely no imagery
    return ee.Algorithms.If(
      imgCount.eq(0),
      
      // If empty, return a dummy feature with an artificially high haze value
      ee.Feature(bounds_geom, {
        'Year': y,
        'Month': m,
        'Start_Date': startDate.format('YYYY-MM-dd'),
        'Window_Label': ee.String(y.format('%d')).cat('-').cat(m.format('%02d')).cat(' No Data'),
        'Mean_Haze_Index': 99999, 
        'Image_Count': 0,
        'system:time_start': startDate.millis() 
      }),
      
      // If images exist, map over them to calculate the Blue Band Haze Proxy
      landsat_month_col.map(function(img) {
        
        // Calculate the mean value of the 'blue' band across the ranch
        // High blue = Clouds/Haze. Low blue = Clear skies.
        var meanBlueRaw = img.select('blue').reduceRegion({
          reducer: ee.Reducer.mean(),
          geometry: bounds_geom,
          scale: 30, // Landsat native resolution
          maxPixels: 1e9
        }).get('blue');
        
        // Failsafe in case reduceRegion returns null
        var meanBlue = ee.Algorithms.If(
          ee.Algorithms.IsEqual(meanBlueRaw, null), 
          99999, 
          meanBlueRaw
        );
        
        // Format a clean label for charting/export
        var imgDate = ee.Date(img.get('system:time_start'));
        var windowLabel = ee.String(y.format('%d')).cat('-')
                  .cat(m.format('%02d')).cat('-')
                  .cat(imgDate.format('dd'));
        
        return ee.Feature(bounds_geom, {
          'Year': y,
          'Month': m,
          'Start_Date': imgDate.format('YYYY-MM-dd'),
          'Window_Label': windowLabel,
          'Mean_Haze_Index': meanBlue, 
          'Image_Count': 1,
          'system:time_start': img.get('system:time_start') 
        });
      })
      // Sort the month's available composites by the lowest Haze Index (clearest day)
      .sort('Mean_Haze_Index')
      .first()
    );
  });
  
  return bestForYear;
});

// Flatten the List of Lists of Features into a single 1D FeatureCollection
var bestWindowsFc = ee.FeatureCollection(allBestWindowsList.flatten())
  .filter(ee.Filter.neq('Image_Count', 0)) // Drop the "No Data" dummy features
  .sort('system:time_start'); // Ensure chronological order


// =========================================================================
// OUTPUTS & EXPORTS
// =========================================================================

// Print the Master Feature Collection to the console
print('All Best Monthly Landsat Composites (Export Data):', bestWindowsFc);

// Loop through our client-side list of years and generate ONE chart per year
clientYears.forEach(function(year) {
  
  // Filter the master collection to just this specific year
  var yearly_fc = bestWindowsFc.filter(ee.Filter.eq('Year', year));
  
  var chart = ui.Chart.feature.byFeature({
    features: yearly_fc,
    xProperty: 'Window_Label',
    yProperties: ['Mean_Haze_Index']
  })
  .setChartType('ColumnChart')
  .setOptions({
    title: 'Landsat Clarity Proxy (Lowest is Best) for ' + year,
    hAxis: {
      title: '8-Day Composite Start Date', 
      slantedText: true, 
      slantedTextAngle: 45
    },
    vAxis: {
      title: 'Mean Blue Reflectance (DN)'
    },
    colors: ['#4575b4'],
    legend: {position: 'none'}
  });
  
  print(chart);
});

// Export the comprehensive master table to Google Drive
Export.table.toDrive({
  collection: bestWindowsFc,
  description: 'SRER_Best_Landsat_Monthly_Composites',
  folder: 'GEE_Downloads',
  fileFormat: 'CSV'
});

// Export the Feature Collection directly to an Earth Engine Asset
Export.table.toAsset({
  collection: bestWindowsFc,
  description: 'Export_Landsat_Cloud_FeatureClass_Asset',
  assetId: 'projects/ee-andrewfullhart/assets/Cloud_FeatureClass_Landsat'
});
