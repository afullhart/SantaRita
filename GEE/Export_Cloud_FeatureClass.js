// =========================================================================
// SETUP & ASSETS
// =========================================================================
var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm');
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');

// Extract the raw Geometry from the first feature directly
var bounds_geom = bounds_fc.first().geometry().bounds();

// =========================================================================
// MONTHLY LANDSAT INVENTORY (1984 - 2025)
// =========================================================================
print('Scanning 40 years of Landsat 5, 7, 8, and 9... This will take a moment.');

var start_year = 1984; 
var end_year = 2025;   

var clientYears = [];
for (var y = start_year; y <= end_year; y++) {
  clientYears.push(y);
}

var years = ee.List(clientYears);
var months = ee.List.sequence(1, 12);

// Map over every Year
var allMonthlyDataList = years.map(function(y) {
  y = ee.Number(y);
  
  // Map over every Month within that Year
  var monthlyData = months.map(function(m) {
    m = ee.Number(m);
    
    var startDate = ee.Date.fromYMD(y, m, 1);
    var endDate = startDate.advance(1, 'month'); 
    
    // 1. Check Landsat 9
    var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
      .filterBounds(bounds_geom).filterDate(startDate, endDate);
      
    // 2. Check Landsat 8
    var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
      .filterBounds(bounds_geom).filterDate(startDate, endDate);
      
    // 3. Check Landsat 7
    var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
      .filterBounds(bounds_geom).filterDate(startDate, endDate);
      
    // 4. Check Landsat 5
    var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
      .filterBounds(bounds_geom).filterDate(startDate, endDate);
      
    // Merge them all together to get the total count of images for this month
    var combined_month_col = l9.merge(l8).merge(l7).merge(l5);
    var imgCount = combined_month_col.size();
    
    // Format the label
    var windowLabel = ee.Algorithms.If(
      imgCount.eq(0),
      ee.String(y.format('%d')).cat('-').cat(m.format('%02d')).cat(' (No Data)'),
      ee.String(y.format('%d')).cat('-').cat(m.format('%02d')).cat(' (Median of ').cat(imgCount.format('%d')).cat(' imgs)')
    );

    return ee.Feature(bounds_geom, {
      'Year': y,
      'Month': m,
      'Start_Date': startDate.format('YYYY-MM-dd'),
      'Window_Label': windowLabel,
      'Image_Count': imgCount,
      'system:time_start': startDate.millis() 
    });
  });
  
  return monthlyData;
});

// Flatten the List of Lists of Features into a single 1D FeatureCollection
var monthlyIndexFc = ee.FeatureCollection(allMonthlyDataList.flatten())
  .filter(ee.Filter.neq('Image_Count', 0)) // Drop the empty months
  .sort('system:time_start'); // Ensure chronological order


// =========================================================================
// OUTPUTS & EXPORTS
// =========================================================================

// Print the Master Feature Collection to the console
print('Monthly Landsat Image Index (Export Data):', monthlyIndexFc);

// Generate a Time-Series Chart showing Data Density over time
var chart = ui.Chart.feature.byFeature({
  features: monthlyIndexFc,
  xProperty: 'system:time_start',
  yProperties: ['Image_Count']
})
.setChartType('ColumnChart')
.setOptions({
  title: 'Landsat Monthly Image Density (1984 - 2025)',
  hAxis: {title: 'Date', format: 'yyyy'},
  vAxis: {title: 'Number of Images per Month'},
  colors: ['#1a9850'],
  legend: {position: 'none'}
});

print('Data Density Chart:', chart);

// Export the comprehensive master table to Google Drive
Export.table.toDrive({
  collection: monthlyIndexFc,
  description: 'SRER_Landsat_Monthly_Image_Index',
  folder: 'GEE_Downloads',
  fileFormat: 'CSV'
});

// Export the Feature Collection directly to your Earth Engine Asset
Export.table.toAsset({
  collection: monthlyIndexFc,
  description: 'Export_Landsat_Monthly_Index_Asset',
  assetId: 'projects/ee-andrewfullhart/assets/Cloud_FeatureClass_Landsat'
});
