// =========================================================================
// SETUP & ASSETS
// =========================================================================
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_geom = bounds_fc.first().geometry().bounds();

// =========================================================================
// CLOUD MASKING LOGIC
// =========================================================================
function maskClouds(qa) {
  var dilatedCloud = qa.bitwiseAnd(1 << 1).eq(0);
  var cirrus = qa.bitwiseAnd(1 << 2).eq(0);
  var cloud = qa.bitwiseAnd(1 << 3).eq(0);
  var shadow = qa.bitwiseAnd(1 << 4).eq(0);
  return dilatedCloud.and(cirrus).and(cloud).and(shadow);
}

// =========================================================================
// MONTHLY LANDSAT INVENTORY (1984 - 2025)
// =========================================================================
print('Scanning and filtering 40 years of Landsat imagery...');

var start_year = 1984; 
var end_year = 2025;   
var years = ee.List.sequence(start_year, end_year);
var months = ee.List.sequence(1, 12);

// Map over every Year
var allMonthlyDataList = years.map(function(y) {
  y = ee.Number(y);
  
  // Map over every Month
  var monthlyData = months.map(function(m) {
    m = ee.Number(m);
    
    var startDate = ee.Date.fromYMD(y, m, 1);
    var endDate = startDate.advance(1, 'month'); 
    
    var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(bounds_geom).filterDate(startDate, endDate);
    var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(bounds_geom).filterDate(startDate, endDate);
    var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(bounds_geom).filterDate(startDate, endDate);
    var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(bounds_geom).filterDate(startDate, endDate);
      
    var combined_month_col = l9.merge(l8).merge(l7).merge(l5);
    
    // --- THE FIX: CALCULATE CLEAR PIXEL FRACTION ---
    var scored_col = combined_month_col.map(function(img) {
      var qa = img.select('QA_PIXEL');
      var clear_mask = maskClouds(qa);
      var native_mask = img.select('SR_B2').mask(); // Identifies data gaps (like L7 stripes)
      
      // Combine masks. unmask(0) forces missing swaths to pull down the average.
      var master_mask = clear_mask.and(native_mask).unmask(0).rename('Quality');
      
      // Calculate percentage of the bounds that are perfectly clear and valid
      var fraction = master_mask.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: bounds_geom,
        scale: 30,
        maxPixels: 1e9
      }).get('Quality');
      
      fraction = ee.Algorithms.If(ee.Algorithms.IsEqual(fraction, null), 0, fraction);
      return img.set('Local_Clear_Fraction', fraction);
    });
    
    // Filter out garbage images (Must cover >70% of the bounds clearly)
    var good_col = scored_col.filter(ee.Filter.gte('Local_Clear_Fraction', 0.70));
    var imgCount = good_col.size();
    
    // Extract the exact system IDs of the surviving images
    var valid_ids = good_col.aggregate_array('system:id').join(',');
    
    var windowLabel = ee.Algorithms.If(
      imgCount.eq(0),
      ee.String(y.format('%d')).cat('-').cat(m.format('%02d')).cat(' (No Data)'),
      ee.String(y.format('%d')).cat('-').cat(m.format('%02d')).cat(' (Median of ').cat(imgCount.format('%d')).cat(' clean imgs)')
    );

    return ee.Feature(bounds_geom, {
      'Year': y,
      'Month': m,
      'Start_Date': startDate.format('YYYY-MM-dd'),
      'Window_Label': windowLabel,
      'Image_Count': imgCount,
      'Valid_IDs': valid_ids, // Saved for the UI script to read
      'system:time_start': startDate.millis() 
    });
  });
  
  return monthlyData;
});

var monthlyIndexFc = ee.FeatureCollection(allMonthlyDataList.flatten())
  .filter(ee.Filter.neq('Image_Count', 0)) 
  .sort('system:time_start'); 

Export.table.toAsset({
  collection: monthlyIndexFc,
  description: 'Export_Landsat_Monthly_Index_Asset',
  assetId: 'projects/ee-andrewfullhart/assets/Cloud_FeatureClass_Landsat'
});
