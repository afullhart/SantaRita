// Takes ~1 hr for the feature collection export to asset

// ========================================================================= 
// SETUP & ASSETS 
// ========================================================================= 
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
// Extract the raw Geometry from the first feature directly 
var bounds_geom = bounds_fc.first().geometry(); 
var v_classified_may = ee.Image('users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_may_2019_5cm');
var v_classified_sep = ee.Image('users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_sep_2019_5cm');
var v_foot_prints = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_drone_footprints');

// Import the ecological states polygons again 
var v_srer_polys = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_ecological_states')
  .map(function(ft){
    return ft.set('area_ha', ft.area(1).divide(10000));
  });

// Safely returns an ee.Geometry bounding box 
var v_extent = bounds_geom.bounds(); 

// ========================================================================= 
// TOPOGRAPHIC PRE-PROCESSING 
// ========================================================================= 
var dem = ee.Image('USGS/3DEP/10m').clip(v_extent);
var terrain_slope = ee.Terrain.slope(dem).multiply(Math.PI / 180); 
var terrain_aspect = ee.Terrain.aspect(dem).multiply(Math.PI / 180);

// ========================================================================= 
// DYNAMIC DATES (PULLED FROM LANDSAT CLOUD ASSET) 
// ========================================================================= 
var cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass_Landsat');

var may_window = ee.Feature(cloud_windows.filter(ee.Filter.eq('Year', 2019)).filter(ee.Filter.eq('Month', 5)).first());
var may_start = ee.String(may_window.get('Start_Date')); 
var may_end   = ee.Date(may_start).advance(1, 'month'); 

var sep_window = ee.Feature(cloud_windows.filter(ee.Filter.eq('Year', 2019)).filter(ee.Filter.eq('Month', 9)).first());
var sep_start = ee.String(sep_window.get('Start_Date')); 
var sep_end   = ee.Date(sep_start).advance(1, 'month');

// ========================================================================= 
// PART 1: STATIC GRID GENERATION (LANDSAT 30m ALIGNMENT) 
// ========================================================================= 
var projLandsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(v_extent)
  .first()
  .select('SR_B2')
  .projection();

var landsat_grid = v_extent.coveringGrid(projLandsat, projLandsat.nominalScale());
var focused_grid = landsat_grid.filterBounds(v_foot_prints);
var valid_area_mask = ee.Image.constant(0).paint(v_foot_prints, 1);

var grid_overlap = valid_area_mask.reduceRegions({
  collection: focused_grid,
  reducer: ee.Reducer.mean(),
  scale: 1,
  crs: projLandsat.crs(),
  tileScale: 4 
});

var final_grid = grid_overlap.filter(ee.Filter.gte('mean', 0.99));

var v_spatial_filter = ee.Filter.intersects({leftField:'.geo', rightField:'.geo', maxError:1});
var v_saveFirstJoin = ee.Join.saveFirst({matchKey: 'poly', outer: true}); 
var v_joined_grids = v_saveFirstJoin.apply(final_grid, v_srer_polys, v_spatial_filter)
  .map(function(ft) {
    var hasPoly = ft.propertyNames().contains('poly');
    var polyFeat = ee.Feature(ee.Algorithms.If(
      hasPoly,
      ft.get('poly'),
      ee.Feature(null, {'Plant_Comm': '','Pasture': '','Transect': '','Utility': '','S_Desc': '','Exclosure': '','area_ha': null})
    ));

    return ft.set({
      'Plant_Comm': polyFeat.get('Plant_Comm'),
      'Pasture': polyFeat.get('Pasture'),
      'Transect': polyFeat.get('Transect'),
      'Utility': polyFeat.get('Utility'),
      'S_Desc': polyFeat.get('S_Desc'),
      'Exclosure': polyFeat.get('Exclosure'),
      'area_ha': polyFeat.get('area_ha'),
      'poly': null 
    });
  });

// ========================================================================= 
// PART 2: EXTRACT LANDSAT & TOPOGRAPHIC PREDICTORS 
// =========================================================================
function maskClouds(qa) {
  var dilatedCloud = qa.bitwiseAnd(1 << 1).eq(0);
  var cirrus = qa.bitwiseAnd(1 << 2).eq(0);
  var cloud = qa.bitwiseAnd(1 << 3).eq(0);
  var shadow = qa.bitwiseAnd(1 << 4).eq(0);
  return dilatedCloud.and(cirrus).and(cloud).and(shadow);
}

function prepOLI(img) {
  var scaled = img.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'])
                  .multiply(0.0000275).add(-0.2); 
                  
  var qa = img.select('QA_PIXEL');
  var qaMask = maskClouds(qa);
  
  var darkMask = scaled.select('SR_B5').gt(0.15); 
  var baseMask = qaMask.and(darkMask);
  var finalMask = baseMask.focal_min({radius: 30, kernelType: 'square', units: 'meters'});

  var sz_num = ee.Number(img.get('SUN_ELEVATION')).subtract(90).multiply(-1).multiply(Math.PI / 180);
  var sa_num = ee.Number(img.get('SUN_AZIMUTH')).multiply(Math.PI / 180);
  
  var cosZ = ee.Image.constant(sz_num.cos());
  var sinZ = ee.Image.constant(sz_num.sin());
  var sa_img = ee.Image.constant(sa_num);
  
  var cosS = terrain_slope.cos();
  var sinS = terrain_slope.sin();
  var cosAzAsp = sa_img.subtract(terrain_aspect).cos();
  var illumination = cosZ.multiply(cosS).add(sinZ.multiply(sinS).multiply(cosAzAsp)).rename('illumination');

  var optical = scaled.rename(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']).addBands(illumination);
  return optical.updateMask(finalMask);
}

function extractLandsatData(startDate, endDate, monthLabel) {
  var ls_ic = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(v_extent)
    .filterDate(startDate, endDate)
    .map(prepOLI);

  var ls_median = ls_ic.median().clip(v_extent); 
  
  // =====================================================================
  // COMPLETELY REPLACE THE IMAGE WITH A SYNTHETIC MEDIAN
  // =====================================================================
  if (monthLabel === 'Sept') {
    var aug31_col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
      .filterBounds(v_extent)
      .filterDate('2019-08-31', '2019-09-01')
      .map(prepOLI); 
      
    var oct02_col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
      .filterBounds(v_extent)
      .filterDate('2019-10-02', '2019-10-03')
      .map(prepOLI);

    ls_median = aug31_col.merge(oct02_col).median().clip(v_extent);
  }
  // =====================================================================

  var optical_bands = ls_median.setDefaultProjection({crs: projLandsat.crs(), scale: projLandsat.nominalScale()});

  var ndvi = optical_bands.normalizedDifference(['NIR', 'Red']).rename('NDVI');
  var bsi = optical_bands.expression(
    '((SWIR1 + Red) - (NIR + Blue)) / ((SWIR1 + Red) + (NIR + Blue))', {
      'Blue':  optical_bands.select('Blue'),
      'Red':   optical_bands.select('Red'),
      'NIR':   optical_bands.select('NIR'),
      'SWIR1': optical_bands.select('SWIR1')
    }).rename('BSI');

  var nbr2 = optical_bands.normalizedDifference(['SWIR1', 'SWIR2']).rename('NBR2');

  var final_img = optical_bands.addBands([
    ndvi, bsi, nbr2,
    terrain_slope.rename('slope'),
    terrain_aspect.rename('aspect')
  ]);

  var bandsToExtract = final_img.select([
    'Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2',
    'NDVI', 'BSI', 'NBR2', 'slope', 'illumination', 'aspect'
  ]);

  var gridWithMonth = v_joined_grids.map(function(feat) {
    return feat.set('Month', monthLabel);
  });

  var extracted = bandsToExtract.reduceRegions({
    collection: gridWithMonth,
    reducer: ee.Reducer.mean(),
    scale: projLandsat.nominalScale(),
    crs: projLandsat.crs(),
    tileScale: 4
  });

  return extracted;
} 

var may_ls_data = extractLandsatData(may_start, may_end, 'May'); 
var sep_ls_data = extractLandsatData(sep_start, sep_end, 'Sept'); 

var may_ls_utm = may_ls_data.map(function(f) { return f.transform('EPSG:32612', 0.05); });
var sep_ls_utm = sep_ls_data.map(function(f) { return f.transform('EPSG:32612', 0.05); });

// ========================================================================= 
// PART 3: DRONE METRICS (BGR, LPI, MFT, HERB/WOODY) 
// =========================================================================
function processMonthMetrics(grid_subset, classified_img) {
  var native_proj = classified_img.projection();
  var binary = classified_img.eq(3).selfMask(); 
  
  var herbMask = classified_img.eq(1).rename('herb');
  var woodyMask = classified_img.eq(2).rename('woody');
  
  // Multiply by pixelArea to calculate raw square meter coverage
  var hwCombined = herbMask.addBands(woodyMask).multiply(ee.Image.pixelArea());
  
  return grid_subset.map(function(ft){
    
    // --- 1. BGR Calculation --- 
    var v_area_image = binary.multiply(ee.Image.pixelArea());
    var v_area_ft = ft.area(0.05);
    var v_area = v_area_image.reduceRegion({
      reducer: ee.Reducer.sum(), geometry: ft.geometry(), scale: 0.05, maxPixels: 1e13
    }).get('classification');
    var v_pct_area = ee.Number(v_area).divide(v_area_ft).multiply(100);

    // --- 2. LPI Calculation --- 
    var patch_vectors = binary.reduceToVectors({
      reducer: ee.Reducer.countEvery(), geometry: ft.geometry(), scale: 0.05, crs: native_proj.crs(),
      geometryType: 'polygon', eightConnected: true, labelProperty: 'class_val', maxPixels: 1e13, tileScale: 16
    });

    var patches_with_area = patch_vectors.map(function(feat){
      return feat.set('area_sqm', feat.area(0.05));
    });
    
    var areas = patches_with_area.aggregate_array('area_sqm');
    var max_area = ee.Algorithms.If(
      areas.length().gt(0), areas.reduce(ee.Reducer.max()), 0
    );

    // Convert raw sq meters into a true percentage index
    var lpi_pct = ee.Number(max_area).divide(v_area_ft).multiply(100);

    // --- 3. MFT Calculation --- 
    function Get_Mean_Fetch(c_image, geom) {
      var solid_img = c_image.unmask(99);
      var obstacles = solid_img.neq(3); 

      var v_distance = obstacles.fastDistanceTransform(2048).sqrt()
        .multiply(0.05).rename("distance");

      var N_PTS = 1000;
      var v_rnd = ee.FeatureCollection.randomPoints({
        region: geom, points: N_PTS, seed: 1234
      });

      var sampled_points = v_distance.reduceRegions({
        collection: v_rnd, reducer: ee.Reducer.first().setOutputs(["distance"]), scale: 0.05
      });

      var valid_points = sampled_points.filter(ee.Filter.notNull(['distance']));
      return ee.Number(valid_points.aggregate_mean('distance'));
    } 
    
    var v_mft = Get_Mean_Fetch(classified_img, ft.geometry());

    // --- 4. HERB, WOODY, & RATIO CALCULATION ---
    var hwAreas = hwCombined.reduceRegion({
      reducer: ee.Reducer.sum(), geometry: ft.geometry(), scale: 0.05, maxPixels: 1e13
    });

    var herbArea = ee.Number(hwAreas.get('herb'));
    var woodyArea = ee.Number(hwAreas.get('woody'));

    // Convert raw square meters into a percentage of the total cell
    var herbPct = herbArea.divide(v_area_ft).multiply(100);
    var woodyPct = woodyArea.divide(v_area_ft).multiply(100);

    var hwRatio = ee.Algorithms.If(
      woodyArea.gt(0), herbArea.divide(woodyArea), null 
    );

    return ft.set(
      'LPI', lpi_pct, 
      'BGR', v_pct_area, 
      'MFT', v_mft, 
      'Herb_Woody_Ratio', hwRatio,
      'Herb_pct', herbPct,
      'Woody_pct', woodyPct
    );
  }); 
}

var final_may = processMonthMetrics(may_ls_utm, v_classified_may); 
var final_sep = processMonthMetrics(sep_ls_utm, v_classified_sep); 
var final_combined = final_may.merge(final_sep);

// ========================================================================= 
// EXPORTS 
// =========================================================================
Export.table.toDrive({
  collection: final_combined,
  description: 'SRER_Training_FeatureClass_Landsat',
  folder: 'GEE_Downloads',
  fileFormat: 'CSV' 
});

Export.table.toAsset({
  collection: final_combined,
  assetId: 'projects/ee-andrewfullhart/assets/SR_landsat_model_grid_utm',
  description: 'FC_asset' 
});
