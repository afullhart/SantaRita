//Takes ~1 hr for the feature collection export to asset (800000.0000 EECU-seconds)

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
// DYNAMIC DATES (PULLED FROM CLOUD ASSET) 
// ========================================================================= 
var cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass_S2');

// Query the asset for May 2019 
var may_window = ee.Feature(cloud_windows.filter(ee.Filter.eq('Year', 2019)).filter(ee.Filter.eq('Month', 5)).first());
var may_start = ee.String(may_window.get('Start_Date')); 
var may_end   = ee.Date(may_start).advance(7, 'day'); // 7-day exclusive filter limit

// Query the asset for Sept 2019 
var sep_window = ee.Feature(cloud_windows.filter(ee.Filter.eq('Year', 2019)).filter(ee.Filter.eq('Month', 9)).first());
var sep_start = ee.String(sep_window.get('Start_Date')); 
var sep_end   = ee.Date(sep_start).advance(7, 'day');

// ========================================================================= 
// PART 1: STATIC GRID GENERATION (FOOTPRINTS ONLY) 
// ========================================================================= 
// Extract projection from a single image in the collection to build the grid 
var projSent2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(v_extent)
  .first()
  .select('B2')
  .projection();

// Generate the base Sentinel-2 pixel grid 
var sent2_grid = v_extent.coveringGrid(projSent2, projSent2.nominalScale());

// Pre-filter the grid bounds strictly by drone footprints 
var focused_grid = sent2_grid.filterBounds(v_foot_prints);

// Create a high-resolution binary mask (1) ONLY where the footprints overlap 
var valid_area_mask = ee.Image.constant(0).paint(v_foot_prints, 1);

// Calculate the exact percentage of overlap for every grid cell 
var grid_overlap = valid_area_mask.reduceRegions({
  collection: focused_grid,
  reducer: ee.Reducer.mean(),
  scale: 1,
  crs: projSent2.crs(),
  tileScale: 4 
});

// Filter for strict 99% overlap with drone footprints (Syncs to exact ArcPy geometry count)
var final_grid = grid_overlap.filter(ee.Filter.gte('mean', 0.99));

// --- OUTER SPATIAL JOIN ---
// Keep all footprints, but fetch polygon attributes if they overlap 
var v_spatial_filter = ee.Filter.intersects({leftField:'.geo', rightField:'.geo', maxError:1});

var v_saveFirstJoin = ee.Join.saveFirst({matchKey: 'poly', outer: true}); 
var v_sent2_joined_grids = v_saveFirstJoin.apply(final_grid, v_srer_polys, v_spatial_filter)
  .map(function(ft) {
    // Check if the join successfully attached a polygon
    var hasPoly = ft.propertyNames().contains('poly');
    
    // Safely extract the polygon if it exists, otherwise provide blank defaults 
    var polyFeat = ee.Feature(ee.Algorithms.If(
      hasPoly,
      ft.get('poly'),
      ee.Feature(null, {
        'Plant_Comm': '',
        'Pasture': '',
        'Transect': '',
        'Utility': '',
        'S_Desc': '',
        'Exclosure': '',
        'area_ha': null
      })
    ));

    return ft.set({
      'Plant_Comm': polyFeat.get('Plant_Comm'),
      'Pasture': polyFeat.get('Pasture'),
      'Transect': polyFeat.get('Transect'),
      'Utility': polyFeat.get('Utility'),
      'S_Desc': polyFeat.get('S_Desc'),
      'Exclosure': polyFeat.get('Exclosure'),
      'area_ha': polyFeat.get('area_ha'),
      'poly': null // Strip the heavy geometry data to keep export clean 
    });
  });

// ========================================================================= 
// PART 2: EXTRACT SENTINEL-2 & TOPOGRAPHIC PREDICTORS 
// =========================================================================
function extractS2Data(startDate, endDate, monthLabel) {
  var sent2_ic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(v_extent)
    .filterDate(startDate, endDate)
    .map(function(img) {
      return img.set('has_aux_bands', img.bandNames().contains('MSK_CLDPRB'));
    })
    .filter(ee.Filter.eq('has_aux_bands', true))
    .map(function(img) {
      var probMask = img.select('MSK_CLDPRB').lt(20);
      var scl = img.select('SCL');
      var sclMask = scl.neq(8).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11)).and(scl.neq(3));
      var blueMask = img.select('B2').lt(2500);
      var masterMask = probMask.and(sclMask).and(blueMask);
      var maskedImg = img.updateMask(masterMask);

      // --- Calculate Illumination (NO DIVISION APPLIED) --- 
      var sz_num = ee.Number(img.get('MEAN_SOLAR_ZENITH_ANGLE')).multiply(Math.PI / 180);
      var sa_num = ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')).multiply(Math.PI / 180);
      
      var cosZ = ee.Image.constant(sz_num.cos());
      var sinZ = ee.Image.constant(sz_num.sin());
      var sa_img = ee.Image.constant(sa_num);
      
      var cosS = terrain_slope.cos();
      var sinS = terrain_slope.sin();
      var cosAzAsp = sa_img.subtract(terrain_aspect).cos();
      
      var illumination = cosZ.multiply(cosS).add(sinZ.multiply(sinS).multiply(cosAzAsp)).rename('illumination');
      
      return maskedImg.addBands(illumination);
    });

  var sent2_median = sent2_ic.median().clip(v_extent); 
  var optical_bands = sent2_median.select(['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12'])
    .multiply(0.0001)
    .setDefaultProjection({crs: projSent2.crs(), scale: projSent2.nominalScale()});

  var ndvi = optical_bands.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var mcari = optical_bands.expression(
    '((B5 - B4) - 0.2 * (B5 - B3)) * (B5 / B4)', {
      'B3': optical_bands.select('B3'),
      'B4': optical_bands.select('B4'),
      'B5': optical_bands.select('B5')
    }).rename('MCARI');

  var bsi = optical_bands.expression(
    '((B11 + B4) - (B8 + B2)) / ((B11 + B4) + (B8 + B2))', {
      'B2':  optical_bands.select('B2'),
      'B4':  optical_bands.select('B4'),
      'B8':  optical_bands.select('B8'),
      'B11': optical_bands.select('B11')
    }).rename('BSI');

  var nbr2 = optical_bands.normalizedDifference(['B11', 'B12']).rename('NBR2');

  var final_img = optical_bands.addBands([
    ndvi, mcari, bsi, nbr2,
    sent2_median.select('illumination'),
    terrain_slope.rename('slope'),
    terrain_aspect.rename('aspect')
  ]);

  var bandsToExtract = final_img.select([
    'B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12',
    'NDVI', 'MCARI', 'BSI', 'NBR2', 'slope', 'illumination', 'aspect'
  ]);

  // Use the Outer Joined grid collection 
  var gridWithMonth = v_sent2_joined_grids.map(function(feat) {
    return feat.set('Month', monthLabel);
  });

  var extracted = bandsToExtract.reduceRegions({
    collection: gridWithMonth,
    reducer: ee.Reducer.mean(),
    scale: projSent2.nominalScale(),
    crs: projSent2.crs(),
    tileScale: 4
  });

  // Return the collection directly to retain null values for 1:1 shapefile alignment
  return extracted;
} 

var may_s2_data = extractS2Data(may_start, may_end, 'May'); 
var sep_s2_data = extractS2Data(sep_start, sep_end, 'Sept'); 
var may_s2_utm = may_s2_data.map(function(f) { return f.transform('EPSG:32612', 0.05); });
var sep_s2_utm = sep_s2_data.map(function(f) { return f.transform('EPSG:32612', 0.05); });

// ========================================================================= 
// PART 3: DRONE METRICS (BGR, LPI, MFT, HERB/WOODY) 
// =========================================================================
function processMonthMetrics(grid_subset, classified_img) {
  var native_proj = classified_img.projection();
  // Masked binary for BGR and LPI (Targeting class 3)
  var binary = classified_img.eq(3).selfMask(); 
  
  // --- NEW: Create masks for Herb (1) and Woody (2) ---
  var herbMask = classified_img.eq(1).rename('herb');
  var woodyMask = classified_img.eq(2).rename('woody');
  var hwCombined = herbMask.addBands(woodyMask);
  
  return grid_subset.map(function(ft){
    
    // --- 1. BGR Calculation --- 
    var v_area_image = binary.multiply(ee.Image.pixelArea());
    var v_area_ft = ft.area(0.05);
    var v_area = v_area_image.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: ft.geometry(),
      scale: 0.05,
      maxPixels: 1e13
    }).get('classification');
    var v_pct_area = ee.Number(v_area).divide(v_area_ft).multiply(100);

    // --- 2. LPI Calculation --- 
    var patch_vectors = binary.reduceToVectors({
      reducer: ee.Reducer.countEvery(),
      geometry: ft.geometry(),
      scale: 0.05,
      crs: native_proj.crs(),
      geometryType: 'polygon',
      eightConnected: true,
      labelProperty: 'class_val',
      maxPixels: 1e13,
      tileScale: 16
    });

    var patches_with_area = patch_vectors.map(function(feat){
      return feat.set('area_sqm', feat.area(0.05));
    });
    
    var areas = patches_with_area.aggregate_array('area_sqm');
    var max_area = ee.Algorithms.If(
      areas.length().gt(0),
      areas.reduce(ee.Reducer.max()),
      0
    );

    // --- 3. MFT Calculation (1000 Random Points + Void Barrier) --- 
    function Get_Mean_Fetch(c_image, geom) {
      // THE VOID FIX: Unmask the image with a dummy value (99).
      var solid_img = c_image.unmask(99);
      var obstacles = solid_img.neq(3); 

      // Calculate pixel distance to nearest obstacle (1), then multiply by 0.05m
      var v_distance = obstacles.fastDistanceTransform(2048).sqrt()
        .multiply(0.05) 
        .rename("distance");

      // Generate 1000 random points strictly within the cell geometry
      var N_PTS = 1000;
      var v_rnd = ee.FeatureCollection.randomPoints({
        region: geom, 
        points: N_PTS, 
        seed: 1234
      });

      // Sample the natively metric distance surface
      var sampled_points = v_distance.reduceRegions({
        collection: v_rnd,
        reducer: ee.Reducer.first().setOutputs(["distance"]),
        scale: 0.05
      });

      // THE MATH FIX: Filter out nulls and use Earth Engine's native mean aggregator
      var valid_points = sampled_points.filter(ee.Filter.notNull(['distance']));
      
      return ee.Number(valid_points.aggregate_mean('distance'));
    } 
    
    var v_mft = Get_Mean_Fetch(classified_img, ft.geometry());

    // --- 4. HERB-TO-WOODY RATIO ---
    // Single reducer pass extracts both herb and woody pixel counts
    var hwCounts = hwCombined.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: ft.geometry(),
      scale: 0.05,
      maxPixels: 1e13
    });

    var herbCount = ee.Number(hwCounts.get('herb'));
    var woodyCount = ee.Number(hwCounts.get('woody'));

    var hwRatio = ee.Algorithms.If(
      woodyCount.gt(0), 
      herbCount.divide(woodyCount), 
      null // Null fallback for division by zero
    );

    return ft.set('LPI', max_area, 'BGR', v_pct_area, 'MFT', v_mft, 'Herb_Woody_Ratio', hwRatio);
  }); 
}

var final_may = processMonthMetrics(may_s2_utm, v_classified_may); 
var final_sep = processMonthMetrics(sep_s2_utm, v_classified_sep); 
var final_combined = final_may.merge(final_sep);

// ========================================================================= 
// EXPORTS 
// =========================================================================
Export.table.toDrive({
  collection: final_combined,
  description: 'SRER_Training_FeatureClass',
  folder: 'GEE_Downloads',
  fileFormat: 'CSV' 
});

Export.table.toAsset({
  collection: final_combined,
  assetId: 'projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm',
  description: 'FC_asset' 
});
