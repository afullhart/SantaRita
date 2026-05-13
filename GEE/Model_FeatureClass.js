// =========================================================================
// SETUP & ASSETS
// =========================================================================

var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');

// Extract the raw Geometry from the first feature directly
var bounds_geom = bounds_fc.first().geometry();

var v_classified_may = ee.Image('users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_may_2019_5cm');
var v_classified_sep = ee.Image('users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_sep_2019_5cm');

var v_foot_prints = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_drone_footprints');

var v_srer_polys = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_ecological_states')
                     .map(function(ft){
                       return ft.set('area_ha', ft.area(1).divide(10000)); 
                     });

// Safely returns an ee.Geometry bounding box
var v_extent = bounds_geom.bounds();

// =========================================================================
// TOPOGRAPHIC PRE-PROCESSING
// =========================================================================
// Calculate terrain properties in radians once for the whole script
var dem = ee.Image('USGS/3DEP/10m').clip(v_extent);
var terrain_slope = ee.Terrain.slope(dem).multiply(Math.PI / 180);
var terrain_aspect = ee.Terrain.aspect(dem).multiply(Math.PI / 180);

// =========================================================================
// DYNAMIC DATES (PULLED FROM CLOUD ASSET)
// =========================================================================
var cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass');

// Query the asset for May 2019
var may_window = ee.Feature(cloud_windows.filter(ee.Filter.eq('Year', 2019)).filter(ee.Filter.eq('Month', 5)).first());
var may_start = ee.String(may_window.get('Start_Date'));
var may_end   = ee.Date(may_start).advance(7, 'day'); // 7-day exclusive filter limit

// Query the asset for Sept 2019
var sep_window = ee.Feature(cloud_windows.filter(ee.Filter.eq('Year', 2019)).filter(ee.Filter.eq('Month', 9)).first());
var sep_start = ee.String(sep_window.get('Start_Date'));
var sep_end   = ee.Date(sep_start).advance(7, 'day');

// =========================================================================
// PART 1: STATIC GRID GENERATION
// =========================================================================

// Extract projection from a single image in the collection to build the grid
var projSent2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(v_extent)
  .first()
  .select('B2')
  .projection();

// Generate the base Sentinel-2 pixel grid
var sent2_grid = v_extent.coveringGrid(projSent2, projSent2.nominalScale());

// Pre-filter the grid so we don't calculate overlap on empty space
var focused_grid = sent2_grid.filterBounds(v_srer_polys).filterBounds(v_foot_prints);

// Create a high-resolution binary mask (1) where the footprints and polygons overlap
var footprint_mask = ee.Image.constant(0).paint(v_foot_prints, 1);
var poly_mask = ee.Image.constant(0).paint(v_srer_polys, 1);
var valid_area_mask = footprint_mask.and(poly_mask);

// Calculate the exact percentage of overlap for every grid cell
var grid_overlap = valid_area_mask.reduceRegions({
  collection: focused_grid,
  reducer: ee.Reducer.mean(),
  scale: 1, 
  crs: projSent2.crs(),
  tileScale: 4
});

// Filter for strict 100% overlap
var final_grid = grid_overlap.filter(ee.Filter.gte('mean', 0.99));

// Transfer the polygon attributes to your perfectly overlapping grid cells
var v_spatial_filter = ee.Filter.intersects({leftField:'.geo', rightField:'.geo', maxError:1});
var v_saveAllJoin = ee.Join.saveAll({matchesKey:'polys'});

var v_sent2_joined_grids = v_saveAllJoin.apply(final_grid, v_srer_polys, v_spatial_filter)
  .map(function (ft){
    var ft1 = ee.Feature(ee.List(ft.get('polys')).get(0));
    return ft.set({
      'Plant_Comm': ft1.get('Plant_Comm'),
      'Pasture': ft1.get('Pasture'),
      'Transect': ft1.get('Transect'),
      'Utility': ft1.get('Utility'),
      'S_Desc': ft1.get('S_Desc'),
      'Exclosure': ft1.get('Exclosure'),
      'area_ha': ft1.get('area_ha'),
      'polys': null 
    });
  });

// =========================================================================
// PART 2: EXTRACT SENTINEL-2 BANDS & INDICES 
// =========================================================================

function extractS2Data(startDate, endDate, monthLabel) {
  var sent2_ic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(v_extent)                     
    .filterDate(startDate, endDate);                

  // Triple Defense Masking + Topographic Cosine Correction
  var sent2_im = sent2_ic
    .map(function(img) {
      // --- Defense 1, 2, 3 ---
      var probMask = img.select('MSK_CLDPRB').lt(20);
      var scl = img.select('SCL');
      var sclMask = scl.neq(8).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11)).and(scl.neq(3));
      var blueMask = img.select('B2').lt(2500);
      var masterMask = probMask.and(sclMask).and(blueMask);
      var maskedImg = img.updateMask(masterMask);

      // --- Topographic Illumination Correction (Cosine) ---
      // Convert sun angles into constant Images to prevent Number vs Image math crashes
      var sz_num = ee.Number(img.get('MEAN_SOLAR_ZENITH_ANGLE')).multiply(Math.PI / 180);
      var sa_num = ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')).multiply(Math.PI / 180);
      
      var cosZ = ee.Image.constant(sz_num.cos());
      var sinZ = ee.Image.constant(sz_num.sin());
      var sa_img = ee.Image.constant(sa_num);
      
      var cosS = terrain_slope.cos();
      var sinS = terrain_slope.sin();
      var cosAzAsp = sa_img.subtract(terrain_aspect).cos();
      
      // Calculate Illumination Condition (IL) securely with Image math
      var illumination = cosZ.multiply(cosS).add(sinZ.multiply(sinS).multiply(cosAzAsp));
      
      // Clamp to avoid extreme overcorrection in deep shadows (prevent dividing by ~0)
      var illu_clamped = illumination.max(0.1); 
      
      var correctionFactor = cosZ.divide(illu_clamped);
      
      // Apply correction only to optical reflectance bands
      var bandsToCorrect = ['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12'];
      
      // THE FIX: Explicitly cast to float so median() can seamlessly stack them!
      var correctedBands = maskedImg.select(bandsToCorrect).multiply(correctionFactor).toFloat();
      
      // Overwrite raw bands with corrected bands
      return maskedImg.addBands(correctedBands, null, true);
    })
    .median()
    .clip(v_extent)
    .select(['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12'])
    .multiply(0.0001) 
    .setDefaultProjection({crs: projSent2.crs(), scale: projSent2.nominalScale()});

  var ndvi = sent2_im.normalizedDifference(['B8', 'B4']).rename('NDVI');

  var mcari = sent2_im.expression(
      '((B5 - B4) - 0.2 * (B5 - B3)) * (B5 / B4)', {
        'B3': sent2_im.select('B3'), 
        'B4': sent2_im.select('B4'), 
        'B5': sent2_im.select('B5')  
  }).rename('MCARI');

  // Calculate Bare Soil Index (BSI)
  var bsi = sent2_im.expression(
      '((B11 + B4) - (B8 + B2)) / ((B11 + B4) + (B8 + B2))', {
        'B2':  sent2_im.select('B2'),
        'B4':  sent2_im.select('B4'),
        'B8':  sent2_im.select('B8'),
        'B11': sent2_im.select('B11')
  }).rename('BSI');

  // Calculate Normalized Burn Ratio 2 (NBR2)
  var nbr2 = sent2_im.normalizedDifference(['B11', 'B12']).rename('NBR2');

  // Add all indices back to the image (Slope removed)
  sent2_im = sent2_im.addBands([ndvi, mcari, bsi, nbr2]);
  
  // Define exactly what to extract so it gets attached to the grid features
  var bandsToExtract = sent2_im.select([
    'B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12', 
    'NDVI', 'MCARI', 'BSI', 'NBR2'
  ]);

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

  // Filter out any cells that fall outside the satellite coverage
  return extracted.filter(ee.Filter.notNull([
    'B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12', 
    'NDVI', 'MCARI', 'BSI', 'NBR2'
  ]));
}

// Pass the variables from the dynamic query into the extraction function
var may_s2_data = extractS2Data(may_start, may_end, 'May');
var sep_s2_data = extractS2Data(sep_start, sep_end, 'Sept');

// Force the grid geometries into WGS84 UTM Zone 12N (Meters)
var may_s2_utm = may_s2_data.map(function(f) { return f.transform('EPSG:32612', 0.05); });
var sep_s2_utm = sep_s2_data.map(function(f) { return f.transform('EPSG:32612', 0.05); });


// =========================================================================
// PART 3: DRONE METRICS (BGR, LPI, MFT)
// =========================================================================

function processMonthMetrics(grid_subset, classified_img) {
  var native_proj = classified_img.projection();
  
  // Define BOTH masks
  var binary = classified_img.eq(3).selfMask();       // Bare Ground Mask
  var obstacles = classified_img.neq(3).selfMask();   // Vegetation/Obstacle Mask

  return grid_subset.map(function(ft){
    
    // --- BGR Calculation ---
    var v_area_image = binary.multiply(ee.Image.pixelArea());
    
    var v_area_ft = ft.area(0.05); 
    
    var v_area = v_area_image.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: ft.geometry(),
      scale: 0.05,
      maxPixels: 1e13
     }).get('classification');
    var v_pct_area = ee.Number(v_area).divide(v_area_ft).multiply(100);

    // --- LPI Calculation ---
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

    // --- MFT Calculation ---
    function Get_Mean_Fetch(obstacle_mask, bare_mask, v_points) {
      var v_distance = obstacle_mask.fastDistanceTransform().sqrt().multiply(ee.Image.pixelArea().sqrt()).rename("distance");
      var fetch_on_bare = v_distance.updateMask(bare_mask);

      v_points = fetch_on_bare.reduceRegions({
        collection: v_points,
        reducer: ee.Reducer.first().setOutputs(["distance"]),
        scale: 0.05
      });
      
      var valid_points = v_points.filter(ee.Filter.notNull(['distance']));
      var raw_mean = valid_points.reduceColumns(ee.Reducer.mean(),['distance']).get('mean');
      var final_mean = ee.Algorithms.If(ee.Algorithms.IsEqual(raw_mean, null), 0, raw_mean);

      return ee.Number(final_mean);
    }

    var N_PTS = 1000;
    var v_rnd = ee.FeatureCollection.randomPoints(ft.geometry(), N_PTS, 1234, 0.05);
    
    var v_nearestMeanValues = Get_Mean_Fetch(obstacles, binary, v_rnd);

    return ft.set('LPI', max_area, 'BGR', v_pct_area, 'MFT', v_nearestMeanValues);
  });
}

// Pass the UTM-transformed data into the metric processor
var final_may = processMonthMetrics(may_s2_utm, v_classified_may);
var final_sep = processMonthMetrics(sep_s2_utm, v_classified_sep);

// Now we merge them for the final export
var final_combined = final_may.merge(final_sep);

// =========================================================================
// EXPORTS
// =========================================================================

Export.table.toDrive({
  collection: final_combined,
  description: 'CSV_table',
  folder: 'GEE_Downloads',
  fileFormat: 'CSV'
});

Export.table.toAsset({
  collection: final_combined,
  assetId: 'projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm',
  description: 'FC_asset'
});


