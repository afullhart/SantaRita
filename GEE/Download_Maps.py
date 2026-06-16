import ee
import math
import time

# =========================================================================
# INITIALIZE EARTH ENGINE
# =========================================================================
# Initialize the Earth Engine module. (If it prompts you, run ee.Authenticate() first)
try:
    ee.Initialize()
    print("Earth Engine Initialized Successfully.")
except Exception as e:
    print("Authentication required. Please run `earthengine authenticate` in your terminal.")
    ee.Authenticate()
    ee.Initialize()

# =========================================================================
# SETUP & ASSETS
# =========================================================================
fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_landsat_model_grid_utm')
bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds')
cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass_Landsat')
bounds_geom = bounds_fc.first().geometry().bounds()

dem = ee.Image('projects/ee-andrewfullhart/assets/SR_10m_DEM_Resampled')
terrain_slope = ee.Terrain.slope(dem).multiply(math.pi / 180)
terrain_aspect = ee.Terrain.aspect(dem).multiply(math.pi / 180)

# =========================================================================
# MODEL TRAINING 
# =========================================================================
print("Setting up Gradient Tree Boost models...")

inputProps = ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'BSI', 'NBR2', 'slope', 'illumination', 'aspect']

# Core Dataset
core_training_fc = fc.filter(ee.Filter.notNull(inputProps + ['Herb_pct', 'Woody_pct']))

# HWR Dataset (Log Transformed)
def add_log_hwr(ft):
    raw_hwr = ee.Number(ft.get('Herb_Woody_Ratio'))
    log_hwr = raw_hwr.add(1).log()
    return ft.set('Log_HWR', log_hwr)

hwr_training_fc = core_training_fc.filter(ee.Filter.notNull(['Herb_Woody_Ratio'])).map(add_log_hwr)

hyperpars = {
    'numberOfTrees': 400, 
    'shrinkage': 0.05, 
    'samplingRate': 0.7, 
    'maxNodes': 32, 
    'loss': 'Huber', 
    'seed': 123
}

model_bgr = ee.Classifier.smileGradientTreeBoost(**hyperpars).setOutputMode('REGRESSION').train(core_training_fc, 'BGR', inputProps)
model_lpi = ee.Classifier.smileGradientTreeBoost(**hyperpars).setOutputMode('REGRESSION').train(core_training_fc, 'LPI', inputProps)
model_mft = ee.Classifier.smileGradientTreeBoost(**hyperpars).setOutputMode('REGRESSION').train(core_training_fc, 'MFT', inputProps)
model_herb = ee.Classifier.smileGradientTreeBoost(**hyperpars).setOutputMode('REGRESSION').train(core_training_fc, 'Herb_pct', inputProps)
model_woody = ee.Classifier.smileGradientTreeBoost(**hyperpars).setOutputMode('REGRESSION').train(core_training_fc, 'Woody_pct', inputProps)
model_hwr = ee.Classifier.smileGradientTreeBoost(**hyperpars).setOutputMode('REGRESSION').train(hwr_training_fc, 'Log_HWR', inputProps)

# =========================================================================
# LANDSAT EXTRACTION PIPELINE
# =========================================================================
projLandsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(bounds_geom).first().select('SR_B2').projection()

def maskClouds(qa):
    # In Python, we use the standard integer bitshift operators directly
    dilatedCloud = qa.bitwiseAnd(1 << 1).eq(0)
    cirrus = qa.bitwiseAnd(1 << 2).eq(0)
    cloud = qa.bitwiseAnd(1 << 3).eq(0)
    shadow = qa.bitwiseAnd(1 << 4).eq(0)
    # Use .And() for Earth Engine logical operations in Python
    return dilatedCloud.And(cirrus).And(cloud).And(shadow)

def prepOLI(img):
    scaled = img.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']).multiply(0.0000275).add(-0.2)
    qaMask = maskClouds(img.select('QA_PIXEL'))
    finalMask = qaMask.And(scaled.select('SR_B5').gt(0.15)).focal_min(radius=30, kernelType='square', units='meters')

    sz_num = ee.Number(img.get('SUN_ELEVATION')).subtract(90).multiply(-1).multiply(math.pi / 180)
    sa_num = ee.Number(img.get('SUN_AZIMUTH')).multiply(math.pi / 180)
    
    illumination = ee.Image.constant(sz_num.cos()).multiply(terrain_slope.cos()).add(ee.Image.constant(sz_num.sin()).multiply(terrain_slope.sin()).multiply(ee.Image.constant(sa_num).subtract(terrain_aspect).cos())).rename('illumination')

    optical = scaled.rename(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']).addBands(illumination)
    return optical.updateMask(finalMask)

def prepTM(img):
    scaled = img.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']).multiply(0.0000275).add(-0.2)
    qaMask = maskClouds(img.select('QA_PIXEL'))
    finalMask = qaMask.And(scaled.select('SR_B4').gt(0.15)).focal_min(radius=30, kernelType='square', units='meters')

    sz_num = ee.Number(img.get('SUN_ELEVATION')).subtract(90).multiply(-1).multiply(math.pi / 180)
    sa_num = ee.Number(img.get('SUN_AZIMUTH')).multiply(math.pi / 180)
    
    illumination = ee.Image.constant(sz_num.cos()).multiply(terrain_slope.cos()).add(ee.Image.constant(sz_num.sin()).multiply(terrain_slope.sin()).multiply(ee.Image.constant(sa_num).subtract(terrain_aspect).cos())).rename('illumination')

    optical = scaled.rename(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']).addBands(illumination)
    return optical.updateMask(finalMask)

def buildLandsatComposite(valid_ids_str, year, month):
    fullIdList = ee.String(valid_ids_str).split(',')
    
    def extract_index(id_str):
        parts = ee.String(id_str).split('/')
        return parts.get(parts.length().subtract(1))
        
    indexList = fullIdList.map(extract_index)

    l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filter(ee.Filter.inList('system:index', indexList)).map(prepOLI)
    l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filter(ee.Filter.inList('system:index', indexList)).map(prepOLI)
    l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filter(ee.Filter.inList('system:index', indexList)).map(prepTM)
    l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filter(ee.Filter.inList('system:index', indexList)).map(prepTM)

    combined = l9.merge(l8).merge(l7).merge(l5)
    ls_median = combined.median().clip(bounds_geom)

    # Synthetic Median Replacement for Sept 2019
    isTargetCond = ee.Number(year).eq(2019).And(ee.Number(month).eq(9))
    
    aug31_col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(bounds_geom).filterDate('2019-08-31', '2019-09-01').map(prepOLI)
    oct02_col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(bounds_geom).filterDate('2019-10-02', '2019-10-03').map(prepOLI)

    synthetic_sept = aug31_col.merge(oct02_col).median().clip(bounds_geom)

    ls_median = ee.Image(ee.Algorithms.If(isTargetCond, synthetic_sept, ls_median))

    optical_bands = ls_median.setDefaultProjection(crs=projLandsat.crs(), scale=projLandsat.nominalScale())

    ndvi = optical_bands.normalizedDifference(['NIR', 'Red']).rename('NDVI')
    bsi = optical_bands.expression(
        '((SWIR1 + Red) - (NIR + Blue)) / ((SWIR1 + Red) + (NIR + Blue))', {
            'Blue':  optical_bands.select('Blue'),
            'Red':   optical_bands.select('Red'),
            'NIR':   optical_bands.select('NIR'),
            'SWIR1': optical_bands.select('SWIR1')
        }).rename('BSI')
    nbr2 = optical_bands.normalizedDifference(['SWIR1', 'SWIR2']).rename('NBR2')

    return optical_bands.addBands([ndvi, bsi, nbr2, terrain_slope.rename('slope'), terrain_aspect.rename('aspect')])

# =========================================================================
# AUTOMATED TASK SUBMISSION
# =========================================================================
print("Fetching valid cloud windows from asset (this takes a moment)...")
# Pull the actual data dictionary of the FeatureCollection from the server to Python
windows_list = cloud_windows.getInfo().get('features', [])

print(f"Found {len(windows_list)} valid monthly windows. Submitting export tasks...")

for window in windows_list:
    props = window['properties']
    year = int(props['Year'])
    month = int(props['Month'])
    valid_ids_str = props['Valid_IDs']
    
    print(f"Queueing task for {year}-{month:02d}...")
    
    ls_img = buildLandsatComposite(valid_ids_str, year, month)

    p_bgr = ls_img.classify(model_bgr).rename('Pred_BGR')
    p_lpi = ls_img.classify(model_lpi).rename('Pred_LPI')
    p_mft = ls_img.classify(model_mft).rename('Pred_MFT')
    p_herb = ls_img.classify(model_herb).rename('Pred_Herb_pct')
    p_woody = ls_img.classify(model_woody).rename('Pred_Woody_pct')
    p_hwr = ls_img.classify(model_hwr).rename('Pred_Log_HWR')

    # Stack all 6 prediction bands into a single image to save drive space and task overhead
    final_pred_img = ee.Image.cat([p_bgr, p_lpi, p_mft, p_herb, p_woody, p_hwr]).toFloat()

    task_name = f'SRER_Predictive_Maps_{year}_{month:02d}'
    
    # Configure the Google Drive export task
    task = ee.batch.Export.image.toDrive(
        image=final_pred_img,
        description=task_name,
        folder='SRER_Predictive_Maps', # Will create a folder in your Drive
        fileNamePrefix=task_name,
        region=bounds_geom,
        scale=30,
        crs='EPSG:32612', # UTM Zone 12N
        maxPixels=1e13
    )
    
    # Start the task on Google's servers
    task.start()
    time.sleep(0.2) # Small delay to prevent hitting the API request rate limit

print("All tasks have been successfully submitted!")
print("You can monitor their progress at: https://code.earthengine.google.com/tasks")
