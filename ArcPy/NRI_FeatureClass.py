import arcpy
import numpy as np
from arcpy.sa import *

# ====================================================================
# CONFIGURATION
# ====================================================================
may_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_May_2019_UTM12N_Mosaic.tif'
sep_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_Sep_2019_UTM12N_Mosaic.tif'
input_plots = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\MyProject1.gdb\SRER_NRI_Plots'
csv_folder = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1'

target_value = 3  # Class value for Bare Ground / Canopy Gap
herb_value = 1    # Class value for Herb Cover
woody_value = 2   # Class value for Woody Cover

arcpy.ResetEnvironments()
arcpy.env.workspace = 'memory'
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('Spatial')

cell_size = float(arcpy.management.GetRasterProperties(may_tiff, 'CELLSIZEX').getOutput(0))

# Helper function to extract gap lengths
def get_gap_lengths(transect_array, gap_val, p_size):
  padded = np.concatenate(([False], (transect_array == gap_val), [False]))
  diffs = np.diff(padded.astype(int))
  starts, ends = np.where(diffs == 1)[0], np.where(diffs == -1)[0]
  return (ends - starts) * p_size

# ====================================================================
# DATABASE PREPARATION & TEMPORAL DUPLICATION
# ====================================================================
print('Preparing feature class schema and duplicating geometries for temporal analysis...')

# 1. Add Month Field
existing_fields = [f.name for f in arcpy.ListFields(input_plots)]
if 'Month' not in existing_fields:
  arcpy.management.AddField(input_plots, 'Month', 'TEXT', field_length=15)

# 2. Add Metric Fields (Added NRI Herb-to-Woody Ratios)
new_fields = [
  ('Exact_BGR_Pct', 'DOUBLE'), ('Exact_LPI_Pct', 'DOUBLE'), 
  ('Exact_Fetch_m', 'DOUBLE'), ('Exact_Herb_Woody_Ratio', 'DOUBLE'), 
  ('Exact_Gap_0_24', 'DOUBLE'), ('Exact_Gap_25_50', 'DOUBLE'), ('Exact_Gap_51_100', 'DOUBLE'),
  ('Exact_Gap_101_200', 'DOUBLE'), ('Exact_Gap_gt_200', 'DOUBLE'), ('ExactTotal_Gap', 'DOUBLE'),
  ('NRI_BGR_0cm_Pct', 'DOUBLE'), ('NRI_BGR_50cm_Pct', 'DOUBLE'), ('NRI_BGR_100cm_Pct', 'DOUBLE'),
  ('NRI_HW_Ratio_0cm', 'DOUBLE'), ('NRI_HW_Ratio_50cm', 'DOUBLE'), ('NRI_HW_Ratio_100cm', 'DOUBLE'),
  ('NRI_Fetch_50cm', 'DOUBLE'), ('NRI_Fetch_100cm', 'DOUBLE'),
  ('NRI_Gap_0_24', 'DOUBLE'), ('NRI_Gap_25_50', 'DOUBLE'),
  ('NRI_Gap_51_100', 'DOUBLE'), ('NRI_Gap_101_200', 'DOUBLE'), ('NRI_Gap_gt_200', 'DOUBLE'),
  ('NRI_Total_Gap', 'DOUBLE')
]

for field_name, field_type in new_fields:
  if field_name not in existing_fields:
    arcpy.management.AddField(input_plots, field_name, field_type)

# 3. Duplicate Rows for May & September
may_count, sep_count = 0, 0
with arcpy.da.SearchCursor(input_plots, ['Month']) as cursor:
  for row in cursor:
    if row[0] == 'May_2019': may_count += 1
    if row[0] == 'Sep_2019': sep_count += 1

if may_count == 0 and sep_count == 0:
  print('  -> First run detected. Generating temporal pairs...')
  sept_shapes = []
  
  with arcpy.da.UpdateCursor(input_plots, ['SHAPE@', 'Month']) as cursor:
    for row in cursor:
      row[1] = 'May_2019'
      sept_shapes.append(row[0])
      cursor.updateRow(row)
      
  with arcpy.da.InsertCursor(input_plots, ['SHAPE@', 'Month']) as cursor:
    for shape in sept_shapes:
      cursor.insertRow([shape, 'Sep_2019'])
else:
  print('  -> Temporal pairs already exist. Proceeding to extraction...')

cursor_fields = ['SHAPE@', 'Month'] + [f[0] for f in new_fields]
total_plots = int(arcpy.management.GetCount(input_plots).getOutput(0))

# ====================================================================
# MAIN PROCESSING LOOP
# ====================================================================
print(f'\nStarting metric extraction for all {total_plots} temporal plots...')

with arcpy.da.UpdateCursor(input_plots, cursor_fields) as cursor:
  for i, row in enumerate(cursor):
    plot_geom = row[0]
    month_val = row[1]
    
    active_tiff = may_tiff if month_val == 'May_2019' else sep_tiff
    print(f'  -> Processing Plot {i + 1} of {total_plots} [{month_val}]...')
    
    arcpy.env.snapRaster = active_tiff
    arcpy.env.cellSize = active_tiff

    try:
      # EXTRACT EXACT PIXELS 
      masked_raster = ExtractByMask(active_tiff, plot_geom)
      main_array = arcpy.RasterToNumPyArray(masked_raster, nodata_to_value=-9999)
      
      valid_mask = main_array != -9999
      total_valid_pixels = np.sum(valid_mask)
      
      if total_valid_pixels == 0:
        print(f'     [Warning] Plot {i + 1} contains no valid raster data. Skipping.')
        continue

      # ==============================================================
      # SUITE 1: EXACT 2D METRICS
      # ==============================================================
      bare_ground_pixels = np.sum(main_array[valid_mask] == target_value)
      exact_bgr = (float(bare_ground_pixels) / float(total_valid_pixels)) * 100

      isolated_class = Con((masked_raster == target_value), 1)
      grouped_patches = RegionGroup(isolated_class, 'EIGHT', 'WITHIN', 'NO_LINK')
      patch_array = arcpy.RasterToNumPyArray(grouped_patches, nodata_to_value=0)
      
      patch_ids, counts = np.unique(patch_array, return_counts=True)
      valid_patch_mask = patch_ids > 0 
      max_patch_pixels = np.max(counts[valid_patch_mask]) if np.any(valid_patch_mask) else 0
      exact_lpi = (float(max_patch_pixels) / float(total_valid_pixels)) * 100

      # HERB-TO-WOODY RATIO
      herb_pixels = np.sum(main_array[valid_mask] == herb_value)
      woody_pixels = np.sum(main_array[valid_mask] == woody_value)
      exact_herb_woody = (float(herb_pixels) / float(woody_pixels)) if woody_pixels > 0 else None

      # Calculate Distance from Vegetation (Obstacles)
      obstacle_class = Con((masked_raster != target_value) & (~IsNull(masked_raster)), 1)
      fetch_dist_raster = EucDistance(obstacle_class)
      
      # EXACT FETCH
      full_fetch_array = arcpy.RasterToNumPyArray(ExtractByMask(fetch_dist_raster, plot_geom), nodata_to_value=-9999)
      valid_full_fetch = full_fetch_array[full_fetch_array >= 0]
      exact_fetch = np.mean(valid_full_fetch) if valid_full_fetch.size > 0 else 0.0

      # ==============================================================
      # SUITE 2: EXHAUSTIVE VIRTUAL TRANSECTS (ALL ROWS/COLS)
      # ==============================================================
      exhaust_gap_arrays = []
      for r in range(main_array.shape[0]):
        exhaust_gap_arrays.append(get_gap_lengths(main_array[r, :], target_value, cell_size))
      for c in range(main_array.shape[1]):
        exhaust_gap_arrays.append(get_gap_lengths(main_array[:, c], target_value, cell_size))
          
      all_exhaust_gaps = np.concatenate(exhaust_gap_arrays)
      total_exhaust_length_m = total_valid_pixels * cell_size * 2
      
      ex_0_24 = (np.sum(all_exhaust_gaps[all_exhaust_gaps < 0.25]) / total_exhaust_length_m) * 100
      ex_25_50 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 0.25) & (all_exhaust_gaps <= 0.50)]) / total_exhaust_length_m) * 100
      ex_51_100 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 0.51) & (all_exhaust_gaps <= 1.00)]) / total_exhaust_length_m) * 100
      ex_101_200 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 1.01) & (all_exhaust_gaps <= 2.00)]) / total_exhaust_length_m) * 100
      ex_gt_200 = (np.sum(all_exhaust_gaps[all_exhaust_gaps > 2.00]) / total_exhaust_length_m) * 100
      ex_total = ex_0_24 + ex_25_50 + ex_51_100 + ex_101_200 + ex_gt_200

      # ==============================================================
      # SUITE 3: SYNTHETIC NRI 1D METRICS (3 SPOKES)
      # ==============================================================
      cx, cy = plot_geom.centroid.X, plot_geom.centroid.Y
      xmin, ymax = masked_raster.extent.XMin, masked_raster.extent.YMax
      
      transect_length_m = 50.0
      start_buffer_m = 5.0
      distances = np.linspace(start_buffer_m, start_buffer_m + transect_length_m, int(transect_length_m / cell_size))
      
      all_nri_transects = []
      all_fetch_50cm = []
      all_fetch_100cm = []
      
      spoke_bgr_pixels = 0
      spoke_bgr_pixels_50cm = 0
      spoke_bgr_pixels_100cm = 0

      spoke_herb_pixels = 0
      spoke_woody_pixels = 0
      spoke_herb_pixels_50cm = 0
      spoke_woody_pixels_50cm = 0
      spoke_herb_pixels_100cm = 0
      spoke_woody_pixels_100cm = 0
      
      total_pins_50cm = 0
      total_pins_100cm = 0
      
      for az in [0, 120, 240]:
        az_rad = np.radians(az)
        dx = distances * np.sin(az_rad)
        dy = distances * np.cos(az_rad)
        
        cols = np.round(((cx + dx) - xmin) / cell_size).astype(int)
        rows = np.round((ymax - (cy + dy)) / cell_size).astype(int)
        
        rows = np.clip(rows, 0, main_array.shape[0] - 1)
        cols = np.clip(cols, 0, main_array.shape[1] - 1)
        
        transect_values = main_array[rows, cols]
        all_nri_transects.append(transect_values)
        
        fetch_transect_values = full_fetch_array[rows, cols]
        
        # --- 0cm Continuous Counts ---
        spoke_bgr_pixels += np.sum(transect_values == target_value)
        spoke_herb_pixels += np.sum(transect_values == herb_value)
        spoke_woody_pixels += np.sum(transect_values == woody_value)
        
        # --- 50cm Interval Measurement ---
        vals_50cm = transect_values[::10]
        spoke_bgr_pixels_50cm += np.sum(vals_50cm == target_value)
        spoke_herb_pixels_50cm += np.sum(vals_50cm == herb_value)
        spoke_woody_pixels_50cm += np.sum(vals_50cm == woody_value)
        total_pins_50cm += len(vals_50cm)
        
        # --- 100cm Interval Measurement ---
        vals_100cm = transect_values[::20]
        spoke_bgr_pixels_100cm += np.sum(vals_100cm == target_value)
        spoke_herb_pixels_100cm += np.sum(vals_100cm == herb_value)
        spoke_woody_pixels_100cm += np.sum(vals_100cm == woody_value)
        total_pins_100cm += len(vals_100cm)
        
        # INTERVAL FETCH 
        all_fetch_50cm.extend(fetch_transect_values[::10])
        all_fetch_100cm.extend(fetch_transect_values[::20])
      
      # --- Calculate 1D BGR percentages ---
      nri_bgr = (float(spoke_bgr_pixels) / (len(distances) * 3)) * 100
      nri_bgr_50cm = (float(spoke_bgr_pixels_50cm) / float(total_pins_50cm)) * 100 if total_pins_50cm > 0 else 0.0
      nri_bgr_100cm = (float(spoke_bgr_pixels_100cm) / float(total_pins_100cm)) * 100 if total_pins_100cm > 0 else 0.0

      # --- Calculate 1D Herb-to-Woody Ratios ---
      nri_hw_0cm = (float(spoke_herb_pixels) / float(spoke_woody_pixels)) if spoke_woody_pixels > 0 else None
      nri_hw_50cm = (float(spoke_herb_pixels_50cm) / float(spoke_woody_pixels_50cm)) if spoke_woody_pixels_50cm > 0 else None
      nri_hw_100cm = (float(spoke_herb_pixels_100cm) / float(spoke_woody_pixels_100cm)) if spoke_woody_pixels_100cm > 0 else None
      
      # --- Calculate 1D Interval Fetch ---
      f50_valid = np.array(all_fetch_50cm)
      f50_valid = f50_valid[f50_valid >= 0] 
      nri_fetch_50cm = float(np.mean(f50_valid)) if f50_valid.size > 0 else 0.0
      
      f100_valid = np.array(all_fetch_100cm)
      f100_valid = f100_valid[f100_valid >= 0] 
      nri_fetch_100cm = float(np.mean(f100_valid)) if f100_valid.size > 0 else 0.0
      
      # --- Calculate Gap Metrics ---
      total_nri_length_m = 3 * transect_length_m
      all_nri_gaps = np.concatenate([get_gap_lengths(t, target_value, cell_size) for t in all_nri_transects])
      
      nri_0_24 = (np.sum(all_nri_gaps[all_nri_gaps < 0.25]) / total_nri_length_m) * 100
      nri_25_50 = (np.sum(all_nri_gaps[(all_nri_gaps >= 0.25) & (all_nri_gaps <= 0.50)]) / total_nri_length_m) * 100
      nri_51_100 = (np.sum(all_nri_gaps[(all_nri_gaps >= 0.51) & (all_nri_gaps <= 1.00)]) / total_nri_length_m) * 100
      nri_101_200 = (np.sum(all_nri_gaps[(all_nri_gaps >= 1.01) & (all_nri_gaps <= 2.00)]) / total_nri_length_m) * 100
      nri_gt_200 = (np.sum(all_nri_gaps[all_nri_gaps > 2.00]) / total_nri_length_m) * 100
      nri_total = nri_0_24 + nri_25_50 + nri_51_100 + nri_101_200 + nri_gt_200

      # ==============================================================
      # WRITE DATA BACK TO ROW (Shifted indices for new metrics)
      # ==============================================================
      row[2] = exact_bgr; row[3] = exact_lpi; row[4] = exact_fetch; row[5] = exact_herb_woody
      row[6] = ex_0_24; row[7] = ex_25_50; row[8] = ex_51_100; row[9] = ex_101_200; row[10] = ex_gt_200; row[11] = ex_total
      
      row[12] = nri_bgr
      row[13] = nri_bgr_50cm
      row[14] = nri_bgr_100cm
      
      row[15] = nri_hw_0cm
      row[16] = nri_hw_50cm
      row[17] = nri_hw_100cm
      
      row[18] = nri_fetch_50cm
      row[19] = nri_fetch_100cm
      
      row[20] = nri_0_24; row[21] = nri_25_50; row[22] = nri_51_100; row[23] = nri_101_200; row[24] = nri_gt_200; row[25] = nri_total
      
      cursor.updateRow(row)

    except Exception as e:
      print(f'     [Error] Failed on Plot {i + 1}: {e}')
    
    finally:
      if 'masked_raster' in locals(): arcpy.management.Delete(masked_raster)
      if 'isolated_class' in locals(): arcpy.management.Delete(isolated_class)
      if 'grouped_patches' in locals(): arcpy.management.Delete(grouped_patches)
      if 'obstacle_class' in locals(): arcpy.management.Delete(obstacle_class)
      if 'fetch_dist_raster' in locals(): arcpy.management.Delete(fetch_dist_raster)

print('\nSuccess! All plot attributes have been updated successfully.')

arcpy.conversion.TableToTable(input_plots, csv_folder, 'SRER_NRI_Plots.csv')
print("Export complete!")

arcpy.CheckInExtension('Spatial')
