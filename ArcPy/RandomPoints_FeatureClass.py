import arcpy
import os
import numpy as np
from scipy.ndimage import label, distance_transform_edt
import concurrent.futures
import multiprocessing

# ====================================================================
# USER INPUTS
# ====================================================================
out_gdb = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\MyProject1.gdb'
output_nri_plots = os.path.join(out_gdb, 'SRER_NRI_Plots_110m')
output_grid_10m = os.path.join(out_gdb, 'SRER_Grid_10m')
output_grid_30m = os.path.join(out_gdb, 'SRER_Grid_30m')

# TIFFs
may_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_May_2019_UTM12N_Mosaic.tif'
sep_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_Sep_2019_UTM12N_Mosaic.tif'

# ====================================================================
# WORKER FUNCTION (Runs on multiple cores)
# ====================================================================
def process_grid_cell(data_packet):
  """Runs isolated on a separate CPU core."""
  oid, season, xmin, ymin, xmax, ymax, tiff_path, c_size, is_circle = data_packet
  target_value = 3 
  
  try:
    lower_left = arcpy.Point(xmin, ymin)
    ncols = int(round((xmax - xmin) / c_size))
    nrows = int(round((ymax - ymin) / c_size))
    
    main_array = arcpy.RasterToNumPyArray(tiff_path, lower_left, ncols, nrows)
    
    # --- CIRCULAR MASKING LOGIC ---
    if is_circle:
      Y, X = np.ogrid[:nrows, :ncols]
      center_r, center_c = nrows / 2.0, ncols / 2.0
      dist_from_center = np.sqrt((X - center_c)**2 + (Y - center_r)**2) * c_size
      
      mask = dist_from_center <= 55.0
      is_bare = (main_array == target_value) & mask
      total_pixels = np.sum(mask) 
      main_array[~mask] = 0 
    else:
      is_bare = (main_array == target_value)
      total_pixels = main_array.size
        
    if total_pixels == 0:
      return oid, season, None

    # --- 1. BGR ---
    bgr_percent = (np.sum(is_bare) / total_pixels) * 100

    # --- 2. LPI ---
    structure = np.ones((3, 3), dtype=int)
    labeled_array, num_features = label(is_bare, structure=structure)
    if num_features > 0:
      counts = np.bincount(labeled_array.ravel())
      max_patch_pixels = np.max(counts[1:]) if len(counts) > 1 else 0
      lpi_percent = (max_patch_pixels / total_pixels) * 100
    else:
      lpi_percent = 0.0

    # --- 3. Mean Fetch ---
    dist_array = distance_transform_edt(is_bare) * c_size
    if is_circle:
      mean_fetch_exact = np.mean(dist_array[mask])
    else:
      mean_fetch_exact = np.mean(dist_array)

    # --- 4. CANOPY GAP FRACTIONS ---
    horizontal_transects = [main_array[i, :] for i in range(nrows)]
    vertical_transects = [main_array[:, j] for j in range(ncols)]
    all_transects = horizontal_transects + vertical_transects
    
    total_transect_length_m = len(all_transects) * (ncols * c_size)
    if is_circle:
      total_transect_length_m = 2 * (total_pixels * c_size)

    all_gap_lengths = []
    for t in all_transects:
      is_gap = (t == target_value)
      padded = np.concatenate(([False], is_gap, [False]))
      diffs = np.diff(padded.astype(int))
      starts = np.where(diffs == 1)[0]
      ends = np.where(diffs == -1)[0]
      all_gap_lengths.extend((ends - starts) * c_size)

    all_gap_lengths = np.array(all_gap_lengths)

    if total_transect_length_m > 0:
      f_0_24   = (np.sum(all_gap_lengths[(all_gap_lengths < 0.25)]) / total_transect_length_m) * 100
      f_25_50  = (np.sum(all_gap_lengths[(all_gap_lengths >= 0.25) & (all_gap_lengths <= 0.50)]) / total_transect_length_m) * 100
      f_51_100 = (np.sum(all_gap_lengths[(all_gap_lengths >= 0.51) & (all_gap_lengths <= 1.00)]) / total_transect_length_m) * 100
      f_101_200= (np.sum(all_gap_lengths[(all_gap_lengths >= 1.01) & (all_gap_lengths <= 2.00)]) / total_transect_length_m) * 100
      f_gt_200 = (np.sum(all_gap_lengths[(all_gap_lengths > 2.00)]) / total_transect_length_m) * 100
    else:
      f_0_24 = f_25_50 = f_51_100 = f_101_200 = f_gt_200 = 0.0

    metrics = [bgr_percent, lpi_percent, mean_fetch_exact, f_0_24, f_25_50, f_51_100, f_101_200, f_gt_200]
    return oid, season, metrics
      
  except Exception as e:
    return oid, season, None

# ====================================================================
# EXECUTION LOGIC
# ====================================================================
def calculate_metrics_for_fc(target_fc, cell_size, is_circle=False):
  fc_name = os.path.basename(target_fc)
  print(f'\n--- Processing 2X Seasonal Metrics for {fc_name} ---')

  # 1. Read all geometries into memory
  print('  -> Reading boundaries and generating May/Sep tasks...')
  shapes_dict = {}
  tasks = []
  
  with arcpy.da.SearchCursor(target_fc, ['OID@', 'SHAPE@']) as cursor:
    for oid, geom in cursor:
      shapes_dict[oid] = geom
      ext = geom.extent
      
      # Snap to 10m grid to prevent GEE drift
      min_x = round(ext.XMin / 10.0) * 10.0
      min_y = round(ext.YMin / 10.0) * 10.0
      max_x = round(ext.XMax / 10.0) * 10.0
      max_y = round(ext.YMax / 10.0) * 10.0

      # Create TWO tasks for every single feature
      tasks.append((oid, 'May', min_x, min_y, max_x, max_y, may_tiff, cell_size, is_circle))
      tasks.append((oid, 'Sep', min_x, min_y, max_x, max_y, sep_tiff, cell_size, is_circle))

  total_tasks = len(tasks)
  print(f'  -> Built {total_tasks} total extraction tasks ({len(shapes_dict)} features x 2 seasons).')
  
  # 2. Multiprocessing
  print(f'  -> Spinning up {multiprocessing.cpu_count()} CPU cores...')
  results_list = []
  with concurrent.futures.ProcessPoolExecutor() as executor:
    for count, result in enumerate(executor.map(process_grid_cell, tasks), 1):
      oid, season, metrics = result
      if metrics:
        results_list.append((oid, season, metrics))
      if count % 1000 == 0:
        print(f'  -> Processed {count} / {total_tasks} tasks...')

  # 3. Create a clean replacement feature class to hold the 2X rows
  print('  -> Constructing final 2X feature class...')
  sr = arcpy.Describe(target_fc).spatialReference
  temp_out = r'memory\temp_metrics_fc'
  
  # Clean up memory if it exists from a previous loop
  if arcpy.Exists(temp_out):
    arcpy.management.Delete(temp_out)
      
  arcpy.management.CreateFeatureclass('memory', 'temp_metrics_fc', 'POLYGON', spatial_reference=sr)
  
  # Add fields
  arcpy.management.AddField(temp_out, 'Season', 'TEXT')
  new_fields = ['BGR_Exact', 'LPI_Exact', 'Fetch_Exact', 'Gap_0_24_Exact', 'Gap_25_50_Exact', 'Gap_51_100_Exact', 'Gap_101_200_Exact', 'Gap_gt_200_Exact']
  for field in new_fields:
    arcpy.management.AddField(temp_out, field, 'DOUBLE')

  # 4. Insert the duplicated rows with their respective metrics
  insert_fields = ['SHAPE@', 'Season'] + new_fields
  with arcpy.da.InsertCursor(temp_out, insert_fields) as icursor:
    for oid, season, metrics in results_list:
      geom = shapes_dict[oid]
      row_data = [geom, season] + metrics
      icursor.insertRow(row_data)

  # 5. Overwrite the original feature class with our new 2X version
  arcpy.management.CopyFeatures(temp_out, target_fc)
  arcpy.management.Delete(temp_out)
  
  print(f'  -> Successfully updated {fc_name} with May and September rows!')

if __name__ == '__main__':
  arcpy.ResetEnvironments() 
  arcpy.env.overwriteOutput = True

  # Retrieve raster cell size from the May TIFF
  print('Reading base raster properties...')
  pixel_size = float(arcpy.management.GetRasterProperties(may_tiff, 'CELLSIZEX').getOutput(0))
  print(f'Detected pixel size: {pixel_size}m')

  # Execute for all three datasets
  calculate_metrics_for_fc(output_grid_10m, pixel_size, is_circle=False)
  calculate_metrics_for_fc(output_grid_30m, pixel_size, is_circle=False)
  calculate_metrics_for_fc(output_nri_plots, pixel_size, is_circle=True)

  print('\nAll Extracted Metrics Processing Complete! Your features now contain 2X rows for seasonal comparison.')
