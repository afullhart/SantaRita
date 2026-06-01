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

may_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_May_2019_UTM12N_Mosaic.tif'
sep_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_Sep_2019_UTM12N_Mosaic.tif'

# --- SAMPLING DESIGN INCREMENTS ---
PT_INCS = [1, 3, 5, 7, 10, 15, 20, 30, 50, 70, 100, 200, 300, 500, 700, 1000, 10000]
LN_INCS = [2, 3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 3000]

# ====================================================================
# WORKER FUNCTION (Runs on multiple cores)
# ====================================================================
def process_grid_cell(data_packet):
  oid, season, xmin, ymin, xmax, ymax, tiff_path, c_size, is_circle = data_packet
  
  target_value = 3  
  herb_value = 1
  woody_value = 2
  
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
      mask = np.ones_like(main_array, dtype=bool)
      is_bare = (main_array == target_value)
      total_pixels = main_array.size
        
    if total_pixels == 0:
      return oid, season, None

    # ==========================================
    # 1. EXACT METRICS CALCULATION
    # ==========================================
    bgr_percent = (np.sum(is_bare) / total_pixels) * 100

    structure = np.ones((3, 3), dtype=int)
    labeled_array, num_features = label(is_bare, structure=structure)
    if num_features > 0:
      counts = np.bincount(labeled_array.ravel())
      max_patch_pixels = np.max(counts[1:]) if len(counts) > 1 else 0
      lpi_percent = (max_patch_pixels / total_pixels) * 100
    else:
      lpi_percent = 0.0

    dist_array = distance_transform_edt(is_bare) * c_size
    mean_fetch_exact = np.mean(dist_array[mask]) if is_circle else np.mean(dist_array)

    # --- HERB-TO-WOODY EXACT RATIO ---
    if is_circle:
      herb_px = np.sum((main_array == herb_value) & mask)
      woody_px = np.sum((main_array == woody_value) & mask)
    else:
      herb_px = np.sum(main_array == herb_value)
      woody_px = np.sum(main_array == woody_value)

    hw_ratio_exact = (float(herb_px) / float(woody_px)) if woody_px > 0 else None

    horizontal_transects = [main_array[i, :] for i in range(nrows)]
    vertical_transects = [main_array[:, j] for j in range(ncols)]
    all_transects = horizontal_transects + vertical_transects
    
    total_transect_length_m = 2 * (total_pixels * c_size) if is_circle else len(all_transects) * (ncols * c_size)
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

    output_metrics = [bgr_percent, lpi_percent, mean_fetch_exact, hw_ratio_exact, f_0_24, f_25_50, f_51_100, f_101_200, f_gt_200]

    # ==========================================
    # 2. RANDOM POINT UNDERSAMPLING 
    # ==========================================
    valid_bare = is_bare[mask]
    valid_fetch = dist_array[mask]
    
    # Establish valid herb/woody matrices for random sampling
    if is_circle:
      valid_herb = (main_array == herb_value)[mask]
      valid_woody = (main_array == woody_value)[mask]
    else:
      valid_herb = (main_array == herb_value).flatten()
      valid_woody = (main_array == woody_value).flatten()

    num_valid_pts = len(valid_bare)

    for n_pts in PT_INCS:
      if num_valid_pts == 0:
        output_metrics.extend([0.0, 0.0, None])
        continue
      
      idx = np.random.choice(num_valid_pts, size=n_pts, replace=True)
      
      # Sample BGR & Fetch
      smpl_bgr = (np.sum(valid_bare[idx]) / n_pts) * 100
      smpl_fetch = np.mean(valid_fetch[idx])
      
      # Sample Herb & Woody to calculate ratio
      smpl_herb_count = np.sum(valid_herb[idx])
      smpl_woody_count = np.sum(valid_woody[idx])
      smpl_hw_ratio = (float(smpl_herb_count) / float(smpl_woody_count)) if smpl_woody_count > 0 else None
      
      output_metrics.extend([smpl_bgr, smpl_fetch, smpl_hw_ratio])

    # ==========================================
    # 3. VIRTUAL TRANSECT UNDERSAMPLING (GAP)
    # ==========================================
    valid_y, valid_x = np.where(mask)
    num_valid_starts = len(valid_y)

    for L_meters in LN_INCS:
      target_px = int(round(L_meters / c_size))
      collected_px = 0
      chunks = []
      
      while collected_px < target_px and num_valid_starts > 0:
        r_idx = np.random.randint(0, num_valid_starts)
        sy, sx = valid_y[r_idx], valid_x[r_idx]
        direction = np.random.randint(0, 4)
        rem_px = target_px - collected_px
        
        if direction == 0: 
          chunk = main_array[sy, sx : min(sx + rem_px, ncols)]
          c_mask = mask[sy, sx : min(sx + rem_px, ncols)]
        elif direction == 1: 
          chunk = main_array[sy : min(sy + rem_px, nrows), sx]
          c_mask = mask[sy : min(sy + rem_px, nrows), sx]
        elif direction == 2: 
          chunk = main_array[sy, max(0, sx - rem_px + 1) : sx + 1][::-1]
          c_mask = mask[sy, max(0, sx - rem_px + 1) : sx + 1][::-1]
        else: 
          chunk = main_array[max(0, sy - rem_px + 1) : sy + 1, sx][::-1]
          c_mask = mask[max(0, sy - rem_px + 1) : sy + 1, sx][::-1]
        
        invalid_idx = np.where(~c_mask)[0]
        if len(invalid_idx) > 0:
          chunk = chunk[:invalid_idx[0]]
          
        if len(chunk) > 0:
          chunks.append(chunk == target_value)
          chunks.append(np.array([False])) 
          collected_px += len(chunk)

      if len(chunks) > 0:
        full_transect = np.concatenate(chunks)
        actual_L = collected_px * c_size
        
        padded_t = np.concatenate(([False], full_transect, [False]))
        diffs_t = np.diff(padded_t.astype(int))
        starts_t = np.where(diffs_t == 1)[0]
        ends_t = np.where(diffs_t == -1)[0]
        gaps_t = (ends_t - starts_t) * c_size
        
        s_0_24   = (np.sum(gaps_t[(gaps_t < 0.25)]) / actual_L) * 100 if actual_L > 0 else 0
        s_25_50  = (np.sum(gaps_t[(gaps_t >= 0.25) & (gaps_t <= 0.50)]) / actual_L) * 100 if actual_L > 0 else 0
        s_51_100 = (np.sum(gaps_t[(gaps_t >= 0.51) & (gaps_t <= 1.00)]) / actual_L) * 100 if actual_L > 0 else 0
        s_101_200= (np.sum(gaps_t[(gaps_t >= 1.01) & (gaps_t <= 2.00)]) / actual_L) * 100 if actual_L > 0 else 0
        s_gt_200 = (np.sum(gaps_t[(gaps_t > 2.00)]) / actual_L) * 100 if actual_L > 0 else 0
        
        output_metrics.extend([s_0_24, s_25_50, s_51_100, s_101_200, s_gt_200])
      else:
        output_metrics.extend([0.0, 0.0, 0.0, 0.0, 0.0])

    return oid, season, output_metrics
      
  except Exception as e:
    return oid, season, None

# ====================================================================
# EXECUTION LOGIC
# ====================================================================
def calculate_metrics_for_fc(target_fc, cell_size, is_circle=False):
  fc_name = os.path.basename(target_fc)
  print(f'\n--- Processing 2X Seasonal Metrics & Undersampling for {fc_name} ---')

  shapes_dict = {}
  tasks = []
  
  with arcpy.da.SearchCursor(target_fc, ['OID@', 'SHAPE@']) as cursor:
    for oid, geom in cursor:
      shapes_dict[oid] = geom
      ext = geom.extent
      
      min_x = round(ext.XMin / 10.0) * 10.0
      min_y = round(ext.YMin / 10.0) * 10.0
      max_x = round(ext.XMax / 10.0) * 10.0
      max_y = round(ext.YMax / 10.0) * 10.0

      tasks.append((oid, 'May', min_x, min_y, max_x, max_y, may_tiff, cell_size, is_circle))
      tasks.append((oid, 'Sep', min_x, min_y, max_x, max_y, sep_tiff, cell_size, is_circle))

  total_tasks = len(tasks)
  print(f'  -> Built {total_tasks} extraction tasks.')
  
  print(f'  -> Spinning up {multiprocessing.cpu_count()} CPU cores...')
  results_list = []
  with concurrent.futures.ProcessPoolExecutor() as executor:
    for count, result in enumerate(executor.map(process_grid_cell, tasks), 1):
      oid, season, metrics = result
      if metrics:
        results_list.append((oid, season, metrics))
      if count % 1000 == 0:
        print(f'  -> Processed {count} / {total_tasks} tasks...')

  print('  -> Constructing final feature class with all simulation fields...')
  sr = arcpy.Describe(target_fc).spatialReference
  temp_out = r'memory\temp_metrics_fc'
  
  if arcpy.Exists(temp_out):
    arcpy.management.Delete(temp_out)
      
  arcpy.management.CreateFeatureclass('memory', 'temp_metrics_fc', 'POLYGON', spatial_reference=sr)
  
  arcpy.management.AddField(temp_out, 'Season', 'TEXT')
  
  # 1. Exact Fields
  new_fields = [
    'BGR_Exact', 'LPI_Exact', 'Fetch_Exact', 'HW_Ratio_Exact',
    'Gap_0_24_Exact', 'Gap_25_50_Exact', 'Gap_51_100_Exact', 'Gap_101_200_Exact', 'Gap_gt_200_Exact'
  ]
  
  # 2. Point Increment Fields (Now includes HW_Ratio_pt)
  for pt in PT_INCS:
    new_fields.extend([f'BGR_pt_{pt}', f'Fetch_pt_{pt}', f'HW_Ratio_pt_{pt}'])
    
  # 3. Line Increment Fields
  for ln in LN_INCS:
    new_fields.extend([
      f'Gap_0_24_L_{ln}', f'Gap_25_50_L_{ln}', f'Gap_51_100_L_{ln}', 
      f'Gap_101_200_L_{ln}', f'Gap_gt_200_L_{ln}'
    ])

  for field in new_fields:
    arcpy.management.AddField(temp_out, field, 'DOUBLE')

  insert_fields = ['SHAPE@', 'Season'] + new_fields
  with arcpy.da.InsertCursor(temp_out, insert_fields) as icursor:
    for oid, season, metrics in results_list:
      geom = shapes_dict[oid]
      row_data = [geom, season] + metrics
      icursor.insertRow(row_data)

  arcpy.management.CopyFeatures(temp_out, target_fc)
  arcpy.management.Delete(temp_out)
  
  print(f'  -> Successfully updated {fc_name}!')

if __name__ == '__main__':
  arcpy.ResetEnvironments() 
  arcpy.env.overwriteOutput = True

  print('Reading base raster properties...')
  pixel_size = float(arcpy.management.GetRasterProperties(may_tiff, 'CELLSIZEX').getOutput(0))
  print(f'Detected pixel size: {pixel_size}m')

  calculate_metrics_for_fc(output_grid_10m, pixel_size, is_circle=False)
  calculate_metrics_for_fc(output_grid_30m, pixel_size, is_circle=False)
  calculate_metrics_for_fc(output_nri_plots, pixel_size, is_circle=True)

  print('\nAll Extracted Metrics & Simulations Processing Complete!')

  csv_out_folder = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data'

  output_grid_10m = os.path.join(out_gdb, 'SRER_Grid_10m')
  output_grid_30m = os.path.join(out_gdb, 'SRER_Grid_30m')
  output_nri_plots = os.path.join(out_gdb, 'SRER_NRI_Plots_110m')

  print('\n--- Exporting Tables to CSV ---')

  export_tasks = {
    output_grid_10m: 'SRER_Grid_10m_Metrics.csv',
    output_grid_30m: 'SRER_Grid_30m_Metrics.csv',
    output_nri_plots: 'SRER_NRI_Plots_110m_Metrics.csv'
  }

  for fc, csv_name in export_tasks.items():
    out_csv_path = os.path.join(csv_out_folder, csv_name)
    
    if arcpy.Exists(out_csv_path):
      arcpy.management.Delete(out_csv_path)
      
    print(f'  -> Exporting {os.path.basename(fc)} to {csv_name}...')
    
    arcpy.conversion.ExportTable(
      in_table=fc,
      out_table=out_csv_path
    )

  print(f'\nSuccess! All metrics exported to: {csv_out_folder}')
  
