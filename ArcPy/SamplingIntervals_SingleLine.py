import os
import arcpy
import numpy as np
from scipy.ndimage import label, distance_transform_edt
import concurrent.futures
import multiprocessing

# ====================================================================
# CONFIGURATION
# ====================================================================
may_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_May_2019_UTM12N_Mosaic.tif'
sep_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_Sep_2019_UTM12N_Mosaic.tif'
csv_folder = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1'
out_gdb = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\MyProject1.gdb'

# Updated suffix to reflect the single line methodology
feature_classes_to_process = {
    os.path.join(out_gdb, 'SRER_NRI_Plots_110m'): 'SRER_NRI_Plots_110m_Single_Line_Convergence.csv',
    os.path.join(out_gdb, 'SRER_Grid_30m'): 'SRER_Grid_30m_Single_Line_Convergence.csv',
    os.path.join(out_gdb, 'SRER_Grid_10m'): 'SRER_Grid_10m_Single_Line_Convergence.csv'
}

target_value = 3  
herb_value = 1    
woody_value = 2   

# The identical physical intervals we are maintaining for the straight line
LINE_INTERVALS_CM = [500, 250, 200, 100, 50, 25, 10, 5]

# ====================================================================
# WORKER FUNCTION
# ====================================================================
def process_polygon(task):
    oid, cx, cy, xmin, ymin, xmax, ymax, active_tiff, cell_size, is_circle = task
    
    try:
        lower_left = arcpy.Point(xmin, ymin)
        ncols = int(round((xmax - xmin) / cell_size))
        nrows = int(round((ymax - ymin) / cell_size))
        
        main_array = arcpy.RasterToNumPyArray(active_tiff, lower_left, ncols, nrows, nodata_to_value=-9999)
        valid_mask = main_array != -9999
        
        if is_circle:
            Y, X = np.ogrid[:nrows, :ncols]
            center_col = (cx - xmin) / cell_size
            center_row = (ymax - cy) / cell_size
            dist_from_center = np.sqrt((X - center_col)**2 + (Y - center_row)**2)
            poly_mask = dist_from_center <= (55.0 / cell_size)
            valid_mask = valid_mask & poly_mask

        total_valid_pixels = np.sum(valid_mask)
        if total_valid_pixels == 0: return oid, None

        # --- EXACT METRICS ---
        bare_ground_pixels = np.sum(main_array[valid_mask] == target_value)
        exact_bgr = (float(bare_ground_pixels) / float(total_valid_pixels)) * 100

        herb_pixels = np.sum(main_array[valid_mask] == herb_value)
        woody_pixels = np.sum(main_array[valid_mask] == woody_value)
        exact_herb = (float(herb_pixels) / float(total_valid_pixels)) * 100
        exact_woody = (float(woody_pixels) / float(total_valid_pixels)) * 100
        exact_hw = (float(herb_pixels) / float(woody_pixels)) if woody_pixels > 0 else None

        is_bare_full = (main_array == target_value)
        dist_array = distance_transform_edt(is_bare_full) * cell_size
        valid_full_fetch = dist_array[valid_mask]
        exact_fetch = np.mean(valid_full_fetch) if valid_full_fetch.size > 0 else 0.0

        results = [exact_bgr, exact_herb, exact_woody, exact_hw, exact_fetch]

        # ==============================================================
        # --- SYNTHETIC SINGLE NORTH-SOUTH TRANSECT ---
        # ==============================================================
        line_length_m = 110.0
        half_length = line_length_m / 2.0
        
        # Generate base continuous pixel distances from -55m to +55m
        distances = np.linspace(-half_length, half_length, int(line_length_m / cell_size))
        
        # Dictionary to store stats for the single line
        line_data = {inv: {'bgr': 0, 'herb': 0, 'woody': 0, 'pins': 0, 'fetch_vals': []} for inv in LINE_INTERVALS_CM}
        
        # North-South orientation: azimuth = 0
        # sin(0) = 0 (dx), cos(0) = 1 (dy)
        dx = distances * 0.0
        dy = distances * 1.0
        
        cols = np.clip(np.round(((cx + dx) - xmin) / cell_size).astype(int), 0, ncols - 1)
        rows = np.clip(np.round((ymax - (cy + dy)) / cell_size).astype(int), 0, nrows - 1)
        
        transect_values = main_array[rows, cols]
        fetch_transect_values = dist_array[rows, cols]
        
        # Step through the continuous transect using fixed field intervals
        for interval in LINE_INTERVALS_CM:
            step = max(1, int((interval / 100.0) / cell_size))
            vals = transect_values[::step]
            
            line_data[interval]['bgr'] += np.sum(vals == target_value)
            line_data[interval]['herb'] += np.sum(vals == herb_value)
            line_data[interval]['woody'] += np.sum(vals == woody_value)
            line_data[interval]['pins'] += len(vals)
            line_data[interval]['fetch_vals'].extend(fetch_transect_values[::step])
        
        # Calculate percentages and append to results
        for interval in LINE_INTERVALS_CM:
            pins = line_data[interval]['pins']
            bgr_pct = (float(line_data[interval]['bgr']) / pins) * 100 if pins > 0 else 0.0
            herb_pct = (float(line_data[interval]['herb']) / pins) * 100 if pins > 0 else 0.0
            woody_pct = (float(line_data[interval]['woody']) / pins) * 100 if pins > 0 else 0.0
            hw_ratio = (float(line_data[interval]['herb']) / float(line_data[interval]['woody'])) if line_data[interval]['woody'] > 0 else None
            
            f_vals = np.array(line_data[interval]['fetch_vals'])
            mean_fetch = float(np.mean(f_vals[f_vals >= 0])) if f_vals.size > 0 else 0.0
            
            results.extend([bgr_pct, herb_pct, woody_pct, hw_ratio, mean_fetch, pins])

        return oid, results
    except Exception as e:
        return oid, None

# ====================================================================
# MAIN EXECUTION BLOCK
# ====================================================================
if __name__ == '__main__':
    arcpy.ResetEnvironments()
    arcpy.CheckOutExtension('Spatial')
    cell_size = float(arcpy.management.GetRasterProperties(may_tiff, 'CELLSIZEX').getOutput(0))
    cores = multiprocessing.cpu_count()

    # Define dynamic fields (Updated prefix to 'Line_')
    new_fields = ['Exact_BGR_Pct', 'Exact_Herb_Pct', 'Exact_Woody_Pct', 'Exact_HW_Ratio', 'Exact_Fetch_m']
    for interval in LINE_INTERVALS_CM:
        new_fields.extend([
            f'Line_BGR_{interval}cm_Pct', f'Line_Herb_{interval}cm_Pct', f'Line_Woody_{interval}cm_Pct',
            f'Line_HW_Ratio_{interval}cm', f'Line_Fetch_{interval}cm', f'Line_TotalPins_{interval}cm'
        ])

    for current_fc, csv_name in feature_classes_to_process.items():
        print(f'\n--- Processing Convergence for: {os.path.basename(current_fc)} ---')
        if not arcpy.Exists(current_fc): continue

        temp_fc = r'memory\temp_convergence'
        if arcpy.Exists(temp_fc): arcpy.management.Delete(temp_fc)
        arcpy.management.CopyFeatures(current_fc, temp_fc)

        existing_fields = [f.name for f in arcpy.ListFields(temp_fc)]
        for field_name in new_fields:
            if field_name not in existing_fields:
                arcpy.management.AddField(temp_fc, field_name, 'DOUBLE')

        tasks = []
        with arcpy.da.SearchCursor(temp_fc, ['OID@', 'SHAPE@', 'Month']) as cursor:
            for oid, geom, month_val in cursor:
                active_tiff = may_tiff if month_val == 'May_2019' else sep_tiff
                ext = geom.extent
                tasks.append((oid, geom.centroid.X, geom.centroid.Y, ext.XMin, ext.YMin, ext.XMax, ext.YMax, active_tiff, cell_size, '110m' in current_fc))

        results_dict = {}
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for count, result in enumerate(executor.map(process_polygon, tasks), 1):
                oid, metrics = result
                if metrics: results_dict[oid] = metrics

        update_fields = ['OID@'] + new_fields
        with arcpy.da.UpdateCursor(temp_fc, update_fields) as cursor:
            for row in cursor:
                oid = row[0]
                if oid in results_dict:
                    row[1:] = results_dict[oid]
                    cursor.updateRow(row)

        arcpy.conversion.TableToTable(temp_fc, csv_folder, csv_name)
        arcpy.management.Delete(temp_fc)
        print(f"Exported: {csv_name}")

    print('\nALL PROCESSING COMPLETE.')
  
