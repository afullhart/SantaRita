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

# --- INPUT FOOTPRINTS ---
in_nri_plots = os.path.join(out_gdb, 'SRER_NRI_Plots_110m')
in_grid_10m = os.path.join(out_gdb, 'SRER_Grid_10m')
in_grid_30m = os.path.join(out_gdb, 'SRER_Grid_30m')

# --- OUTPUT FEATURE CLASSES ---
out_nri_plots = os.path.join(out_gdb, 'SRER_NRI_Plots_110m_RandomLines')
out_grid_10m = os.path.join(out_gdb, 'SRER_Grid_10m_RandomLines')
out_grid_30m = os.path.join(out_gdb, 'SRER_Grid_30m_RandomLines')

may_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_May_2019_UTM12N_Mosaic.tif'
sep_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_Sep_2019_UTM12N_Mosaic.tif'

# --- SAMPLING DESIGN INCREMENTS ---
PT_INCS = [3, 5, 7, 10, 15, 20, 30, 50, 70, 100, 200, 300, 500, 700, 1000, 5000, 10000]
LN_INCS = [3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 1500, 3000]

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
        
        # Convert shadow (4) to herb cover (1)
        main_array[main_array == 4] = herb_value
        
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

        if is_circle:
            herb_px = np.sum((main_array == herb_value) & mask)
            woody_px = np.sum((main_array == woody_value) & mask)
        else:
            herb_px = np.sum(main_array == herb_value)
            woody_px = np.sum(main_array == woody_value)

        herb_pct_exact = (float(herb_px) / float(total_pixels)) * 100
        woody_pct_exact = (float(woody_px) / float(total_pixels)) * 100
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
            f_101_200= (np.sum(all_gap_lengths[(all_gap_lengths >= 1.01) & (all_gaps_lengths <= 2.00)]) / total_transect_length_m) * 100
            f_gt_200 = (np.sum(all_gap_lengths[(all_gap_lengths > 2.00)]) / total_transect_length_m) * 100
        else:
            f_0_24 = f_25_50 = f_51_100 = f_101_200 = f_gt_200 = 0.0

        output_metrics = [
            bgr_percent, lpi_percent, mean_fetch_exact, hw_ratio_exact, 
            herb_pct_exact, woody_pct_exact,
            f_0_24, f_25_50, f_51_100, f_101_200, f_gt_200
        ]

        # ==========================================
        # 2. INDEPENDENT TRANSECT HELPER FUNCTION
        # ==========================================
        valid_y, valid_x = np.where(mask)
        num_valid_starts = len(valid_y)
        
        step = int(round(1.0 / c_size))
        if step < 1: 
            step = 1

        def get_independent_transect(target_px):
            """Generates an independent, randomized virtual transect of target_px length."""
            collected_px = 0
            chunks_main = []
            chunks_dist = []
            
            while collected_px < target_px and num_valid_starts > 0:
                r_idx = np.random.randint(0, num_valid_starts)
                sy, sx = valid_y[r_idx], valid_x[r_idx]
                direction = np.random.randint(0, 4)
                rem_px = target_px - collected_px
                
                if direction == 0: 
                    c_main = main_array[sy, sx : min(sx + rem_px, ncols)]
                    c_dist = dist_array[sy, sx : min(sx + rem_px, ncols)]
                    c_mask = mask[sy, sx : min(sx + rem_px, ncols)]
                elif direction == 1: 
                    c_main = main_array[sy : min(sy + rem_px, nrows), sx]
                    c_dist = dist_array[sy : min(sy + rem_px, nrows), sx]
                    c_mask = mask[sy : min(sy + rem_px, nrows), sx]
                elif direction == 2: 
                    c_main = main_array[sy, max(0, sx - rem_px + 1) : sx + 1][::-1]
                    c_dist = dist_array[sy, max(0, sx - rem_px + 1) : sx + 1][::-1]
                    c_mask = mask[sy, max(0, sx - rem_px + 1) : sx + 1][::-1]
                else: 
                    c_main = main_array[max(0, sy - rem_px + 1) : sy + 1, sx][::-1]
                    c_dist = dist_array[max(0, sy - rem_px + 1) : sy + 1, sx][::-1]
                    c_mask = mask[max(0, sy - rem_px + 1) : sy + 1, sx][::-1]
                    
                invalid_idx = np.where(~c_mask)[0]
                if len(invalid_idx) > 0:
                    c_main = c_main[:invalid_idx[0]]
                    c_dist = c_dist[:invalid_idx[0]]
                    
                if len(c_main) > 0:
                    chunks_main.append(c_main)
                    chunks_dist.append(c_dist)
                    collected_px += len(c_main)
                    
            return chunks_main, chunks_dist

        # ==========================================
        # 3. LINE-BASED 1-METER POINT SAMPLING
        # ==========================================
        for n_pts in PT_INCS:
            req_px = n_pts * step
            
            chunks_main, chunks_dist = get_independent_transect(req_px)
            
            if req_px == 0 or len(chunks_main) == 0:
                output_metrics.extend([0.0, 0.0, None, 0.0, 0.0])
                continue
                
            full_main = np.concatenate(chunks_main)
            full_dist = np.concatenate(chunks_dist)
            
            # Step through the independent line array at exactly 1m intervals
            idx = np.arange(0, len(full_main), step)
            sampled_main = full_main[idx]
            sampled_dist = full_dist[idx]
            
            actual_pts = len(sampled_main)
            if actual_pts == 0:
                output_metrics.extend([0.0, 0.0, None, 0.0, 0.0])
                continue
                
            smpl_bgr = (np.sum(sampled_main == target_value) / actual_pts) * 100
            smpl_fetch = np.mean(sampled_dist)
            
            smpl_herb_count = np.sum(sampled_main == herb_value)
            smpl_woody_count = np.sum(sampled_main == woody_value)
            
            smpl_hw_ratio = (float(smpl_herb_count) / float(smpl_woody_count)) if smpl_woody_count > 0 else None
            smpl_herb_pct = (float(smpl_herb_count) / actual_pts) * 100
            smpl_woody_pct = (float(smpl_woody_count) / actual_pts) * 100
            
            output_metrics.extend([smpl_bgr, smpl_fetch, smpl_hw_ratio, smpl_herb_pct, smpl_woody_pct])

        # ==========================================
        # 4. CONTINUOUS LINE-INTERCEPT (CGF)
        # ==========================================
        for L_meters in LN_INCS:
            req_px = int(round(L_meters / c_size))
            
            chunks_main, _ = get_independent_transect(req_px)
            
            all_gaps = []
            actual_L = 0
            
            # Loop over chunks, padding each individually to strictly prevent bridging
            for chunk in chunks_main:
                is_gap = (chunk == target_value)
                padded = np.concatenate(([False], is_gap, [False]))
                diffs = np.diff(padded.astype(int))
                starts = np.where(diffs == 1)[0]
                ends = np.where(diffs == -1)[0]
                all_gaps.extend((ends - starts) * c_size)
                actual_L += len(chunk) * c_size
                
            all_gaps = np.array(all_gaps)
            
            if actual_L > 0:
                s_0_24   = (np.sum(all_gaps[(all_gaps < 0.25)]) / actual_L) * 100
                s_25_50  = (np.sum(all_gaps[(all_gaps >= 0.25) & (all_gaps <= 0.50)]) / actual_L) * 100
                s_51_100 = (np.sum(all_gaps[(all_gaps >= 0.51) & (all_gaps <= 1.00)]) / actual_L) * 100
                s_101_200= (np.sum(all_gaps[(all_gaps >= 1.01) & (all_gaps <= 2.00)]) / actual_L) * 100
                s_gt_200 = (np.sum(all_gaps[(all_gaps > 2.00)]) / actual_L) * 100
                output_metrics.extend([s_0_24, s_25_50, s_51_100, s_101_200, s_gt_200])
            else:
                output_metrics.extend([0.0, 0.0, 0.0, 0.0, 0.0])

        return oid, season, output_metrics
        
    except Exception as e:
        return oid, season, None

# ====================================================================
# EXECUTION LOGIC
# ====================================================================
def calculate_metrics_for_fc(in_fc, out_fc, cell_size, is_circle=False):
    fc_name = os.path.basename(out_fc)
    print(f'\n--- Processing Seasonal Metrics & Independent Line Undersampling for {fc_name} ---')

    shapes_dict = {}
    tasks = []
    
    with arcpy.da.SearchCursor(in_fc, ['OID@', 'SHAPE@']) as cursor:
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

    print('  -> Constructing final feature class with all line-sampled fields...')
    
    sr = arcpy.Describe(in_fc).spatialReference
    temp_out = r'memory\temp_metrics_fc'
    
    if arcpy.Exists(temp_out):
        arcpy.management.Delete(temp_out)
            
    arcpy.management.CreateFeatureclass('memory', 'temp_metrics_fc', 'POLYGON', spatial_reference=sr)
    
    arcpy.management.AddField(temp_out, 'Season', 'TEXT')
    
    new_fields = [
        'BGR_Exact', 'LPI_Exact', 'Fetch_Exact', 'HW_Ratio_Exact', 'Herb_Pct_Exact', 'Woody_Pct_Exact',
        'Gap_0_24_Exact', 'Gap_25_50_Exact', 'Gap_51_100_Exact', 'Gap_101_200_Exact', 'Gap_gt_200_Exact'
    ]
    
    for pt in PT_INCS:
        new_fields.extend([f'BGR_pt_{pt}', f'Fetch_pt_{pt}', f'HW_Ratio_pt_{pt}', f'Herb_Pct_pt_{pt}', f'Woody_Pct_pt_{pt}'])
        
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

    arcpy.management.CopyFeatures(temp_out, out_fc)
    arcpy.management.Delete(temp_out)
    
    print(f'  -> Successfully created and updated {fc_name}!')

if __name__ == '__main__':
    arcpy.ResetEnvironments() 
    arcpy.env.overwriteOutput = True

    print('Reading base raster properties...')
    pixel_size = float(arcpy.management.GetRasterProperties(may_tiff, 'CELLSIZEX').getOutput(0))
    print(f'Detected pixel size: {pixel_size}m')

    calculate_metrics_for_fc(in_grid_10m, out_grid_10m, pixel_size, is_circle=False)
    calculate_metrics_for_fc(in_grid_30m, out_grid_30m, pixel_size, is_circle=False)
    calculate_metrics_for_fc(in_nri_plots, out_nri_plots, pixel_size, is_circle=True)

    print('\nAll Extracted Metrics & Independent Line Processing Complete!')

    csv_out_folder = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data'
    print('\n--- Exporting Tables to CSV ---')

    export_tasks = {
        out_grid_10m: 'SRER_Grid_10m_RandomLines.csv',
        out_grid_30m: 'SRER_Grid_30m_RandomLines.csv',
        out_nri_plots: 'SRER_NRI_Plots_110m_RandomLines.csv'
    }

    for fc_path, csv_name in export_tasks.items():
        out_csv_path = os.path.join(csv_out_folder, csv_name)
        
        if arcpy.Exists(out_csv_path):
            arcpy.management.Delete(out_csv_path)
            
        print(f'  -> Exporting {os.path.basename(fc_path)} to {csv_name}...')
        
        arcpy.conversion.ExportTable(
            in_table=fc_path,
            out_table=out_csv_path
        )

    print(f'\nSuccess! All metrics exported to: {csv_out_folder}')
    
