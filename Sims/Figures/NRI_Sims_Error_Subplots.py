import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import concurrent.futures
import multiprocessing
from scipy.ndimage import distance_transform_edt, gaussian_filter, label

# ====================================================================
# 1. HELPER FUNCTIONS
# ====================================================================
def get_gap_lengths(transect_array, gap_val, p_size):
    valid_transect = transect_array[transect_array != -9999]
    if len(valid_transect) == 0:
        return np.array([])
    padded = np.concatenate(([False], (valid_transect == gap_val), [False]))
    diffs = np.diff(padded.astype(np.int8))
    starts, ends = np.where(diffs == 1)[0], np.where(diffs == -1)[0]
    return (ends - starts) * p_size

# ====================================================================
# 2. ISOLATED WORKER FUNCTION
# ====================================================================
def process_simulation_iteration(task):
    target_bg, plot_radius, hub_radius, cell_size = task
    
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))

    x = np.arange(-plot_radius, plot_radius + (cell_size * 0.1), cell_size)
    y = np.arange(-plot_radius, plot_radius + (cell_size * 0.1), cell_size)
    grid_size = len(x)
    X, Y = np.meshgrid(x, y)
    
    dist_from_center = np.sqrt(X**2 + Y**2)
    valid_mask = dist_from_center <= plot_radius
    total_valid_pixels = np.sum(valid_mask)
    
    angles = [0, 120, 240]
    spoke_length_m = plot_radius - hub_radius
    distances = np.arange(0, spoke_length_m + cell_size, cell_size)
    
    spoke_indices = []
    for a in angles:
        rad = np.radians(a)
        dx = (hub_radius * np.sin(rad)) + distances * np.sin(rad)
        dy = (hub_radius * np.cos(rad)) + distances * np.cos(rad)
        
        cols = np.round((dx - (-plot_radius)) / cell_size).astype(int)
        rows = np.round((plot_radius - dy) / cell_size).astype(int)
        cols = np.clip(cols, 0, grid_size - 1)
        rows = np.clip(rows, 0, grid_size - 1)
        spoke_indices.append((rows, cols))
        
    large_shrub_value = 1
    small_shrub_value = 2
    target_value = 3
    
    total_exhaust_length_m = total_valid_pixels * cell_size * 2
    total_nri_length_m = 3 * spoke_length_m

    # ==========================================
    # Z-INDEXED HYBRID LOGIC 
    # ==========================================
    organic_weight = 0.60
    organic_bg_pct = target_bg * organic_weight
    
    initial_veg_coverage = (100.0 - target_bg) / (100.0 - organic_bg_pct)
    lambda_target = -np.log(1 - initial_veg_coverage) if initial_veg_coverage < 0.99 else 5.0
    
    alpha_large = 0.75
    lambda_large = lambda_target * alpha_large
    lambda_small = lambda_target * (1 - alpha_large)
    
    r_large_bnds = (0.3, 2.0)
    r_small_bnds = (0.10, 0.30)
    
    mean_area_large = np.pi * ((r_large_bnds[0]**2 + r_large_bnds[0]*r_large_bnds[1] + r_large_bnds[1]**2) / 3)
    mean_area_small = np.pi * ((r_small_bnds[0]**2 + r_small_bnds[0]*r_small_bnds[1] + r_small_bnds[1]**2) / 3)
    
    plot_area = np.pi * plot_radius**2
    num_large = int(lambda_large * (plot_area / mean_area_large))
    num_small = int(lambda_small * (plot_area / mean_area_small))
    
    # ==========================================
    # DECOUPLED SPATIAL PROCESSES
    # ==========================================
    cluster_spread = np.interp(target_bg, [0, 100], [8.0, 1.0])
    shrubs_per_cluster = int(np.interp(target_bg, [0, 100], [2, 10]))
    
    num_parents = max(1, num_large // shrubs_per_cluster) if num_large > 0 else 0
    
    large_mask = np.zeros((grid_size, grid_size), dtype=bool)
    small_mask = np.zeros((grid_size, grid_size), dtype=bool)

    shadow_angle = np.radians(135)
    shadow_dx = np.cos(shadow_angle)
    shadow_dy = np.sin(shadow_angle)

    if num_parents > 0:
        gen_radius = plot_radius + r_large_bnds[1] + cluster_spread
        r_parents = gen_radius * np.sqrt(np.random.uniform(0, 1, num_parents))
        theta_parents = np.random.uniform(0, 2 * np.pi, num_parents)
        px_arr = r_parents * np.cos(theta_parents)
        py_arr = r_parents * np.sin(theta_parents)

        if num_large > 0:
            r_large_arr = np.random.uniform(r_large_bnds[0], r_large_bnds[1], num_large)
            for rad in r_large_arr:
                parent_idx = np.random.randint(0, num_parents)
                cx = np.random.normal(px_arr[parent_idx], cluster_spread)
                cy = np.random.normal(py_arr[parent_idx], cluster_spread)
                
                c_min = max(0, int((cx - rad * 1.5 + plot_radius) / cell_size))
                c_max = min(grid_size, int((cx + rad * 1.5 + plot_radius) / cell_size) + 1)
                r_min = max(0, int((cy - rad * 1.5 + plot_radius) / cell_size))
                r_max = min(grid_size, int((cy + rad * 1.5 + plot_radius) / cell_size) + 1)
                if c_min >= c_max or r_min >= r_max: continue
                
                sub_X, sub_Y = X[r_min:r_max, c_min:c_max], Y[r_min:r_max, c_min:c_max]
                dx = sub_X - cx
                dy = sub_Y - cy
                r_grid = np.sqrt(dx**2 + dy**2)
                
                dot_prod = (dx * shadow_dx + dy * shadow_dy) / (r_grid + 1e-5)
                shadow_intensity = np.clip(dot_prod, 0, 1)

                r_core = rad * 0.6
                noise_core = np.random.uniform(-0.15, 0.15, dx.shape)
                core_mask = (r_grid <= r_core * (1.0 + noise_core))
                
                lobe_mask = np.zeros(dx.shape, dtype=bool)
                num_lobes = np.random.randint(3, 6)
                for _ in range(num_lobes):
                    angle = np.random.uniform(0, 2 * np.pi)
                    offset = rad * np.random.uniform(0.3, 0.7)
                    ox = offset * np.cos(angle)
                    oy = offset * np.sin(angle)
                    lobe_rad = rad * np.random.uniform(0.4, 0.8)
                    noise_lobe = np.random.uniform(-0.20, 0.20, dx.shape)
                    lobe_mask |= (np.sqrt((dx - ox)**2 + (dy - oy)**2) <= lobe_rad * (1.0 + noise_lobe))
                
                norm_dist = np.clip((r_grid - r_core) / rad, 0, 1)
                porosity_chance = 0.10 + 0.80 * shadow_intensity * (norm_dist ** 1.5)
                fringe_keep_mask = np.random.rand(*dx.shape) > porosity_chance
                lobe_mask &= fringe_keep_mask
                
                large_mask[r_min:r_max, c_min:c_max] |= (core_mask | lobe_mask)

    # ==========================================
    # THOMAS CLUSTER IMPLEMENTATION FOR SMALL SHRUBS
    # ==========================================
    if num_small > 0:
        shrubs_per_cluster_small = int(np.interp(target_bg, [0, 100], [5, 20]))
        num_parents_small = max(1, num_small // shrubs_per_cluster_small)
        cluster_spread_small = cluster_spread * 0.80 

        if num_parents_small > 0:
            gen_radius_small = plot_radius + r_small_bnds[1] + cluster_spread_small
            r_parents_s = gen_radius_small * np.sqrt(np.random.uniform(0, 1, num_parents_small))
            theta_parents_s = np.random.uniform(0, 2 * np.pi, num_parents_small)
            px_arr_s = r_parents_s * np.cos(theta_parents_s)
            py_arr_s = r_parents_s * np.sin(theta_parents_s)

            r_small_arr = np.random.uniform(r_small_bnds[0], r_small_bnds[1], num_small)
            for rad in r_small_arr:
                parent_idx = np.random.randint(0, num_parents_small)
                cx = np.random.normal(px_arr_s[parent_idx], cluster_spread_small)
                cy = np.random.normal(py_arr_s[parent_idx], cluster_spread_small)
                
                c_min = max(0, int((cx - rad * 1.5 + plot_radius) / cell_size))
                c_max = min(grid_size, int((cx + rad * 1.5 + plot_radius) / cell_size) + 1)
                r_min = max(0, int((cy - rad * 1.5 + plot_radius) / cell_size))
                r_max = min(grid_size, int((cy + rad * 1.5 + plot_radius) / cell_size) + 1)
                if c_min >= c_max or r_min >= r_max: continue
                
                sub_X, sub_Y = X[r_min:r_max, c_min:c_max], Y[r_min:r_max, c_min:c_max]
                dx = sub_X - cx
                dy = sub_Y - cy
                r_grid = np.sqrt(dx**2 + dy**2)
                
                r_core = rad * 0.6
                noise_core = np.random.uniform(-0.15, 0.15, dx.shape)
                core_mask = (r_grid <= r_core * (1.0 + noise_core))
                
                lobe_mask = np.zeros(dx.shape, dtype=bool)
                num_lobes = np.random.randint(2, 5)
                for _ in range(num_lobes):
                    angle = np.random.uniform(0, 2 * np.pi)
                    offset = rad * np.random.uniform(0.3, 0.7)
                    ox = offset * np.cos(angle)
                    oy = offset * np.sin(angle)
                    lobe_rad = rad * np.random.uniform(0.3, 0.8)
                    noise_lobe = np.random.uniform(-0.20, 0.20, dx.shape)
                    lobe_mask |= (np.sqrt((dx - ox)**2 + (dy - oy)**2) <= lobe_rad * (1.0 + noise_lobe))
                
                small_mask[r_min:r_max, c_min:c_max] |= (core_mask | lobe_mask)

    # ==========================================
    # SEPARATE WOODY SHRUBS FROM HERBACEOUS COVER (VECTORIZED)
    # ==========================================
    if target_bg <= 50:
        prob_dark = np.interp(target_bg, [5, 50], [0.05, 0.45])
    else:
        prob_dark = np.interp(target_bg, [50, 100], [0.45, 1.0])

    labeled_large, num_large_features = label(large_mask)
    
    # Vectorized probability assignment to bypass the slow for-loop
    keep_mask = np.random.rand(num_large_features + 1) <= prob_dark
    keep_mask[0] = False  # Always drop the background (label 0)
    
    # Map the booleans directly back to the 2D grid in one fast step
    true_shrub_mask = keep_mask[labeled_large]

    small_mask |= (large_mask & ~true_shrub_mask)
    large_mask = true_shrub_mask

    # ==========================================
    # ASSEMBLE GRID & FRACTAL NOISE PUNCHING
    # ==========================================
    main_array = np.full((grid_size, grid_size), target_value, dtype=int)
    main_array[small_mask] = small_shrub_value
    main_array[large_mask] = large_shrub_value

    target_bare_pixels = int(total_valid_pixels * (target_bg / 100.0))
    current_bare_pixels = np.sum((main_array == target_value) & valid_mask)
    pixels_to_punch = target_bare_pixels - current_bare_pixels
    
    if pixels_to_punch > 0:
        punchable_mask = small_mask & ~large_mask & valid_mask
        punchable_count = np.sum(punchable_mask)
        
        if pixels_to_punch >= punchable_count:
            main_array[punchable_mask] = target_value
        else:
            raw_fine = np.random.rand(grid_size, grid_size)
            raw_coarse = np.random.rand(grid_size, grid_size)
            
            blurred_fine = gaussian_filter(raw_fine, sigma=(1.00 / cell_size))
            blurred_coarse = gaussian_filter(raw_coarse, sigma=(3.00 / cell_size))
            
            fractal_noise = (0.95 * blurred_coarse) + (0.05 * blurred_fine)
            
            punchable_noise = fractal_noise[punchable_mask]
            fraction_to_punch = pixels_to_punch / punchable_count
            threshold = np.percentile(punchable_noise, (1 - fraction_to_punch) * 100)
            
            hole_mask = punchable_mask & (fractal_noise >= threshold)
            main_array[hole_mask] = target_value

    main_array[~valid_mask] = -9999 
    
    is_bare_full = (main_array == target_value)
    bare_pixels = np.sum(is_bare_full[valid_mask])
    true_bg_pct = (bare_pixels / total_valid_pixels) * 100

    true_herb_pct = (np.sum(main_array[valid_mask] == small_shrub_value) / total_valid_pixels) * 100
    true_woody_pct = (np.sum(main_array[valid_mask] == large_shrub_value) / total_valid_pixels) * 100
    
    dist_array = distance_transform_edt(is_bare_full) * cell_size
    valid_full_fetch = dist_array[valid_mask]
    exact_fetch = np.mean(valid_full_fetch) if valid_full_fetch.size > 0 else 0.0
    
    is_bare_masked = np.where(valid_mask, is_bare_full, False)
    
    padded_rows = np.pad(is_bare_masked, ((0, 0), (1, 1)), constant_values=False)
    diff_rows = np.diff(padded_rows.astype(np.int8), axis=1)
    row_starts = np.where(diff_rows == 1)
    row_ends = np.where(diff_rows == -1)
    row_gaps = (row_ends[1] - row_starts[1]) * cell_size
    
    padded_cols = np.pad(is_bare_masked.T, ((0, 0), (1, 1)), constant_values=False)
    diff_cols = np.diff(padded_cols.astype(np.int8), axis=1)
    col_starts = np.where(diff_cols == 1)
    col_ends = np.where(diff_cols == -1)
    col_gaps = (col_ends[1] - col_starts[1]) * cell_size
    
    all_exhaust_gaps = np.concatenate([row_gaps, col_gaps])
    
    ex_0_24 = (np.sum(all_exhaust_gaps[all_exhaust_gaps < 0.25]) / total_exhaust_length_m) * 100
    ex_25_50 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 0.25) & (all_exhaust_gaps <= 0.50)]) / total_exhaust_length_m) * 100
    ex_51_100 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 0.51) & (all_exhaust_gaps <= 1.00)]) / total_exhaust_length_m) * 100
    ex_101_200 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 1.01) & (all_exhaust_gaps <= 2.00)]) / total_exhaust_length_m) * 100
    ex_gt_200 = (np.sum(all_exhaust_gaps[all_exhaust_gaps > 2.00]) / total_exhaust_length_m) * 100
    
    res = {
        'BG_Bin': target_bg,
        'True_BG_Pct': true_bg_pct,
        'True_Herb_Pct': true_herb_pct,
        'True_Woody_Pct': true_woody_pct,
        'Exact_Fetch': exact_fetch,
        'Exact_Gap_0_24': ex_0_24,
        'Exact_Gap_25_50': ex_25_50,
        'Exact_Gap_51_100': ex_51_100,
        'Exact_Gap_101_200': ex_101_200,
        'Exact_Gap_gt_200': ex_gt_200
    }

    all_nri_transects = []
    all_fetch_vals = []
    for rows, cols in spoke_indices:
        all_nri_transects.append(main_array[rows, cols])
        all_fetch_vals.append(dist_array[rows, cols])
        
    bg_intervals = {'0cm': 1, '25cm': 5, '50cm': 10, '100cm': 20, '200cm': 40}
    
    for label_str, step in bg_intervals.items():
        if step == 1:
            all_pixels = np.concatenate([t[:-1] for t in all_nri_transects])
            valid_pixels = all_pixels[all_pixels != -9999] 
            total_nri_pixels = len(valid_pixels)
            sampled_bg = (np.sum(valid_pixels == target_value) / total_nri_pixels) * 100 if total_nri_pixels > 0 else 0.0
            sampled_herb = (np.sum(valid_pixels == small_shrub_value) / total_nri_pixels) * 100 if total_nri_pixels > 0 else 0.0
            sampled_woody = (np.sum(valid_pixels == large_shrub_value) / total_nri_pixels) * 100 if total_nri_pixels > 0 else 0.0
        else:
            sampled_pixels = []
            for t in all_nri_transects:
                sampled_pixels.extend(t[:-1:step])
                
            sampled_pixels_arr = np.array(sampled_pixels)
            valid_sampled_arr = sampled_pixels_arr[sampled_pixels_arr != -9999] 
            total_valid_sampled = len(valid_sampled_arr)
            sampled_bg = (np.sum(valid_sampled_arr == target_value) / total_valid_sampled) * 100 if total_valid_sampled > 0 else 0.0
            sampled_herb = (np.sum(valid_sampled_arr == small_shrub_value) / total_valid_sampled) * 100 if total_valid_sampled > 0 else 0.0
            sampled_woody = (np.sum(valid_sampled_arr == large_shrub_value) / total_valid_sampled) * 100 if total_valid_sampled > 0 else 0.0
            
        res[f'BG_{label_str}_Error'] = sampled_bg - true_bg_pct
        res[f'Herb_{label_str}_Error'] = sampled_herb - true_herb_pct
        res[f'Woody_{label_str}_Error'] = sampled_woody - true_woody_pct
        
    fetch_intervals = {'0cm': 1, '25cm': 5, '50cm': 10, '100cm': 20, '200cm': 40}
    for label_str, step in fetch_intervals.items():
        fetch_sampled = np.concatenate([f[:-1:step] for f in all_fetch_vals]) 
        sampled_mean = np.mean(fetch_sampled[fetch_sampled >= 0]) if fetch_sampled.size > 0 else 0.0
        res[f'Fetch_{label_str}_Error'] = sampled_mean - exact_fetch
    
    all_nri_gaps = np.concatenate([get_gap_lengths(t, target_value, cell_size) for t in all_nri_transects])
    
    valid_nri_pixel_count = np.sum(np.concatenate(all_nri_transects) != -9999)
    adjusted_nri_length_m = valid_nri_pixel_count * cell_size
    
    res['Gap_0_24_0cm_Error'] = ((np.sum(all_nri_gaps[all_nri_gaps < 0.25]) / adjusted_nri_length_m) * 100) - ex_0_24
    res['Gap_25_50_0cm_Error'] = ((np.sum(all_nri_gaps[(all_nri_gaps >= 0.25) & (all_nri_gaps <= 0.50)]) / adjusted_nri_length_m) * 100) - ex_25_50
    res['Gap_51_100_0cm_Error'] = ((np.sum(all_nri_gaps[(all_nri_gaps >= 0.51) & (all_nri_gaps <= 1.00)]) / adjusted_nri_length_m) * 100) - ex_51_100
    res['Gap_101_200_0cm_Error'] = ((np.sum(all_nri_gaps[(all_nri_gaps >= 1.01) & (all_nri_gaps <= 2.00)]) / adjusted_nri_length_m) * 100) - ex_101_200
    res['Gap_gt_200_0cm_Error'] = ((np.sum(all_nri_gaps[all_nri_gaps > 2.00]) / adjusted_nri_length_m) * 100) - ex_gt_200
    
    return res

# ====================================================================
# 3. MAIN EXECUTION BLOCK 
# ====================================================================
if __name__ == '__main__':
    num_iterations = 190000 
    plot_radius = 55
    hub_radius = 5
    cell_size = 0.05
    
    target_bgs = np.arange(5, 100, 5)
    iters_per_bin = max(1, num_iterations // len(target_bgs))
    
    tasks = []
    for target_bg in target_bgs:
        for _ in range(iters_per_bin):
            tasks.append((target_bg, plot_radius, hub_radius, cell_size))
            
    total_tasks = len(tasks)
    cores = max(1, multiprocessing.cpu_count() - 3)
    
    print(f"{'='*50}")
    print(f"Starting Multiprocessing Simulation")
    print(f"Target Bins: {len(target_bgs)}")
    print(f"Iterations per Bin: {iters_per_bin}")
    print(f"Total Iterations: {total_tasks}")
    print(f"Detected CPU Cores: {cores}")
    print(f"{'='*50}")

    results = []
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        for count, result in enumerate(executor.map(process_simulation_iteration, tasks), 1):
            results.append(result)
            if count % 100 == 0 or count == total_tasks:
                print(f"  -> Processed {count} / {total_tasks} simulations...")

    print("\nProcessing complete. Generating statistics and plots...")

    df = pd.DataFrame(results)
    
    def calc_mre(error, exact):
        m_exact = np.mean(exact)
        return 0.0 if m_exact == 0 else (np.mean(np.abs(error)) / m_exact) * 100

    def aggregate_metrics(x):
        d = {
            'True_BG_Mean': np.mean(x['True_BG_Pct']),
            'True_Herb_Mean': np.mean(x['True_Herb_Pct']),
            'True_Woody_Mean': np.mean(x['True_Woody_Pct']),
            'Exact_Fetch_Mean': np.mean(x['Exact_Fetch'])
        }
        for gap in ['Gap_0_24', 'Gap_25_50', 'Gap_51_100', 'Gap_101_200', 'Gap_gt_200']:
            d[f'Exact_{gap}_Mean'] = np.mean(x[f'Exact_{gap}'])
        
        for scale in ['0cm', '25cm', '50cm', '100cm', '200cm']:
            d[f'BG_{scale}_Val'] = np.mean(x[f'BG_{scale}_Error'] + x['True_BG_Pct'])
            d[f'BG_{scale}_MAE'] = np.mean(np.abs(x[f'BG_{scale}_Error']))
            d[f'BG_{scale}_MRE'] = calc_mre(x[f'BG_{scale}_Error'], x['True_BG_Pct'])
            d[f'BG_{scale}_Bias'] = np.mean(x[f'BG_{scale}_Error'])

            d[f'Herb_{scale}_Val'] = np.mean(x[f'Herb_{scale}_Error'] + x['True_Herb_Pct'])
            d[f'Herb_{scale}_MAE'] = np.mean(np.abs(x[f'Herb_{scale}_Error']))
            d[f'Herb_{scale}_MRE'] = calc_mre(x[f'Herb_{scale}_Error'], x['True_Herb_Pct'])
            d[f'Herb_{scale}_Bias'] = np.mean(x[f'Herb_{scale}_Error'])

            d[f'Woody_{scale}_Val'] = np.mean(x[f'Woody_{scale}_Error'] + x['True_Woody_Pct'])
            d[f'Woody_{scale}_MAE'] = np.mean(np.abs(x[f'Woody_{scale}_Error']))
            d[f'Woody_{scale}_MRE'] = calc_mre(x[f'Woody_{scale}_Error'], x['True_Woody_Pct'])
            d[f'Woody_{scale}_Bias'] = np.mean(x[f'Woody_{scale}_Error'])
            
        for scale in ['0cm', '25cm', '50cm', '100cm', '200cm']:
            d[f'Fetch_{scale}_Val'] = np.mean(x[f'Fetch_{scale}_Error'] + x['Exact_Fetch'])
            d[f'Fetch_{scale}_MAE'] = np.mean(np.abs(x[f'Fetch_{scale}_Error']))
            d[f'Fetch_{scale}_MRE'] = calc_mre(x[f'Fetch_{scale}_Error'], x['Exact_Fetch'])
            d[f'Fetch_{scale}_Bias'] = np.mean(x[f'Fetch_{scale}_Error'])
            
        for gap in ['Gap_0_24', 'Gap_25_50', 'Gap_51_100', 'Gap_101_200', 'Gap_gt_200']:
            d[f'{gap}_0cm_Val'] = np.mean(x[f'{gap}_0cm_Error'] + x[f'Exact_{gap}'])
            d[f'{gap}_0cm_MAE'] = np.mean(np.abs(x[f'{gap}_0cm_Error']))
            d[f'{gap}_0cm_MRE'] = calc_mre(x[f'{gap}_0cm_Error'], x[f'Exact_{gap}'])
            d[f'{gap}_0cm_Bias'] = np.mean(x[f'{gap}_0cm_Error'])
            
        return pd.Series(d)

    mae_df = df.groupby('BG_Bin').apply(aggregate_metrics).reset_index()
    mae_df = mae_df.sort_values('True_BG_Mean').reset_index(drop=True)

    scale_colors = {
        '0cm': 'black',
        '25cm': 'forestgreen',
        '50cm': 'dodgerblue',
        '100cm': 'darkorange',
        '200cm': 'crimson'
    }

    # ====================================================================
    # PLOTTING - FIGURE 1: Original 4x9 Configuration
    # ====================================================================
    plot_config = [
        ('BG', 'Total Bare Ground (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
        ('Fetch', 'Mean Fetch (m)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
        ('Gap_0_24', 'Canopy Gap 0-24cm (%)', ['0cm'], 'steelblue'),
        ('Gap_25_50', 'Canopy Gap 25-50cm (%)', ['0cm'], 'cadetblue'),
        ('Gap_51_100', 'Canopy Gap 51-100cm (%)', ['0cm'], 'mediumseagreen'),
        ('Gap_101_200', 'Canopy Gap 101-200cm (%)', ['0cm'], 'darkorange'),
        ('Gap_gt_200', 'Canopy Gap >200cm (%)', ['0cm'], 'firebrick'),
        ('Herb', 'Herbaceous Cover (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
        ('Woody', 'Woody Cover (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None)
    ]

    fig, axes = plt.subplots(4, 9, figsize=(36, 16), constrained_layout=True)
    fig.suptitle("Hybrid Simulation Metrics (Values, MAE, MRE, Bias) Across Bare Ground Gradient", fontsize=20, weight='bold')
    
    for col, (prefix, title, scales, col_color) in enumerate(plot_config):
        ax_val = axes[0, col]
        ax_mae = axes[1, col]
        ax_mre = axes[2, col]
        ax_bias = axes[3, col]
        
        if prefix == 'BG':
            ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_BG_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
        elif prefix == 'Fetch':
            ax_val.plot(mae_df['True_BG_Mean'], mae_df['Exact_Fetch_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
        elif prefix.startswith('Gap'):
            ax_val.plot(mae_df['True_BG_Mean'], mae_df[f'Exact_{prefix}_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
        elif prefix == 'Herb':
            ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_Herb_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
        elif prefix == 'Woody':
            ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_Woody_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)

        for scale in scales:
            var_base = f'{prefix}_{scale}'
            color = scale_colors[scale] if len(scales) > 1 else col_color
            label_name = f'{scale} Point' if scale != '0cm' else '0cm Continuous'
            line_kws = {'marker': 'o', 'color': color, 'linestyle': '-', 'linewidth': 2.5, 'markersize': 7}
            
            ax_val.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Val'], label=label_name, **line_kws)
            ax_mae.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MAE'], **line_kws)
            ax_mre.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MRE'], **line_kws)
            ax_bias.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Bias'], **line_kws)
            
        for ax in [ax_val, ax_mae, ax_mre, ax_bias]:
            ax.set_xlim(0, 100)
            ax.set_xticks(np.arange(0, 81, 20)) 
            
        ax_val.set_title(title, fontsize=15, pad=12)
        ax_val.grid(True, alpha=0.3)
        ax_val.legend(loc='upper left', fontsize=10)
        
        if col == 0: ax_val.set_ylabel("Sampled Value", fontsize=13)
        
        ax_mae.grid(True, alpha=0.3)
        if col == 0: ax_mae.set_ylabel("MAE", fontsize=13)
            
        ax_mre.grid(True, alpha=0.3)
        if col == 0: ax_mre.set_ylabel("MRE (%)", fontsize=13)
            
        ax_bias.axhline(0, color='gray', linestyle='--', linewidth=1.5)
        ax_bias.grid(True, alpha=0.3)
        ax_bias.set_xlabel("True Bare Ground (%)", fontsize=12)
        if col == 0: ax_bias.set_ylabel("Mean Bias", fontsize=13)

    # ====================================================================
    # PLOTTING - FIGURE 2: COVER METRICS (BG, HERB, WOODY) - 3 Rows x 4 Cols
    # ====================================================================
    fig2_config = [
        ('BG', 'Total Bare Ground (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
        ('Herb', 'Herbaceous Cover (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
        ('Woody', 'Woody Cover (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None)
    ]

    fig2, axes2 = plt.subplots(3, 4, figsize=(16, 12), constrained_layout=True)
    fig2.suptitle("Hybrid Simulation Cover Metrics Across Bare Ground Gradient", fontsize=20, weight='bold')

    for row, (prefix, title, scales, col_color) in enumerate(fig2_config):
        ax_val = axes2[row, 0]
        ax_mae = axes2[row, 1]
        ax_mre = axes2[row, 2]
        ax_bias = axes2[row, 3]

        if prefix == 'BG':
            ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_BG_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
        elif prefix == 'Herb':
            ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_Herb_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
        elif prefix == 'Woody':
            ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_Woody_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)

        for scale in scales:
            var_base = f'{prefix}_{scale}'
            color = scale_colors[scale] if len(scales) > 1 else col_color
            label_name = f'{scale} Point' if scale != '0cm' else '0cm Continuous'
            line_kws = {'marker': 'o', 'color': color, 'linestyle': '-', 'linewidth': 2.5, 'markersize': 7}
            
            ax_val.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Val'], label=label_name, **line_kws)
            ax_mae.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MAE'], **line_kws)
            ax_mre.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MRE'], **line_kws)
            ax_bias.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Bias'], **line_kws)

        for ax in [ax_val, ax_mae, ax_mre, ax_bias]:
            ax.set_xlim(0, 100)
            ax.set_xticks(np.arange(0, 81, 20)) 
            
        ax_val.set_ylabel(title, fontsize=13, weight='bold') 
        ax_val.grid(True, alpha=0.3)
        
        if row == 0:
            ax_val.legend(loc='upper left', fontsize=10)

        ax_mae.grid(True, alpha=0.3)
        ax_mre.grid(True, alpha=0.3)

        ax_bias.axhline(0, color='gray', linestyle='--', linewidth=1.5)
        ax_bias.grid(True, alpha=0.3)

        if row == 0:
            ax_val.set_title("Sampled Value", fontsize=15, pad=12)
            ax_mae.set_title("MAE", fontsize=15, pad=12)
            ax_mre.set_title("MRE (%)", fontsize=15, pad=12)
            ax_bias.set_title("Mean Bias", fontsize=15, pad=12)

        if row == 2:
            for ax in [ax_val, ax_mae, ax_mre, ax_bias]:
                ax.set_xlabel("True Bare Ground (%)", fontsize=12)

    # ====================================================================
    # PLOTTING - FIGURE 3: SPATIAL METRICS (FETCH, GAPS) - 6 Rows x 4 Cols
    # ====================================================================
    fig3_config = [
        ('Fetch', 'Mean Fetch (m)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
        ('Gap_0_24', 'Canopy Gap\n0-24cm (%)', ['0cm'], 'steelblue'),
        ('Gap_25_50', 'Canopy Gap\n25-50cm (%)', ['0cm'], 'cadetblue'),
        ('Gap_51_100', 'Canopy Gap\n51-100cm (%)', ['0cm'], 'mediumseagreen'),
        ('Gap_101_200', 'Canopy Gap\n101-200cm (%)', ['0cm'], 'darkorange'),
        ('Gap_gt_200', 'Canopy Gap\n>200cm (%)', ['0cm'], 'firebrick')
    ]

    fig3, axes3 = plt.subplots(6, 4, figsize=(16, 24), constrained_layout=True)
    fig3.suptitle("Hybrid Simulation Spatial Metrics Across Bare Ground Gradient", fontsize=20, weight='bold')

    for row, (prefix, title, scales, col_color) in enumerate(fig3_config):
        ax_val = axes3[row, 0]
        ax_mae = axes3[row, 1]
        ax_mre = axes3[row, 2]
        ax_bias = axes3[row, 3]

        if prefix == 'Fetch':
            ax_val.plot(mae_df['True_BG_Mean'], mae_df['Exact_Fetch_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
        elif prefix.startswith('Gap'):
            ax_val.plot(mae_df['True_BG_Mean'], mae_df[f'Exact_{prefix}_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)

        for scale in scales:
            var_base = f'{prefix}_{scale}'
            color = scale_colors[scale] if len(scales) > 1 else col_color
            label_name = f'{scale} Point' if scale != '0cm' else '0cm Continuous'
            line_kws = {'marker': 'o', 'color': color, 'linestyle': '-', 'linewidth': 2.5, 'markersize': 7}
            
            ax_val.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Val'], label=label_name, **line_kws)
            ax_mae.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MAE'], **line_kws)
            ax_mre.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MRE'], **line_kws)
            ax_bias.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Bias'], **line_kws)

        for ax in [ax_val, ax_mae, ax_mre, ax_bias]:
            ax.set_xlim(0, 100)
            ax.set_xticks(np.arange(0, 81, 20)) 
            
        ax_val.set_ylabel(title, fontsize=13, weight='bold') 
        ax_val.grid(True, alpha=0.3)
        
        # Determine legend placement based on whether the data trends up or down
        last_var_base = f'{prefix}_{scales[-1]}'
        start_val = mae_df[f'{last_var_base}_Val'].iloc[0]
        end_val = mae_df[f'{last_var_base}_Val'].iloc[-1]
        leg_loc = 'upper left' if end_val >= start_val else 'upper right'
        
        ax_val.legend(loc=leg_loc, fontsize=10)

        ax_mae.grid(True, alpha=0.3)
        ax_mre.grid(True, alpha=0.3)

        ax_bias.axhline(0, color='gray', linestyle='--', linewidth=1.5)
        ax_bias.grid(True, alpha=0.3)

        if row == 0:
            ax_val.set_title("Sampled Value", fontsize=15, pad=12)
            ax_mae.set_title("MAE", fontsize=15, pad=12)
            ax_mre.set_title("MRE (%)", fontsize=15, pad=12)
            ax_bias.set_title("Mean Bias", fontsize=15, pad=12)

        if row == 5:
            for ax in [ax_val, ax_mae, ax_mre, ax_bias]:
                ax.set_xlabel("True Bare Ground (%)", fontsize=12)

    # ====================================================================
    # SAVING AND SHOWING PLOTS
    # ====================================================================
    script_dir = os.path.dirname(os.path.abspath(__file__))
    img_path1 = os.path.join(script_dir, 'Hybrid_Simulation_Metrics_MultiScale.svg')
    img_path2 = os.path.join(script_dir, 'Hybrid_Simulation_Cover_Metrics.svg')
    img_path3 = os.path.join(script_dir, 'Hybrid_Simulation_Spatial_Metrics_Rows.svg')
    csv_path = os.path.join(script_dir, 'Hybrid_Simulation_Metrics_MultiScale.csv')
    
    fig.savefig(img_path1, format='svg', dpi=300, bbox_inches='tight')
    fig2.savefig(img_path2, format='svg', dpi=300, bbox_inches='tight')
    fig3.savefig(img_path3, format='svg', dpi=300, bbox_inches='tight')
    
    plt.show()
    
