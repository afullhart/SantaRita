import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import concurrent.futures
import multiprocessing
from scipy.ndimage import distance_transform_edt, gaussian_filter

# ====================================================================
# 1. HELPER FUNCTIONS
# ====================================================================
def get_gap_lengths(transect_array, gap_val, p_size):
    padded = np.concatenate(([False], (transect_array == gap_val), [False]))
    diffs = np.diff(padded.astype(int))
    starts, ends = np.where(diffs == 1)[0], np.where(diffs == -1)[0]
    return (ends - starts) * p_size

# ====================================================================
# 2. ISOLATED WORKER FUNCTION
# ====================================================================
def process_simulation_iteration(task):
    target_bg, plot_radius, hub_radius, cell_size = task
    
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))

    # EXACT grid coordinates to prevent sub-pixel aliasing drift along the sampling rays
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
        
    target_value = 3
    shrub_value = 2
    total_exhaust_length_m = total_valid_pixels * cell_size * 2
    total_nri_length_m = 3 * spoke_length_m

    # ==========================================
    # Z-INDEXED HYBRID LOGIC 
    # ==========================================
    organic_weight = 0.60
    organic_bg_pct = target_bg * organic_weight
    
    initial_veg_coverage = (100.0 - target_bg) / (100.0 - organic_bg_pct)
    lambda_target = -np.log(1 - initial_veg_coverage) if initial_veg_coverage < 0.99 else 5.0
    
    # Locked at 0.25 to ensure small vegetation ALWAYS persists across the gradient
    alpha_large = 0.25
    lambda_large = lambda_target * alpha_large
    lambda_small = lambda_target * (1 - alpha_large)
    
    r_large_bnds = (0.5, 3.0)
    r_small_bnds = (0.10, 0.30)
    
    mean_area_large = np.pi * ((r_large_bnds[0]**2 + r_large_bnds[0]*r_large_bnds[1] + r_large_bnds[1]**2) / 3)
    mean_area_small = np.pi * ((r_small_bnds[0]**2 + r_small_bnds[0]*r_small_bnds[1] + r_small_bnds[1]**2) / 3)
    
    plot_area = np.pi * plot_radius**2
    num_large = int(lambda_large * (plot_area / mean_area_large))
    num_small = int(lambda_small * (plot_area / mean_area_small))
    
    large_mask = np.zeros((grid_size, grid_size), dtype=bool)
    if num_large > 0:
        gen_radius_large = plot_radius + r_large_bnds[1]
        r_large_arr = np.random.uniform(r_large_bnds[0], r_large_bnds[1], num_large)
        r_dist_L = gen_radius_large * np.sqrt(np.random.uniform(0, 1, num_large))
        theta_arr_L = np.random.uniform(0, 2 * np.pi, num_large)
        cx_arr_L = r_dist_L * np.cos(theta_arr_L)
        cy_arr_L = r_dist_L * np.sin(theta_arr_L)
        
        for cx, cy, rad in zip(cx_arr_L, cy_arr_L, r_large_arr):
            c_min = max(0, int((cx - rad + plot_radius) / cell_size))
            c_max = min(grid_size, int((cx + rad + plot_radius) / cell_size) + 1)
            r_min = max(0, int((cy - rad + plot_radius) / cell_size))
            r_max = min(grid_size, int((cy + rad + plot_radius) / cell_size) + 1)
            if c_min >= c_max or r_min >= r_max: continue
            
            sub_X, sub_Y = X[r_min:r_max, c_min:c_max], Y[r_min:r_max, c_min:c_max]
            inside = ((sub_X - cx)**2 + (sub_Y - cy)**2) <= rad**2
            large_mask[r_min:r_max, c_min:c_max] |= inside

    small_mask = np.zeros((grid_size, grid_size), dtype=bool)
    if num_small > 0:
        gen_radius_small = plot_radius + r_small_bnds[1]
        r_small_arr = np.random.uniform(r_small_bnds[0], r_small_bnds[1], num_small)
        r_dist_S = gen_radius_small * np.sqrt(np.random.uniform(0, 1, num_small))
        theta_arr_S = np.random.uniform(0, 2 * np.pi, num_small)
        cx_arr_S = r_dist_S * np.cos(theta_arr_S)
        cy_arr_S = r_dist_S * np.sin(theta_arr_S)
        
        for cx, cy, rad in zip(cx_arr_S, cy_arr_S, r_small_arr):
            c_min = max(0, int((cx - rad + plot_radius) / cell_size))
            c_max = min(grid_size, int((cx + rad + plot_radius) / cell_size) + 1)
            r_min = max(0, int((cy - rad + plot_radius) / cell_size))
            r_max = min(grid_size, int((cy + rad + plot_radius) / cell_size) + 1)
            if c_min >= c_max or r_min >= r_max: continue
            
            sub_X, sub_Y = X[r_min:r_max, c_min:c_max], Y[r_min:r_max, c_min:c_max]
            inside = ((sub_X - cx)**2 + (sub_Y - cy)**2) <= rad**2
            small_mask[r_min:r_max, c_min:c_max] |= inside

    main_array = np.full((grid_size, grid_size), target_value, dtype=int)
    main_array[small_mask] = shrub_value
    main_array[large_mask] = shrub_value

    target_bare_pixels = int(total_valid_pixels * (target_bg / 100.0))
    current_bare_pixels = np.sum((main_array == target_value) & valid_mask)
    pixels_to_punch = target_bare_pixels - current_bare_pixels
    
    if pixels_to_punch > 0:
        punchable_mask = small_mask & ~large_mask & valid_mask
        punchable_count = np.sum(punchable_mask)
        
        if pixels_to_punch >= punchable_count:
            main_array[punchable_mask] = target_value
        else:
            # MULTI-SCALE FRACTAL NOISE
            # Combines 3.00m macro-voids with 0.50m jagged micro-edges
            raw_fine = np.random.rand(grid_size, grid_size)
            raw_coarse = np.random.rand(grid_size, grid_size)
            
            blurred_fine = gaussian_filter(raw_fine, sigma=(0.50 / cell_size))
            blurred_coarse = gaussian_filter(raw_coarse, sigma=(3.00 / cell_size))
            
            fractal_noise = (0.85 * blurred_coarse) + (0.15 * blurred_fine)
            
            punchable_noise = fractal_noise[punchable_mask]
            fraction_to_punch = pixels_to_punch / punchable_count
            threshold = np.percentile(punchable_noise, (1 - fraction_to_punch) * 100)
            
            hole_mask = punchable_mask & (fractal_noise >= threshold)
            main_array[hole_mask] = target_value

    main_array[~valid_mask] = -9999 
    
    is_bare_full = (main_array == target_value)
    bare_pixels = np.sum(is_bare_full[valid_mask])
    true_bg_pct = (bare_pixels / total_valid_pixels) * 100
    
    dist_array = distance_transform_edt(is_bare_full) * cell_size
    valid_full_fetch = dist_array[valid_mask]
    exact_fetch = np.mean(valid_full_fetch) if valid_full_fetch.size > 0 else 0.0
    
    is_bare_masked = np.where(valid_mask, is_bare_full, False)
    
    padded_rows = np.pad(is_bare_masked, ((0, 0), (1, 1)), constant_values=False)
    diff_rows = np.diff(padded_rows.astype(int), axis=1)
    row_starts = np.where(diff_rows == 1)
    row_ends = np.where(diff_rows == -1)
    row_gaps = (row_ends[1] - row_starts[1]) * cell_size
    
    padded_cols = np.pad(is_bare_masked, ((1, 1), (0, 0)), constant_values=False)
    diff_cols = np.diff(padded_cols.astype(int), axis=0)
    col_starts = np.where(diff_cols == 1)
    col_ends = np.where(diff_cols == -1)
    col_gaps = (col_ends[0] - col_starts[0]) * cell_size
    
    all_exhaust_gaps = np.concatenate([row_gaps, col_gaps])
    
    ex_0_24 = (np.sum(all_exhaust_gaps[all_exhaust_gaps < 0.25]) / total_exhaust_length_m) * 100
    ex_25_50 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 0.25) & (all_exhaust_gaps <= 0.50)]) / total_exhaust_length_m) * 100
    ex_51_100 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 0.51) & (all_exhaust_gaps <= 1.00)]) / total_exhaust_length_m) * 100
    ex_101_200 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 1.01) & (all_exhaust_gaps <= 2.00)]) / total_exhaust_length_m) * 100
    ex_gt_200 = (np.sum(all_exhaust_gaps[all_exhaust_gaps > 2.00]) / total_exhaust_length_m) * 100
    
    res = {
        'BG_Bin': target_bg,
        'True_BG_Pct': true_bg_pct,
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
    
    # ========================================================================
    # EXPLICIT INTERVAL ASSUMPTION & FOOTPRINT INFLATION
    # ========================================================================
    for label, step in bg_intervals.items():
        if step == 1:
            total_nri_pixels = sum(len(t) for t in all_nri_transects)
            nri_bare_pixels = sum(np.sum(t == target_value) for t in all_nri_transects)
            sampled_bg = (nri_bare_pixels / total_nri_pixels) * 100 if total_nri_pixels > 0 else 0.0
        else:
            # CONSTANT FOOTPRINT: Set to a static 1-pixel buffer for all scales.
            # This applies an identical measurement bias window to the 25cm, 50cm, 
            # 100cm, and 200cm points, tightly clustering the bias lines together!
            footprint_radius = 1 
            stretch_value = target_value if true_bg_pct <= 50 else shrub_value
            
            sampled_pixels_inflated = []
            for t in all_nri_transects:
                for i in range(0, len(t), step):
                    start = max(0, i - footprint_radius)
                    end = min(len(t), i + footprint_radius + 1)
                    if stretch_value in t[start:end]:
                        sampled_pixels_inflated.append(stretch_value)
                    else:
                        sampled_pixels_inflated.append(t[i])
                        
            sampled_pixels_arr = np.array(sampled_pixels_inflated)
            sampled_bg = (np.sum(sampled_pixels_arr == target_value) / len(sampled_pixels_arr)) * 100 if len(sampled_pixels_arr) > 0 else 0.0
            
        res[f'BG_{label}_Error'] = sampled_bg - true_bg_pct
        
    fetch_intervals = {'25cm': 5, '50cm': 10, '100cm': 20, '200cm': 40}
    for label, step in fetch_intervals.items():
        fetch_sampled = np.concatenate([f[::step] for f in all_fetch_vals]) 
        sampled_mean = np.mean(fetch_sampled[fetch_sampled >= 0]) if fetch_sampled.size > 0 else 0.0
        res[f'Fetch_{label}_Error'] = sampled_mean - exact_fetch
    
    all_nri_gaps = np.concatenate([get_gap_lengths(t, target_value, cell_size) for t in all_nri_transects])
    res['Gap_0_24_0cm_Error'] = ((np.sum(all_nri_gaps[all_nri_gaps < 0.25]) / total_nri_length_m) * 100) - ex_0_24
    res['Gap_25_50_0cm_Error'] = ((np.sum(all_nri_gaps[(all_nri_gaps >= 0.25) & (all_nri_gaps <= 0.50)]) / total_nri_length_m) * 100) - ex_25_50
    res['Gap_51_100_0cm_Error'] = ((np.sum(all_nri_gaps[(all_nri_gaps >= 0.51) & (all_nri_gaps <= 1.00)]) / total_nri_length_m) * 100) - ex_51_100
    res['Gap_101_200_0cm_Error'] = ((np.sum(all_nri_gaps[(all_nri_gaps >= 1.01) & (all_nri_gaps <= 2.00)]) / total_nri_length_m) * 100) - ex_101_200
    res['Gap_gt_200_0cm_Error'] = ((np.sum(all_nri_gaps[all_nri_gaps > 2.00]) / total_nri_length_m) * 100) - ex_gt_200
    
    return res

# ====================================================================
# 3. MAIN EXECUTION BLOCK 
# ====================================================================
if __name__ == '__main__':
    num_iterations = 19000 
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
    cores = multiprocessing.cpu_count() - 3
    
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
            # Prints progress to the console every 100 iterations
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
            'Exact_Fetch_Mean': np.mean(x['Exact_Fetch'])
        }
        
        for scale in ['0cm', '25cm', '50cm', '100cm', '200cm']:
            d[f'BG_{scale}_MAE'] = np.mean(np.abs(x[f'BG_{scale}_Error']))
            d[f'BG_{scale}_MRE'] = calc_mre(x[f'BG_{scale}_Error'], x['True_BG_Pct'])
            d[f'BG_{scale}_Bias'] = np.mean(x[f'BG_{scale}_Error'])
            
        for scale in ['25cm', '50cm', '100cm', '200cm']:
            d[f'Fetch_{scale}_MAE'] = np.mean(np.abs(x[f'Fetch_{scale}_Error']))
            d[f'Fetch_{scale}_MRE'] = calc_mre(x[f'Fetch_{scale}_Error'], x['Exact_Fetch'])
            d[f'Fetch_{scale}_Bias'] = np.mean(x[f'Fetch_{scale}_Error'])
            
        for gap in ['Gap_0_24', 'Gap_25_50', 'Gap_51_100', 'Gap_101_200', 'Gap_gt_200']:
            d[f'{gap}_0cm_MAE'] = np.mean(np.abs(x[f'{gap}_0cm_Error']))
            d[f'{gap}_0cm_MRE'] = calc_mre(x[f'{gap}_0cm_Error'], x[f'Exact_{gap}'])
            d[f'{gap}_0cm_Bias'] = np.mean(x[f'{gap}_0cm_Error'])
            
        return pd.Series(d)

    # Calculate aggregations and SORT the dataframe to ensure lines don't zig-zag
    mae_df = df.groupby('BG_Bin').apply(aggregate_metrics).reset_index()
    mae_df = mae_df.sort_values('True_BG_Mean').reset_index(drop=True)

    print("\n" + "="*50)
    print("Average Exact Mean Fetch per Bare Ground Bin")
    print("="*50)
    for _, row in mae_df.iterrows():
        print(f"Target BG: {int(row['BG_Bin'])}%  ->  Mean Fetch: {row['Exact_Fetch_Mean']:.4f} m")
    print("="*50 + "\n")

    plot_config = [
        ('BG', 'Total Bare Ground (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
        ('Fetch', 'Mean Fetch (m)', ['25cm', '50cm', '100cm', '200cm'], None),
        ('Gap_0_24', 'Canopy Gap 0-24cm (%)', ['0cm'], 'steelblue'),
        ('Gap_25_50', 'Canopy Gap 25-50cm (%)', ['0cm'], 'cadetblue'),
        ('Gap_51_100', 'Canopy Gap 51-100cm (%)', ['0cm'], 'mediumseagreen'),
        ('Gap_101_200', 'Canopy Gap 101-200cm (%)', ['0cm'], 'darkorange'),
        ('Gap_gt_200', 'Canopy Gap >200cm (%)', ['0cm'], 'firebrick')
    ]

    scale_colors = {
        '0cm': 'black',
        '25cm': 'forestgreen',
        '50cm': 'dodgerblue',
        '100cm': 'darkorange',
        '200cm': 'crimson'
    }

    fig, axes = plt.subplots(3, 7, figsize=(28, 12), constrained_layout=True)
    fig.suptitle("Hybrid Simulation Metrics (MAE, MRE, Bias) Across Bare Ground Gradient", fontsize=20, weight='bold')
    
    for col, (prefix, title, scales, col_color) in enumerate(plot_config):
        ax_mae = axes[0, col]
        ax_mre = axes[1, col]
        ax_bias = axes[2, col]
        
        for scale in scales:
            var_base = f'{prefix}_{scale}'
            color = scale_colors[scale] if len(scales) > 1 else col_color
            label = f'{scale} Point' if scale != '0cm' else '0cm Continuous'
            line_kws = {'marker': 'o', 'color': color, 'linestyle': '-', 'linewidth': 2.5, 'markersize': 7}
            
            ax_mae.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MAE'], label=label, **line_kws)
            ax_mre.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MRE'], **line_kws)
            ax_bias.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Bias'], **line_kws)
            
        for ax in [ax_mae, ax_mre, ax_bias]:
            ax.set_xlim(0, 100)
            ax.set_xticks(np.arange(0, 81, 20)) 
            
        ax_mae.set_title(title, fontsize=15, pad=12)
        ax_mae.grid(True, alpha=0.3)
        ax_mae.legend(loc='upper left', fontsize=10) 
        
        if col == 0: ax_mae.set_ylabel("MAE", fontsize=13)
            
        ax_mre.grid(True, alpha=0.3)
        if col == 0: ax_mre.set_ylabel("MRE (%)", fontsize=13)
            
        ax_bias.axhline(0, color='gray', linestyle='--', linewidth=1.5)
        ax_bias.grid(True, alpha=0.3)
        ax_bias.set_xlabel("True Bare Ground (%)", fontsize=12)
        if col == 0: ax_bias.set_ylabel("Mean Bias", fontsize=13)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    img_path = os.path.join(script_dir, 'Hybrid_Simulation_Metrics_MultiScale.png')
    csv_path = os.path.join(script_dir, 'Hybrid_Simulation_Metrics_MultiScale.csv')
    
    plt.savefig(img_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    mae_df.to_csv(csv_path, index=False)
    print(f"\nResults saved to:\n  -> {img_path}\n  -> {csv_path}")
