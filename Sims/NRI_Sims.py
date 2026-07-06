import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import concurrent.futures
import multiprocessing
from shapely.geometry import Point, Polygon, LineString
from shapely.ops import unary_union
from scipy.ndimage import distance_transform_edt
from matplotlib.path import Path

# ====================================================================
# 1. HELPER FUNCTIONS
# ====================================================================
def polygon_to_path(polygon):
    """Converts a Shapely Polygon to a matplotlib Path for fast rasterization."""
    vertices = list(polygon.exterior.coords)
    codes = [Path.MOVETO] + [Path.LINETO] * (len(vertices) - 2) + [Path.CLOSEPOLY]
    for interior in polygon.interiors:
        interior_vertices = list(interior.coords)
        vertices.extend(interior_vertices)
        codes.extend([Path.MOVETO] + [Path.LINETO] * (len(interior_vertices) - 2) + [Path.CLOSEPOLY])
    return Path(vertices, codes)

def get_gap_lengths(transect_array, gap_val, p_size):
    """Extracts lengths of continuous sequences matching gap_val."""
    padded = np.concatenate(([False], (transect_array == gap_val), [False]))
    diffs = np.diff(padded.astype(int))
    starts, ends = np.where(diffs == 1)[0], np.where(diffs == -1)[0]
    return (ends - starts) * p_size

# ==========================================
# Geometry Generation Functions
# ==========================================
def create_irregular_shrub(center_x, center_y, base_radius, num_points=12):
    """Generates an irregular polygon to represent a shrub."""
    angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
    points = []
    for angle in angles:
        r = base_radius * np.random.uniform(0.6, 1.4)
        points.append((center_x + r * np.cos(angle), center_y + r * np.sin(angle)))
    return Polygon(points)

def generate_random_shrubs(num_shrubs=60, plot_radius=55): 
    """Generates a list of randomly placed shrub polygons within the plot."""
    shrubs = []
    max_shrub_radius = 5.0
    generation_radius = plot_radius + max_shrub_radius

    for _ in range(num_shrubs):
        radius = np.random.uniform(1.0, 5.0)
        r = generation_radius * np.sqrt(np.random.uniform(0, 1))
        theta = np.random.uniform(0, 2 * np.pi)

        x = r * np.cos(theta)
        y = r * np.sin(theta)
        shrubs.append(create_irregular_shrub(x, y, radius))
    return shrubs

# ====================================================================
# 2. ISOLATED WORKER FUNCTION
# ====================================================================
def process_simulation_iteration(task):
    target_bg, plot_radius, hub_radius, cell_size = task
    
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))

    # Grid Setup
    grid_size = int((plot_radius * 2) / cell_size)
    x = np.linspace(-plot_radius, plot_radius, grid_size)
    y = np.linspace(-plot_radius, plot_radius, grid_size)
    X, Y = np.meshgrid(x, y)
    
    dist_from_center = np.sqrt(X**2 + Y**2)
    valid_mask = dist_from_center <= plot_radius
    total_valid_pixels = np.sum(valid_mask)
    
    # Spoke Setup
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
    
    E_r_squared = (1**2 + 1*5 + 5**2) / 3
    mean_patch_area = np.pi * E_r_squared
    area_ratio = (np.pi * plot_radius**2) / mean_patch_area

    # Shape Generation
    is_inverted = target_bg <= 50
    target_coverage = target_bg / 100.0 if is_inverted else 1 - (target_bg / 100.0)
    
    num_shapes = int(-area_ratio * np.log(1 - target_coverage))
    shapes = generate_random_shrubs(num_shapes, plot_radius)
    shape_union = unary_union(shapes)
    
    # Rasterization
    main_array = np.full((grid_size, grid_size), target_value, dtype=int)
    
    if not shape_union.is_empty:
        geoms = [shape_union] if shape_union.geom_type == 'Polygon' else shape_union.geoms
        points_flat = np.column_stack((X.flatten(), Y.flatten()))
        for geom in geoms:
            if geom.is_empty: continue
            path = polygon_to_path(geom)
            inside = path.contains_points(points_flat).reshape(grid_size, grid_size)
            main_array[inside] = shrub_value
            
    if is_inverted:
        main_array = np.where(main_array == target_value, shrub_value, target_value)
        
    main_array[~valid_mask] = -9999 
    
    # Exact 2D Metrics
    bare_pixels = np.sum(main_array[valid_mask] == target_value)
    true_bg_pct = (bare_pixels / total_valid_pixels) * 100
    
    is_bare_full = (main_array == target_value)
    dist_array = distance_transform_edt(is_bare_full) * cell_size
    valid_full_fetch = dist_array[valid_mask]
    exact_fetch = np.mean(valid_full_fetch) if valid_full_fetch.size > 0 else 0.0
    
    exhaust_gap_arrays = []
    masked_for_gaps = np.where(valid_mask, main_array, -9999)
    for r in range(grid_size):
        exhaust_gap_arrays.append(get_gap_lengths(masked_for_gaps[r, :], target_value, cell_size))
    for c in range(grid_size):
        exhaust_gap_arrays.append(get_gap_lengths(masked_for_gaps[:, c], target_value, cell_size))
        
    all_exhaust_gaps = np.concatenate(exhaust_gap_arrays)
    ex_0_24 = (np.sum(all_exhaust_gaps[all_exhaust_gaps < 0.25]) / total_exhaust_length_m) * 100
    ex_25_50 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 0.25) & (all_exhaust_gaps <= 0.50)]) / total_exhaust_length_m) * 100
    ex_51_100 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 0.51) & (all_exhaust_gaps <= 1.00)]) / total_exhaust_length_m) * 100
    ex_101_200 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 1.01) & (all_exhaust_gaps <= 2.00)]) / total_exhaust_length_m) * 100
    ex_gt_200 = (np.sum(all_exhaust_gaps[all_exhaust_gaps > 2.00]) / total_exhaust_length_m) * 100
    
    # Synthetic NRI Metrics
    all_nri_transects = []
    all_fetch_vals = []
    
    for rows, cols in spoke_indices:
        all_nri_transects.append(main_array[rows, cols])
        all_fetch_vals.append(dist_array[rows, cols])
        
    total_nri_pixels = sum(len(t) for t in all_nri_transects)
    nri_bare_pixels = sum(np.sum(t == target_value) for t in all_nri_transects)
    sampled_bg_pct = (nri_bare_pixels / total_nri_pixels) * 100 if total_nri_pixels > 0 else 0.0
        
    fetch_sampled = np.concatenate([f[::5] for f in all_fetch_vals]) 
    nri_fetch_25cm = np.mean(fetch_sampled[fetch_sampled >= 0]) if fetch_sampled.size > 0 else 0.0
    
    all_nri_gaps = np.concatenate([get_gap_lengths(t, target_value, cell_size) for t in all_nri_transects])
    nri_0_24 = (np.sum(all_nri_gaps[all_nri_gaps < 0.25]) / total_nri_length_m) * 100
    nri_25_50 = (np.sum(all_nri_gaps[(all_nri_gaps >= 0.25) & (all_nri_gaps <= 0.50)]) / total_nri_length_m) * 100
    nri_51_100 = (np.sum(all_nri_gaps[(all_nri_gaps >= 0.51) & (all_nri_gaps <= 1.00)]) / total_nri_length_m) * 100
    nri_101_200 = (np.sum(all_nri_gaps[(all_nri_gaps >= 1.01) & (all_nri_gaps <= 2.00)]) / total_nri_length_m) * 100
    nri_gt_200 = (np.sum(all_nri_gaps[all_nri_gaps > 2.00]) / total_nri_length_m) * 100
    
    # Export Exact values in addition to errors for proper MRE calculation
    return {
        'BG_Bin': target_bg,
        'True_BG_Pct': true_bg_pct,
        'BG_Error': sampled_bg_pct - true_bg_pct,
        
        'Exact_Fetch': exact_fetch,
        'Fetch_Error': nri_fetch_25cm - exact_fetch,
        
        'Exact_Gap_0_24': ex_0_24,
        'Gap_0_24_Error': nri_0_24 - ex_0_24,
        
        'Exact_Gap_25_50': ex_25_50,
        'Gap_25_50_Error': nri_25_50 - ex_25_50,
        
        'Exact_Gap_51_100': ex_51_100,
        'Gap_51_100_Error': nri_51_100 - ex_51_100,
        
        'Exact_Gap_101_200': ex_101_200,
        'Gap_101_200_Error': nri_101_200 - ex_101_200,
        
        'Exact_Gap_gt_200': ex_gt_200,
        'Gap_gt_200_Error': nri_gt_200 - ex_gt_200
    }

# ====================================================================
# 3. MAIN EXECUTION BLOCK 
# ====================================================================
if __name__ == '__main__':
    num_iterations = 15000 # Total number of simulations per variable.
    plot_radius = 55
    hub_radius = 5
    cell_size = 0.05
    
    target_bgs = np.arange(15, 90, 5)
    iters_per_bin = max(1, num_iterations // len(target_bgs))
    
    tasks = []
    for target_bg in target_bgs:
        for _ in range(iters_per_bin):
            tasks.append((target_bg, plot_radius, hub_radius, cell_size))
            
    total_tasks = len(tasks)
    cores = multiprocessing.cpu_count() - 2
    
    print(f"{'='*50}")
    print(f"Starting Multiprocessing Simulation")
    print(f"Target Bins: {len(target_bgs)}")
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

    # ====================================================================
    # 4. DATA AGGREGATION & PLOTTING (Updated for 3x7 Grid)
    # ====================================================================
    df = pd.DataFrame(results)
    
    def calc_mre(error, exact):
        """Helper to calculate Mean Relative Error stably at the bin level."""
        m_exact = np.mean(exact)
        return 0.0 if m_exact == 0 else (np.mean(np.abs(error)) / m_exact) * 100

    mae_df = df.groupby('BG_Bin').apply(lambda x: pd.Series({
        'True_BG_Mean': np.mean(x['True_BG_Pct']),
        
        'Fetch_MAE': np.mean(np.abs(x['Fetch_Error'])),
        'Fetch_MRE': calc_mre(x['Fetch_Error'], x['Exact_Fetch']),
        'Fetch_Bias': np.mean(x['Fetch_Error']),
        
        'Gap_0_24_MAE': np.mean(np.abs(x['Gap_0_24_Error'])),
        'Gap_0_24_MRE': calc_mre(x['Gap_0_24_Error'], x['Exact_Gap_0_24']),
        'Gap_0_24_Bias': np.mean(x['Gap_0_24_Error']),
        
        'Gap_25_50_MAE': np.mean(np.abs(x['Gap_25_50_Error'])),
        'Gap_25_50_MRE': calc_mre(x['Gap_25_50_Error'], x['Exact_Gap_25_50']),
        'Gap_25_50_Bias': np.mean(x['Gap_25_50_Error']),
        
        'Gap_51_100_MAE': np.mean(np.abs(x['Gap_51_100_Error'])),
        'Gap_51_100_MRE': calc_mre(x['Gap_51_100_Error'], x['Exact_Gap_51_100']),
        'Gap_51_100_Bias': np.mean(x['Gap_51_100_Error']),
        
        'Gap_101_200_MAE': np.mean(np.abs(x['Gap_101_200_Error'])),
        'Gap_101_200_MRE': calc_mre(x['Gap_101_200_Error'], x['Exact_Gap_101_200']),
        'Gap_101_200_Bias': np.mean(x['Gap_101_200_Error']),
        
        'Gap_gt_200_MAE': np.mean(np.abs(x['Gap_gt_200_Error'])),
        'Gap_gt_200_MRE': calc_mre(x['Gap_gt_200_Error'], x['Exact_Gap_gt_200']),
        'Gap_gt_200_Bias': np.mean(x['Gap_gt_200_Error']),

        'BG_MAE': np.mean(np.abs(x['BG_Error'])),
        'BG_MRE': calc_mre(x['BG_Error'], x['True_BG_Pct']),
        'BG_Bias': np.mean(x['BG_Error'])
    })).reset_index()

    # Create a 3x7 grid with a wide layout to accommodate 21 subplots
    fig, axes = plt.subplots(3, 7, figsize=(28, 12), constrained_layout=True)
    fig.suptitle("Simulation Metrics (MAE, MRE, Bias) Across Bare Ground Gradient", fontsize=20, weight='bold')
    
    # NEW ORDER: Total Bare Ground moved to index 0
    vars_to_plot = [
        ('BG', 'Total Bare Ground', 'black'),
        ('Fetch', 'Mean Fetch (m)', 'indigo'),
        ('Gap_0_24', 'Canopy Gap 0-24cm', 'steelblue'),
        ('Gap_25_50', 'Canopy Gap 25-50cm', 'cadetblue'),
        ('Gap_51_100', 'Canopy Gap 51-100cm', 'mediumseagreen'),
        ('Gap_101_200', 'Canopy Gap 101-200cm', 'darkorange'),
        ('Gap_gt_200', 'Canopy Gap >200cm', 'firebrick')
    ]
    
    for col, (prefix, title, color) in enumerate(vars_to_plot):
        
        # --- ROW 0: MAE ---
        ax_mae = axes[0, col]
        ax_mae.plot(mae_df['True_BG_Mean'], mae_df[f'{prefix}_MAE'], marker='o', color=color, linestyle='-', linewidth=2.5, markersize=7)
        ax_mae.set_title(title, fontsize=15, pad=12)
        ax_mae.axvline(50, color='gray', linestyle=':', linewidth=2, label='Phase Shift (50%)' if col==0 else "")
        ax_mae.grid(True, alpha=0.3)
        if col == 0:
            ax_mae.set_ylabel("MAE", fontsize=13)
            ax_mae.legend(fontsize=10)
            
        # --- ROW 1: MRE (%) ---
        ax_mre = axes[1, col]
        ax_mre.plot(mae_df['True_BG_Mean'], mae_df[f'{prefix}_MRE'], marker='o', color=color, linestyle='-', linewidth=2.5, markersize=7)
        ax_mre.axvline(50, color='gray', linestyle=':', linewidth=2)
        ax_mre.grid(True, alpha=0.3)
        if col == 0:
            ax_mre.set_ylabel("MRE (%)", fontsize=13)
            
        # --- ROW 2: Bias ---
        ax_bias = axes[2, col]
        ax_bias.plot(mae_df['True_BG_Mean'], mae_df[f'{prefix}_Bias'], marker='o', color=color, linestyle='-', linewidth=2.5, markersize=7)
        ax_bias.axvline(50, color='gray', linestyle=':', linewidth=2)
        ax_bias.axhline(0, color='gray', linestyle='--', linewidth=1.5) # Zero reference line for Bias
        ax_bias.grid(True, alpha=0.3)
        ax_bias.set_xlabel("True Bare Ground (%)", fontsize=12)
        if col == 0:
            ax_bias.set_ylabel("Mean Bias", fontsize=13)

    # Save the unified figure using absolute paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    img_path = os.path.join(script_dir, 'Full_Simulation_Metrics.png')
    csv_path = os.path.join(script_dir, 'Full_Simulation_Metrics.csv')
    
    plt.savefig(img_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    # Save raw data for external review
    mae_df.to_csv(csv_path, index=False)
    print(f"\nResults saved to:\n  -> {img_path}\n  -> {csv_path}")
