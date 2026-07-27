import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import os

# ====================================================================
# FOLDERS & FILES
# ====================================================================
data_folder = r'C:\Users\andre\ScatterPlots'
out_folder = r'C:\Users\andre\ScatterPlots\Random_Sampling'

csv_files = {
    '10m_Grid': 'SRER_Grid_10m_Random.csv',
    '30m_Grid': 'SRER_Grid_30m_Random.csv',
    '110m_NRI': 'SRER_NRI_Plots_110m_Random.csv'
}

# --- HELPER FUNCTION FOR DUAL-AXIS SLOPE & R2 PLOTTING ---
def plot_slope_r2(ax, x_vals, slopes, r2_vals, title, x_label):
    color_slope = 'tab:blue'
    color_r2 = 'tab:orange'
    
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel(x_label, fontsize=12)
    
    # Left Axis: OLS Slope
    ax.set_ylabel('OLS Slope (Exact vs Sampled)', color=color_slope, fontsize=12)
    ax.plot(x_vals, slopes, marker='o', linestyle='-', color=color_slope, label='Regression Slope')
    ax.tick_params(axis='y', labelcolor=color_slope)
    
    # Add target line for perfect 1:1 Relationship
    ax.axhline(y=1.0, color='red', linestyle='--', linewidth=2, alpha=0.7, label='Perfect 1.0 Slope')
    
    # Format X-Axis (Logarithmic)
    ax.set_xscale('log')
    ax.set_xticks(x_vals)
    ax.set_xticklabels([str(x) for x in x_vals], rotation=45)
    ax.grid(True, which="major", ls="--", alpha=0.5)
    
    # Right Axis: R-squared
    ax2 = ax.twinx()
    ax2.set_ylabel('R-Squared ($R^2$)', color=color_r2, fontsize=12)
    ax2.plot(x_vals, r2_vals, marker='s', linestyle='-', color=color_r2, alpha=0.8, label='$R^2$ Value')
    ax2.tick_params(axis='y', labelcolor=color_r2)
    
    # Combined legend (placed carefully to avoid crowding)
    lines_1, labels_1 = ax.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax2.legend(lines_1 + lines_2, labels_1 + labels_2, loc='lower right')

# --- HELPER TO EXTRACT SLOPE AND R2 ---
def get_ols_metrics(df, x_col, y_col):
    subset = df[[x_col, y_col]].dropna()
    if len(subset) < 2:
        return np.nan, np.nan
    
    X = subset[x_col].values
    Y = subset[y_col].values
    
    X_sm = sm.add_constant(X)
    model = sm.OLS(Y, X_sm).fit()
    
    # slope is at index 1 (index 0 is the intercept)
    return model.params[1], model.rsquared

# ====================================================================
# MAIN PROCESSING LOOP
# ====================================================================
for scale_name, csv_filename in csv_files.items():
    csv_path = os.path.join(data_folder, csv_filename)
    
    if not os.path.exists(csv_path):
        print(f"Skipping {scale_name}: Could not find {csv_path}")
        continue
        
    print(f"\n--- Calculating Slopes for {scale_name} ---")
    df = pd.read_csv(csv_path)

    # 1. Dynamically find sampling increments
    pt_cols = [c for c in df.columns if c.startswith('BGR_pt_')]
    if not pt_cols:
        continue
    points = sorted([int(c.split('_')[-1]) for c in pt_cols])

    ln_cols = [c for c in df.columns if c.startswith('Gap_0_24_L_')]
    lines = sorted([int(c.split('_')[-1]) for c in ln_cols])

    # 2. Storage for calculated metrics
    bgr_slopes, bgr_r2 = [], []
    herb_slopes, herb_r2 = [], []
    woody_slopes, woody_r2 = [], []
    fetch_slopes, fetch_r2 = [], []
    hw_slopes, hw_r2 = [], []

    print(" -> Processing Point Data...")
    for pt in points:
        # BGR
        s, r = get_ols_metrics(df, f'BGR_pt_{pt}', 'BGR_Exact')
        bgr_slopes.append(s); bgr_r2.append(r)
        
        # Herb Pct
        s, r = get_ols_metrics(df, f'Herb_Pct_pt_{pt}', 'Herb_Pct_Exact')
        herb_slopes.append(s); herb_r2.append(r)

        # Woody Pct
        s, r = get_ols_metrics(df, f'Woody_Pct_pt_{pt}', 'Woody_Pct_Exact')
        woody_slopes.append(s); woody_r2.append(r)

        # Fetch
        s, r = get_ols_metrics(df, f'Fetch_pt_{pt}', 'Fetch_Exact')
        fetch_slopes.append(s); fetch_r2.append(r)

        # Herb-to-Woody Ratio
        s, r = get_ols_metrics(df, f'HW_Ratio_pt_{pt}', 'HW_Ratio_Exact')
        hw_slopes.append(s); hw_r2.append(r)

    # 3. Storage for Canopy Gaps
    gap_cats = {
        'Gap_0_24': 'Canopy Gap: 0-24 cm',
        'Gap_25_50': 'Canopy Gap: 25-50 cm',
        'Gap_51_100': 'Canopy Gap: 51-100 cm',
        'Gap_101_200': 'Canopy Gap: 101-200 cm',
        'Gap_gt_200': 'Canopy Gap: >200 cm'
    }
    
    print(" -> Processing Line Data...")
    gap_data = {}
    for cat, title in gap_cats.items():
        g_slopes, g_r2 = [], []
        exact_col = f'{cat}_Exact'
        
        for ln in lines:
            s, r = get_ols_metrics(df, f'{cat}_L_{ln}', exact_col)
            g_slopes.append(s)
            g_r2.append(r)
            
        gap_data[cat] = {'title': title, 'slopes': g_slopes, 'r2': g_r2}

    # 4. Generate the 10-Subplot Figure
    fig, axes = plt.subplots(nrows=10, ncols=1, figsize=(12, 50))

    plot_slope_r2(axes[0], points, bgr_slopes, bgr_r2, 
                  'Bare Ground Percentage', 'Number of Randomly Sampled Points (Log Scale)')
    plot_slope_r2(axes[1], points, herb_slopes, herb_r2, 
                  'Herb Cover Percentage', 'Number of Randomly Sampled Points (Log Scale)')
    plot_slope_r2(axes[2], points, woody_slopes, woody_r2, 
                  'Woody Cover Percentage', 'Number of Randomly Sampled Points (Log Scale)')
    plot_slope_r2(axes[3], points, fetch_slopes, fetch_r2, 
                  'Mean Fetch', 'Number of Randomly Sampled Points (Log Scale)')
    plot_slope_r2(axes[4], points, hw_slopes, hw_r2, 
                  'Herb-to-Woody Ratio', 'Number of Randomly Sampled Points (Log Scale)')

    ax_idx = 5
    for cat in gap_cats.keys():
        plot_slope_r2(axes[ax_idx], lines, gap_data[cat]['slopes'], gap_data[cat]['r2'], 
                      gap_data[cat]['title'], 'Virtual Transect Length in Meters (Log Scale)')
        ax_idx += 1

    # 5. Final Formatting and Export
    fig.suptitle(f'Convergence of OLS Slope to 1:1 ({scale_name.replace("_", " ")})', 
                 fontsize=20, fontweight='bold', y=0.99)
                 
    fig.tight_layout(pad=4.0)  

    output_img_path = os.path.join(out_folder, f'{scale_name}_Slope_Convergence.png')
    plt.savefig(output_img_path, dpi=300)
    plt.close() 
    
    print(f"Saved slope convergence figure to: {output_img_path}")

print("\nAll slope plotting complete!")
