import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ====================================================================
# 1. FOLDER SETUP & DATA LOADING
# ====================================================================
base_dir = r'C:\Users\andre'
csv_path = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\SRER_NRI_Plots_110m.csv'

os.makedirs(base_dir, exist_ok=True)

# Load the drone imagery results
df = pd.read_csv(csv_path)

# ====================================================================
# 2. BINNING & ERROR CALCULATION
# ====================================================================
# Create 5% bins for the Exact Bare Ground Percentage (0-100)
bins = np.arange(0, 105, 5)
df['BG_Bin'] = pd.cut(df['Exact_BGR_Pct'], bins=bins, include_lowest=True, right=False)

# Calculate Errors (Sampled - Exact)
bg_scales = ['0cm', '25cm', '50cm', '100cm', '200cm']
for scale in bg_scales:
    df[f'BG_{scale}_Error'] = df[f'NRI_BGR_{scale}_Pct'] - df['Exact_BGR_Pct']

fetch_scales = ['25cm', '50cm', '100cm', '200cm']
for scale in fetch_scales:
    df[f'Fetch_{scale}_Error'] = df[f'NRI_Fetch_{scale}'] - df['Exact_Fetch_m']

gap_cols = ['Gap_0_24', 'Gap_25_50', 'Gap_51_100', 'Gap_101_200', 'Gap_gt_200']
for gap in gap_cols:
    df[f'{gap}_0cm_Error'] = df[f'NRI_{gap}'] - df[f'Exact_{gap}']

# ====================================================================
# 3. DATA AGGREGATION
# ====================================================================
def calc_mre(error_series, exact_series):
    m_exact = exact_series.mean()
    if m_exact == 0 or pd.isna(m_exact):
        return np.nan
    return (error_series.abs().mean() / m_exact) * 100

results = []

# Group by the 5% bins and calculate metrics
for bin_val, group in df.groupby('BG_Bin', observed=False):
    if group.empty:
        continue
    
    d = {'True_BG_Mean': group['Exact_BGR_Pct'].mean()}
    
    # Aggregate BG
    for scale in bg_scales:
        d[f'BG_{scale}_MAE'] = group[f'BG_{scale}_Error'].abs().mean()
        d[f'BG_{scale}_MRE'] = calc_mre(group[f'BG_{scale}_Error'], group['Exact_BGR_Pct'])
        d[f'BG_{scale}_Bias'] = group[f'BG_{scale}_Error'].mean()
        
    # Aggregate Fetch
    for scale in fetch_scales:
        d[f'Fetch_{scale}_MAE'] = group[f'Fetch_{scale}_Error'].abs().mean()
        d[f'Fetch_{scale}_MRE'] = calc_mre(group[f'Fetch_{scale}_Error'], group['Exact_Fetch_m'])
        d[f'Fetch_{scale}_Bias'] = group[f'Fetch_{scale}_Error'].mean()
        
    # Aggregate Gaps
    for gap in gap_cols:
        d[f'{gap}_0cm_MAE'] = group[f'{gap}_0cm_Error'].abs().mean()
        d[f'{gap}_0cm_MRE'] = calc_mre(group[f'{gap}_0cm_Error'], group[f'Exact_{gap}'])
        d[f'{gap}_0cm_Bias'] = group[f'{gap}_0cm_Error'].mean()
        
    results.append(d)

# Create final dataframe and drop any bins that ended up empty or NaN
mae_df = pd.DataFrame(results).dropna(subset=['True_BG_Mean']).sort_values('True_BG_Mean')

# ====================================================================
# 4. PLOTTING
# ====================================================================
# Dynamic Plotting Configuration identical to simulation script
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
fig.suptitle("Drone Imagery Metrics (MAE, MRE, Bias) Across Bare Ground Gradient", fontsize=20, weight='bold')

for col, (prefix, title, scales, col_color) in enumerate(plot_config):
    ax_mae = axes[0, col]
    ax_mre = axes[1, col]
    ax_bias = axes[2, col]
    
    for scale in scales:
        var_base = f'{prefix}_{scale}'
        
        # Use multi-scale dictionary if there are multiple lines. Otherwise, use the unique column color.
        color = scale_colors[scale] if len(scales) > 1 else col_color
        
        # Explicitly label the '0cm' scales appropriately for the legend
        label = f'{scale} Point' if scale != '0cm' else '0cm Continuous'
        
        # Line styles
        line_kws = {'marker': 'o', 'color': color, 'linestyle': '-', 'linewidth': 2.5, 'markersize': 7}
        
        ax_mae.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MAE'], label=label, **line_kws)
        ax_mre.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MRE'], **line_kws)
        ax_bias.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Bias'], **line_kws)
        
    # Format MAE Row
    ax_mae.set_title(title, fontsize=15, pad=12)
    ax_mae.axvline(50, color='gray', linestyle=':', linewidth=2)
    ax_mae.set_xlim(0, 95)  # Enforce x-axis limits
    ax_mae.grid(True, alpha=0.3)
    ax_mae.legend(loc='upper left', fontsize=10) 
    if col == 0: ax_mae.set_ylabel("MAE", fontsize=13)
        
    # Format MRE Row
    ax_mre.axvline(50, color='gray', linestyle=':', linewidth=2)
    ax_mre.set_xlim(0, 95)  # Enforce x-axis limits
    ax_mre.grid(True, alpha=0.3)
    if col == 0: ax_mre.set_ylabel("MRE (%)", fontsize=13)
        
    # Format Bias Row
    ax_bias.axvline(50, color='gray', linestyle=':', linewidth=2)
    ax_bias.axhline(0, color='gray', linestyle='--', linewidth=1.5)
    ax_bias.set_xlim(0, 95)  # Enforce x-axis limits
    ax_bias.grid(True, alpha=0.3)
    ax_bias.set_xlabel("True Bare Ground (%)", fontsize=12)
    if col == 0: ax_bias.set_ylabel("Mean Bias", fontsize=13)

# Save Outputs
img_path = os.path.join(base_dir, 'Drone_Metrics_MultiScale.png')
csv_path_out = os.path.join(base_dir, 'Drone_Metrics_MultiScale.csv')

plt.savefig(img_path, dpi=300, bbox_inches='tight')
plt.show()

mae_df.to_csv(csv_path_out, index=False)
print(f"Results saved to:\n  -> {img_path}\n  -> {csv_path_out}")
