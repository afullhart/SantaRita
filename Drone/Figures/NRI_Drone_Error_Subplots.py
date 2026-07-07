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
    # FILTER: Skip bins with fewer than 3 plots to prevent volatile edge-case spikes
    if len(group) < 3:
        continue
    
    # Extract the sample size of the current bin
    d = {
        'True_BG_Mean': group['Exact_BGR_Pct'].mean(),
        'Sample_Size': len(group) 
    }
    
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
# 4. HELPER: ROBUST Y-AXIS AUTOSCALING
# ====================================================================
def autoscale_y_robust(ax, margin=0.05, force_zero=False):
    """Autoscales the y-axis tightly based on the actual lines plotted on the axes."""
    min_y = np.inf
    max_y = -np.inf
    has_data = False
    
    # Extract data directly from the plotted matplotlib lines to avoid indexing issues
    for line in ax.lines:
        y_data = line.get_ydata()
        valid_y = y_data[np.isfinite(y_data)]
        if len(valid_y) > 0:
            min_y = min(min_y, valid_y.min())
            max_y = max(max_y, valid_y.max())
            has_data = True
            
    if has_data:
        if force_zero:
            # Ensure 0 is included within the data limits bounds
            min_y = min(min_y, 0)
            max_y = max(max_y, 0)
            
        y_range = max_y - min_y
        if y_range == 0:
            y_range = 1.0 # fallback if a line is completely flat
            
        ax.set_ylim(min_y - (y_range * margin), max_y + (y_range * margin))

# ====================================================================
# 5. PLOTTING
# ====================================================================
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
        color = scale_colors[scale] if len(scales) > 1 else col_color
        label = f'{scale} Point' if scale != '0cm' else '0cm Continuous'
        line_kws = {'marker': 'o', 'color': color, 'linestyle': '-', 'linewidth': 2.5, 'markersize': 7}
        
        # Plot the actual data lines
        ax_mae.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MAE'], label=label, **line_kws)
        ax_mre.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MRE'], **line_kws)
        ax_bias.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Bias'], **line_kws)
        
    # Apply strict Y-Axis autoscaling strictly to the plotted lines before adding structural lines
    autoscale_y_robust(ax_mae, force_zero=False)
    autoscale_y_robust(ax_mre, force_zero=False)
    autoscale_y_robust(ax_bias, force_zero=True) # Bias MUST include the 0 reference line!
        
    # Format MAE Row
    ax_mae.set_title(title, fontsize=15, pad=12)
    ax_mae.axvline(50, color='gray', linestyle=':', linewidth=2)
    ax_mae.set_xlim(0, 95)  
    ax_mae.grid(True, alpha=0.3)
    ax_mae.legend(loc='upper left', fontsize=10) 
    
    # Annotate sample sizes exclusively in the first MAE subplot
    if col == 0: 
        ax_mae.set_ylabel("MAE", fontsize=13)
        for _, row in mae_df.iterrows():
            ax_mae.text(row['True_BG_Mean'], 0.02, f"n={int(row['Sample_Size'])}", 
                        transform=ax_mae.get_xaxis_transform(), 
                        fontsize=9, color='black', ha='center', va='bottom', rotation=90, # Added rotation here
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='lightgray'))
        
    # Format MRE Row
    ax_mre.axvline(50, color='gray', linestyle=':', linewidth=2)
    ax_mre.set_xlim(0, 95)  
    ax_mre.grid(True, alpha=0.3)
    if col == 0: ax_mre.set_ylabel("MRE (%)", fontsize=13)
        
    # Format Bias Row
    ax_bias.axvline(50, color='gray', linestyle=':', linewidth=2)
    ax_bias.axhline(0, color='gray', linestyle='--', linewidth=1.5)
    ax_bias.set_xlim(0, 95)  
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
