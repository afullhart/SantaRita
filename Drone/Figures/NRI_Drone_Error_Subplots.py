import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ====================================================================
# 1. FOLDER SETUP & DATA LOADING
# ====================================================================
base_dir = r'C:\Users\andre'
csv_path = r'C:\Users\andre\ScatterPlots\SRER_NRI_Plots_110m.csv'

os.makedirs(base_dir, exist_ok=True)

# Load the drone imagery results
df = pd.read_csv(csv_path)

# ====================================================================
# 2. BINNING & AGGREGATION
# ====================================================================
# Create 5% bins for the Exact Bare Ground Percentage (0-100)
bins = np.arange(0, 105, 5)
df['BG_Bin'] = pd.cut(df['Exact_BGR_Pct'], bins=bins, include_lowest=True, right=False)

cover_scales = ['0cm', '25cm', '50cm', '100cm', '200cm']
fetch_scales = ['0cm', '25cm', '50cm', '100cm', '200cm']
gap_cols = ['Gap_0_24', 'Gap_25_50', 'Gap_51_100', 'Gap_101_200', 'Gap_gt_200']

results = []

# Group by the 5% bins and aggregate the values (Errors removed for Value-Only Plot)
for bin_val, group in df.groupby('BG_Bin', observed=False):
    # FILTER: Skip bins with fewer than 3 plots to prevent volatile edge-case spikes
    if len(group) < 3:
        continue
    
    # Extract the sample size and true means of the current bin
    d = {
        'True_BG_Mean': group['Exact_BGR_Pct'].mean(),
        'True_Fetch_Mean': group['Exact_Fetch_m'].mean(),
        'Sample_Size': len(group) 
    }
    
    # Extract LPI, Herbaceous, and Woody means
    if 'Exact_LPI_Pct' in group.columns:
        d['True_LPI_Mean'] = group['Exact_LPI_Pct'].mean()
    if 'Exact_Herb_Pct' in group.columns:
        d['True_Herb_Mean'] = group['Exact_Herb_Pct'].mean()
    if 'Exact_Woody_Pct' in group.columns:
        d['True_Woody_Mean'] = group['Exact_Woody_Pct'].mean()
        
    for gap in gap_cols:
        if f'Exact_{gap}' in group.columns:
            d[f'True_{gap}_Mean'] = group[f'Exact_{gap}'].mean()
    
    # Aggregate Sampled Cover Values
    for scale in cover_scales:
        if f'NRI_BGR_{scale}_Pct' in group.columns:
            d[f'BG_{scale}_Val'] = group[f'NRI_BGR_{scale}_Pct'].mean()
        if f'NRI_Herb_{scale}_Pct' in group.columns:
            d[f'Herb_{scale}_Val'] = group[f'NRI_Herb_{scale}_Pct'].mean()
        if f'NRI_Woody_{scale}_Pct' in group.columns:
            d[f'Woody_{scale}_Val'] = group[f'NRI_Woody_{scale}_Pct'].mean()
            
    # Aggregate Fetch
    for scale in fetch_scales:
        if f'NRI_Fetch_{scale}' in group.columns:
            d[f'Fetch_{scale}_Val'] = group[f'NRI_Fetch_{scale}'].mean()
            
    # Aggregate Canopy Gaps
    for gap in gap_cols:
        if f'NRI_{gap}' in group.columns:
            d[f'{gap}_0cm_Val'] = group[f'NRI_{gap}'].mean()
        
    results.append(d)

# Create final dataframe and drop any bins that ended up empty or NaN
val_df = pd.DataFrame(results).dropna(subset=['True_BG_Mean']).sort_values('True_BG_Mean')

# ====================================================================
# 3. HELPER: ROBUST Y-AXIS AUTOSCALING
# ====================================================================
def autoscale_y_robust(ax, margin=0.05, force_zero=False):
    """Autoscales the y-axis tightly based on the actual lines plotted on the axes."""
    min_y = np.inf
    max_y = -np.inf
    has_data = False
    
    for line in ax.lines:
        y_data = line.get_ydata()
        valid_y = y_data[np.isfinite(y_data)]
        if len(valid_y) > 0:
            min_y = min(min_y, valid_y.min())
            max_y = max(max_y, valid_y.max())
            has_data = True
            
    if has_data:
        if force_zero:
            min_y = min(min_y, 0)
            max_y = max(max_y, 0)
            
        y_range = max_y - min_y
        if y_range == 0:
            y_range = 1.0 # fallback if a line is completely flat
            
        ax.set_ylim(min_y - (y_range * margin), max_y + (y_range * margin))

# ====================================================================
# 4. PLOTTING CONFIGURATION (2x5 Grid)
# ====================================================================
plot_config = [
    # Top Row: BGR, HP, WP, MF, LPI
    ('BG', 'Total Bare Ground (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
    ('Herb', 'Herbaceous Cover (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
    ('Woody', 'Woody Cover (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
    ('Fetch', 'Mean Fetch (m)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
    ('LPI', 'Largest Patch Index (%)', [], 'black'), # Empty scales list skips the sampled line plot

    # Bottom Row: CGF fractions
    ('Gap_0_24', 'Canopy Gap 0-24 cm (%)', ['0cm'], 'steelblue'),
    ('Gap_25_50', 'Canopy Gap 25-50 cm (%)', ['0cm'], 'cadetblue'),
    ('Gap_51_100', 'Canopy Gap 51-100 cm (%)', ['0cm'], 'mediumseagreen'),
    ('Gap_101_200', 'Canopy Gap 101-200 cm (%)', ['0cm'], 'darkorange'),
    ('Gap_gt_200', 'Canopy Gap +200 cm (%)', ['0cm'], 'firebrick')
]

scale_colors = {
    '0cm': 'black',
    '25cm': 'forestgreen',
    '50cm': 'dodgerblue',
    '100cm': 'darkorange',
    '200cm': 'crimson'
}

# ====================================================================
# 5. FIGURE GENERATION
# ====================================================================
fig, axes = plt.subplots(2, 5, figsize=(28, 10), constrained_layout=True)
fig.suptitle("Drone Imagery Metrics vs True Bare Ground (%)", fontsize=20, weight='bold')

for idx, (prefix, title, scales, col_color) in enumerate(plot_config):
    row = idx // 5
    col = idx % 5
    ax = axes[row, col]
    
    # Route the correct exact true value column
    if prefix == 'BG': true_col = 'True_BG_Mean'
    elif prefix == 'Fetch': true_col = 'True_Fetch_Mean'
    elif prefix == 'LPI': true_col = 'True_LPI_Mean'
    elif prefix == 'Herb': true_col = 'True_Herb_Mean'
    elif prefix == 'Woody': true_col = 'True_Woody_Mean'
    elif prefix.startswith('Gap'): true_col = f'True_{prefix}_Mean'

    # Plot Exact Reference Line
    if true_col in val_df.columns:
        # LPI stands alone, so setting its exact line to a bold solid color makes it stand out
        linestyle = '-' if prefix == 'LPI' else '--'
        linecolor = 'black' if prefix == 'LPI' else 'gray'
        ax.plot(val_df['True_BG_Mean'], val_df[true_col], color=linecolor, linestyle=linestyle, linewidth=2, label='Exact True Value', zorder=1)
        
    # Plot Sampled Data Lines
    for scale in scales:
        var_base = f'{prefix}_{scale}'
        if f'{var_base}_Val' not in val_df.columns:
            continue
            
        color = scale_colors[scale] if len(scales) > 1 else col_color
        label = f'{scale} Point' if scale != '0cm' else '0cm Continuous'
        line_kws = {'marker': 'o', 'color': color, 'linestyle': '-', 'linewidth': 2.5, 'markersize': 7}
        
        ax.plot(val_df['True_BG_Mean'], val_df[f'{var_base}_Val'], label=label, **line_kws)
        
    # Apply strict Y-Axis autoscaling
    autoscale_y_robust(ax, force_zero=False)
        
    # Formatting
    ax.set_title(title, fontsize=15, pad=12)
    ax.axvline(50, color='gray', linestyle=':', linewidth=2)
    ax.set_xlim(0, 100)
    ax.set_xticks(np.arange(0, 81, 20))  
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=10)
    
    # Dynamic Axes Labeling based on the 2x5 Grid
    if col == 0: 
        if row == 0: ax.set_ylabel("Cover / Score", fontsize=13)
        else: ax.set_ylabel("Gap Fraction (%)", fontsize=13)
            
    if row == 1:
        ax.set_xlabel("True Bare Ground (%)", fontsize=12)

    # Re-apply Sample Size annotations on the bottom of the top-left plot (BGR)
    if idx == 0: 
        for _, r in val_df.iterrows():
            ax.text(r['True_BG_Mean'], 0.02, f"n={int(r['Sample_Size'])}", 
                        transform=ax.get_xaxis_transform(), 
                        fontsize=9, color='black', ha='center', va='bottom', rotation=90, 
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='lightgray'))

# Save Outputs
img_path = os.path.join(base_dir, 'Drone_Metrics_Values_2x5.png')
csv_path_out = os.path.join(base_dir, 'Drone_Metrics_Values_2x5.csv')

plt.savefig(img_path, dpi=300, bbox_inches='tight')
plt.show()

val_df.to_csv(csv_path_out, index=False)
print(f"Results saved to:\n  -> {img_path}\n  -> {csv_path_out}")
