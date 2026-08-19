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
# 2. BINNING & ERROR CALCULATION
# ====================================================================
# Create 5% bins for the Exact Bare Ground Percentage (0-100)
bins = np.arange(0, 105, 5)
df['BG_Bin'] = pd.cut(df['Exact_BGR_Pct'], bins=bins, include_lowest=True, right=False)

cover_scales = ['0cm', '25cm', '50cm', '100cm', '200cm']

# Calculate Errors (Sampled - Exact) for Bare Ground
for scale in cover_scales:
    df[f'BG_{scale}_Error'] = df[f'NRI_BGR_{scale}_Pct'] - df['Exact_BGR_Pct']

# Calculate Errors for Herbaceous and Woody
for scale in cover_scales:
    if f'NRI_Herb_{scale}_Pct' in df.columns and 'Exact_Herb_Pct' in df.columns:
        df[f'Herb_{scale}_Error'] = df[f'NRI_Herb_{scale}_Pct'] - df['Exact_Herb_Pct']
    if f'NRI_Woody_{scale}_Pct' in df.columns and 'Exact_Woody_Pct' in df.columns:
        df[f'Woody_{scale}_Error'] = df[f'NRI_Woody_{scale}_Pct'] - df['Exact_Woody_Pct']

# Calculate Errors for Fetch
fetch_scales = ['0cm', '25cm', '50cm', '100cm', '200cm']
for scale in fetch_scales:
    df[f'Fetch_{scale}_Error'] = df[f'NRI_Fetch_{scale}'] - df['Exact_Fetch_m']

# Calculate Errors for Canopy Gaps
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
    
    # Extract the sample size, true means, and exact geometric center of the bin
    d = {
        'True_BG_Mean': group['Exact_BGR_Pct'].mean(),
        'True_Fetch_Mean': group['Exact_Fetch_m'].mean(),
        'Sample_Size': len(group),
        'Bin_Center': bin_val.mid
    }
    
    if 'Exact_LPI_Pct' in group.columns:
        d['True_LPI_Mean'] = group['Exact_LPI_Pct'].mean()
    if 'Exact_Herb_Pct' in group.columns:
        d['True_Herb_Mean'] = group['Exact_Herb_Pct'].mean()
    if 'Exact_Woody_Pct' in group.columns:
        d['True_Woody_Mean'] = group['Exact_Woody_Pct'].mean()
        
    for gap in gap_cols:
        if f'Exact_{gap}' in group.columns:
            d[f'True_{gap}_Mean'] = group[f'Exact_{gap}'].mean()
    
    # Aggregate BG
    for scale in cover_scales:
        d[f'BG_{scale}_Val'] = group[f'NRI_BGR_{scale}_Pct'].mean()
        d[f'BG_{scale}_MAE'] = group[f'BG_{scale}_Error'].abs().mean()
        d[f'BG_{scale}_MRE'] = calc_mre(group[f'BG_{scale}_Error'], group['Exact_BGR_Pct'])
        d[f'BG_{scale}_Bias'] = group[f'BG_{scale}_Error'].mean()
        
    # Aggregate Herbaceous
    for scale in cover_scales:
        if f'Herb_{scale}_Error' in group.columns:
            d[f'Herb_{scale}_Val'] = group[f'NRI_Herb_{scale}_Pct'].mean()
            d[f'Herb_{scale}_MAE'] = group[f'Herb_{scale}_Error'].abs().mean()
            d[f'Herb_{scale}_MRE'] = calc_mre(group[f'Herb_{scale}_Error'], group['Exact_Herb_Pct'])
            d[f'Herb_{scale}_Bias'] = group[f'Herb_{scale}_Error'].mean()

    # Aggregate Woody
    for scale in cover_scales:
        if f'Woody_{scale}_Error' in group.columns:
            d[f'Woody_{scale}_Val'] = group[f'NRI_Woody_{scale}_Pct'].mean()
            d[f'Woody_{scale}_MAE'] = group[f'Woody_{scale}_Error'].abs().mean()
            d[f'Woody_{scale}_MRE'] = calc_mre(group[f'Woody_{scale}_Error'], group['Exact_Woody_Pct'])
            d[f'Woody_{scale}_Bias'] = group[f'Woody_{scale}_Error'].mean()
        
    # Aggregate Fetch
    for scale in fetch_scales:
        if f'NRI_Fetch_{scale}' in group.columns:
            d[f'Fetch_{scale}_Val'] = group[f'NRI_Fetch_{scale}'].mean()
            d[f'Fetch_{scale}_MAE'] = group[f'Fetch_{scale}_Error'].abs().mean()
            d[f'Fetch_{scale}_MRE'] = calc_mre(group[f'Fetch_{scale}_Error'], group['Exact_Fetch_m'])
            d[f'Fetch_{scale}_Bias'] = group[f'Fetch_{scale}_Error'].mean()
        
    # Aggregate Gaps
    for gap in gap_cols:
        if f'{gap}_0cm_Error' in group.columns:
            d[f'{gap}_0cm_Val'] = group[f'NRI_{gap}'].mean()
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
# 5. GLOBAL PLOTTING CONFIGURATION
# ====================================================================
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
plot_config_1 = [
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

fig1, axes1 = plt.subplots(4, 9, figsize=(36, 16), constrained_layout=True)
fig1.suptitle("Drone Imagery Metrics (Values, MAE, MRE, Bias) Across Bare Ground Gradient", fontsize=20, weight='bold')

for col, (prefix, title, scales, col_color) in enumerate(plot_config_1):
    ax_val = axes1[0, col]
    ax_mae = axes1[1, col]
    ax_mre = axes1[2, col]
    ax_bias = axes1[3, col]
    
    # Plot Exact Reference Lines on the Value Row
    if prefix == 'BG':
        ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_BG_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
    elif prefix == 'Fetch':
        ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_Fetch_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
    elif prefix.startswith('Gap'):
        ax_val.plot(mae_df['True_BG_Mean'], mae_df[f'True_{prefix}_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
    elif prefix == 'Herb' and 'True_Herb_Mean' in mae_df.columns:
        ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_Herb_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
    elif prefix == 'Woody' and 'True_Woody_Mean' in mae_df.columns:
        ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_Woody_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
        
    for scale in scales:
        var_base = f'{prefix}_{scale}'
        if f'{var_base}_Val' not in mae_df.columns:
            continue
            
        color = scale_colors[scale] if len(scales) > 1 else col_color
        label_name = f'{scale} Point' if scale != '0cm' else '0cm Continuous'
        line_kws = {'marker': 'o', 'color': color, 'linestyle': '-', 'linewidth': 2.5, 'markersize': 7}
        
        ax_val.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Val'], label=label_name, **line_kws)
        ax_mae.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MAE'], **line_kws)
        ax_mre.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MRE'], **line_kws)
        ax_bias.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Bias'], **line_kws)
        
    autoscale_y_robust(ax_val, force_zero=False)
    autoscale_y_robust(ax_mae, force_zero=False)
    autoscale_y_robust(ax_mre, force_zero=False)
    autoscale_y_robust(ax_bias, force_zero=True) 
        
    ax_val.set_title(title, fontsize=15, pad=12)
    ax_val.axvline(50, color='gray', linestyle=':', linewidth=2)
    ax_val.set_xlim(0, 100)
    ax_val.set_xticks(np.arange(0, 81, 20))  
    ax_val.grid(True, alpha=0.3)
    ax_val.legend(loc='upper right', fontsize=10)
    
    if col == 0: 
        ax_val.set_ylabel("Sampled Value", fontsize=13)
        
        # Add a secondary axis for the sample size bar chart using exact bin intervals
        ax_twin = ax_val.twinx()
        ax_twin.bar(mae_df['Bin_Center'], mae_df['Sample_Size'], width=4.5, color='gray', alpha=0.25, zorder=0)
        ax_twin.set_ylim(0, 75) # Scale y-axis up to 75 so the max bars stretch closer to the legend
        ax_twin.set_yticks([0, 10, 20, 30, 40])
        # Place ticks inside to prevent layout collision; push labels slightly inside too
        ax_twin.tick_params(axis='y', direction='in', pad=-15, labelsize=9)
        for label in ax_twin.get_yticklabels():
            label.set_horizontalalignment('right')
            
        ax_twin.set_ylabel("Bin Sample Size", fontsize=10, rotation=270, labelpad=15)
            
        # Ensure bars stay behind the main plot lines
        ax_val.set_zorder(ax_twin.get_zorder() + 1)
        ax_val.patch.set_visible(False)
    
    ax_mae.axvline(50, color='gray', linestyle=':', linewidth=2)
    ax_mae.set_xlim(0, 100)
    ax_mae.set_xticks(np.arange(0, 81, 20))  
    ax_mae.grid(True, alpha=0.3)
    if col == 0: ax_mae.set_ylabel("MAE", fontsize=13)
        
    ax_mre.axvline(50, color='gray', linestyle=':', linewidth=2)
    ax_mre.set_xlim(0, 100)
    ax_mre.set_xticks(np.arange(0, 81, 20))  
    ax_mre.grid(True, alpha=0.3)
    if col == 0: ax_mre.set_ylabel("MRE (%)", fontsize=13)
        
    ax_bias.axvline(50, color='gray', linestyle=':', linewidth=2)
    ax_bias.axhline(0, color='gray', linestyle='--', linewidth=1.5)
    ax_bias.set_xlim(0, 100)  
    ax_bias.set_xticks(np.arange(0, 81, 20))
    ax_bias.grid(True, alpha=0.3)
    ax_bias.set_xlabel("True Bare Ground (%)", fontsize=12)
    if col == 0: ax_bias.set_ylabel("Mean Bias", fontsize=13)

# ====================================================================
# PLOTTING - FIGURE 2: UNIQUE 2x5 DRONE VALUES GRID
# ====================================================================
plot_config_2 = [
    ('BG', 'Total Bare\nGround (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
    ('Herb', 'Herbaceous\nCover (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
    ('Woody', 'Woody\nCover (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
    ('Fetch', 'Mean\nFetch (m)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
    ('LPI', 'Largest Patch\nIndex (%)', [], 'black'), 
    ('Gap_0_24', 'Canopy Gap\n0-24 cm (%)', ['0cm'], 'steelblue'),
    ('Gap_25_50', 'Canopy Gap\n25-50 cm (%)', ['0cm'], 'cadetblue'),
    ('Gap_51_100', 'Canopy Gap\n51-100 cm (%)', ['0cm'], 'mediumseagreen'),
    ('Gap_101_200', 'Canopy Gap\n101-200 cm (%)', ['0cm'], 'darkorange'),
    ('Gap_gt_200', 'Canopy Gap\n+200 cm (%)', ['0cm'], 'firebrick')
]

plt.rcParams.update({
    'figure.constrained_layout.use': True,
    'figure.constrained_layout.w_pad': 0.02,
    'figure.constrained_layout.h_pad': 0.02,
    'figure.constrained_layout.wspace': 0.02,
    'figure.constrained_layout.hspace': 0.02,
})

fig2, axes2 = plt.subplots(2, 5, figsize=(11, 5))
fig2.suptitle("Drone Imagery Metrics vs True Bare Ground (%)", fontsize=20, weight='bold')

for idx, (prefix, title, scales, col_color) in enumerate(plot_config_2):
    row = idx // 5
    col = idx % 5
    ax = axes2[row, col]
    
    letter = chr(ord('a') + idx)
    
    if prefix == 'BG': true_col = 'True_BG_Mean'
    elif prefix == 'Fetch': true_col = 'True_Fetch_Mean'
    elif prefix == 'LPI': true_col = 'True_LPI_Mean'
    elif prefix == 'Herb': true_col = 'True_Herb_Mean'
    elif prefix == 'Woody': true_col = 'True_Woody_Mean'
    elif prefix.startswith('Gap'): true_col = f'True_{prefix}_Mean'

    if true_col in mae_df.columns:
        linestyle = '-' if prefix == 'LPI' else '--'
        linecolor = 'black' if prefix == 'LPI' else 'gray'
        ax.plot(mae_df['True_BG_Mean'], mae_df[true_col], color=linecolor, linestyle=linestyle, linewidth=2, label='Exact Value', zorder=1)
        
    for scale in scales:
        var_base = f'{prefix}_{scale}'
        if f'{var_base}_Val' not in mae_df.columns:
            continue
            
        color = scale_colors[scale] if len(scales) > 1 else col_color
        label = f'{scale}' if scale != '0cm' else '0cm Cont.'
        line_kws = {'marker': 'o', 'color': color, 'linestyle': '-', 'linewidth': 2.5, 'markersize': 6}
        
        ax.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Val'], label=label, **line_kws)
        
    autoscale_y_robust(ax, force_zero=False)
        
    title_str = r"$\mathbf{" + letter + r".}$ " + title
    ax.set_title(title_str, fontsize=14, pad=8) 
    ax.axvline(50, color='gray', linestyle=':', linewidth=2)
    ax.set_xlim(0, 100)
    ax.set_xticks(np.arange(0, 81, 20)) 
    
    ax.tick_params(axis='both', which='major', labelsize=7.5) 
    ax.grid(True, alpha=0.3)
    
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc='upper right', fontsize=7.5, handlelength=1.0, borderpad=0.2, labelspacing=0.1, handletextpad=0.4, framealpha=0.7) 
            
    if row == 1:
        if col == 2:
            ax.set_xlabel("True Bare Ground (%)", fontsize=15, weight='bold') 
        else:
            ax.set_xlabel("")

fig2.supylabel("Fraction % (or, m for Mean Fetch)", fontsize=14, weight='bold')

# ====================================================================
# PLOTTING - FIGURE 3: COVER METRICS (BG, HERB, WOODY)
# ====================================================================
plt.rcdefaults() # Reset rcParams to avoid constrained_layout squishing from Figure 2

fig3_config = [
    ('BG', 'Total Bare Ground (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
    ('Herb', 'Herbaceous Cover (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
    ('Woody', 'Woody Cover (%)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None)
]

fig3, axes3 = plt.subplots(3, 4, figsize=(16, 12), constrained_layout=True)
fig3.suptitle("Drone Imagery Cover Metrics Across Bare Ground Gradient", fontsize=20, weight='bold')

for row, (prefix, title, scales, col_color) in enumerate(fig3_config):
    ax_val = axes3[row, 0]
    ax_mae = axes3[row, 1]
    ax_mre = axes3[row, 2]
    ax_bias = axes3[row, 3]

    if prefix == 'BG':
        ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_BG_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
    elif prefix == 'Herb' and 'True_Herb_Mean' in mae_df.columns:
        ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_Herb_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
    elif prefix == 'Woody' and 'True_Woody_Mean' in mae_df.columns:
        ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_Woody_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)

    for scale in scales:
        var_base = f'{prefix}_{scale}'
        if f'{var_base}_Val' not in mae_df.columns:
            continue
            
        color = scale_colors[scale] if len(scales) > 1 else col_color
        label_name = f'{scale} Point' if scale != '0cm' else '0cm Continuous'
        line_kws = {'marker': 'o', 'color': color, 'linestyle': '-', 'linewidth': 2.5, 'markersize': 7}
        
        ax_val.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Val'], label=label_name, **line_kws)
        ax_mae.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MAE'], **line_kws)
        ax_mre.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MRE'], **line_kws)
        ax_bias.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Bias'], **line_kws)

    autoscale_y_robust(ax_val, force_zero=False)
    autoscale_y_robust(ax_mae, force_zero=False)
    autoscale_y_robust(ax_mre, force_zero=False)
    autoscale_y_robust(ax_bias, force_zero=True) 

    for ax in [ax_val, ax_mae, ax_mre, ax_bias]:
        ax.set_xlim(0, 100)
        ax.set_xticks(np.arange(0, 81, 20)) 
        ax.axvline(50, color='gray', linestyle=':', linewidth=2)
        
    ax_val.set_ylabel(title, fontsize=13, weight='bold') 
    ax_val.grid(True, alpha=0.3)
    
    if row == 0:
        ax_val.legend(loc='upper right', fontsize=10)
        
        # Add a secondary axis for the sample size bar chart using exact bin intervals
        ax_twin = ax_val.twinx()
        ax_twin.bar(mae_df['Bin_Center'], mae_df['Sample_Size'], width=4.5, color='gray', alpha=0.25, zorder=0)
        ax_twin.set_ylim(0, 75) # Scale y-axis up to 75 so the max bars stretch closer to the legend
        ax_twin.set_yticks([0, 10, 20, 30, 40])
        # Place ticks inside to prevent layout collision; push labels slightly inside too
        ax_twin.tick_params(axis='y', direction='in', pad=-15, labelsize=9)
        for label in ax_twin.get_yticklabels():
            label.set_horizontalalignment('right')
            
        ax_twin.set_ylabel("Bin Sample Size", fontsize=10, rotation=270, labelpad=15)
            
        # Ensure bars stay behind the main plot lines
        ax_val.set_zorder(ax_twin.get_zorder() + 1)
        ax_val.patch.set_visible(False)

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
# PLOTTING - FIGURE 4: SPATIAL METRICS (FETCH, GAPS)
# ====================================================================
fig4_config = [
    ('Fetch', 'Mean Fetch (m)', ['0cm', '25cm', '50cm', '100cm', '200cm'], None),
    ('Gap_0_24', 'Canopy Gap\n0-24cm (%)', ['0cm'], 'steelblue'),
    ('Gap_25_50', 'Canopy Gap\n25-50cm (%)', ['0cm'], 'cadetblue'),
    ('Gap_51_100', 'Canopy Gap\n51-100cm (%)', ['0cm'], 'mediumseagreen'),
    ('Gap_101_200', 'Canopy Gap\n101-200cm (%)', ['0cm'], 'darkorange'),
    ('Gap_gt_200', 'Canopy Gap\n>200cm (%)', ['0cm'], 'firebrick')
]

fig4, axes4 = plt.subplots(6, 4, figsize=(16, 24), constrained_layout=True)
fig4.suptitle("Drone Imagery Spatial Metrics Across Bare Ground Gradient", fontsize=24, weight='bold')

for row, (prefix, title, scales, col_color) in enumerate(fig4_config):
    ax_val = axes4[row, 0]
    ax_mae = axes4[row, 1]
    ax_mre = axes4[row, 2]
    ax_bias = axes4[row, 3]

    if prefix == 'Fetch':
        ax_val.plot(mae_df['True_BG_Mean'], mae_df['True_Fetch_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)
    elif prefix.startswith('Gap'):
        ax_val.plot(mae_df['True_BG_Mean'], mae_df[f'True_{prefix}_Mean'], color='gray', linestyle='--', linewidth=2, label='Exact True Value', zorder=1)

    for scale in scales:
        var_base = f'{prefix}_{scale}'
        if f'{var_base}_Val' not in mae_df.columns:
            continue
            
        color = scale_colors[scale] if len(scales) > 1 else col_color
        label_name = f'{scale} Point' if scale != '0cm' else '0cm Continuous'
        line_kws = {'marker': 'o', 'color': color, 'linestyle': '-', 'linewidth': 2.5, 'markersize': 7}
        
        ax_val.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Val'], label=label_name, **line_kws)
        ax_mae.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MAE'], **line_kws)
        ax_mre.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_MRE'], **line_kws)
        ax_bias.plot(mae_df['True_BG_Mean'], mae_df[f'{var_base}_Bias'], **line_kws)

    autoscale_y_robust(ax_val, force_zero=False)
    autoscale_y_robust(ax_mae, force_zero=False)
    autoscale_y_robust(ax_mre, force_zero=False)
    autoscale_y_robust(ax_bias, force_zero=True) 

    for ax in [ax_val, ax_mae, ax_mre, ax_bias]:
        ax.set_xlim(0, 100)
        ax.set_xticks(np.arange(0, 81, 20)) 
        ax.axvline(50, color='gray', linestyle=':', linewidth=2)
        ax.tick_params(axis='both', labelsize=12)
        
    ax_val.set_ylabel(title, fontsize=16, weight='bold') 
    ax_val.grid(True, alpha=0.3)
    ax_val.legend(loc='upper right', fontsize=12)

    ax_mae.grid(True, alpha=0.3)
    ax_mre.grid(True, alpha=0.3)

    ax_bias.axhline(0, color='gray', linestyle='--', linewidth=1.5)
    ax_bias.grid(True, alpha=0.3)

    if row == 0:
        ax_val.set_title("Sampled Value", fontsize=18, pad=12)
        ax_mae.set_title("MAE", fontsize=18, pad=12)
        ax_mre.set_title("MRE (%)", fontsize=18, pad=12)
        ax_bias.set_title("Mean Bias", fontsize=18, pad=12)

    if row == 5:
        for ax in [ax_val, ax_mae, ax_mre, ax_bias]:
            ax.set_xlabel("True Bare Ground (%)", fontsize=15)

# ====================================================================
# 6. EXPORT
# ====================================================================
img_path1 = os.path.join(base_dir, 'Drone_Metrics_MultiScale.svg')
img_path2 = os.path.join(base_dir, 'Drone_Metrics_Values_2x5.svg')
img_path3 = os.path.join(base_dir, 'Drone_Cover_Metrics.svg')
img_path4 = os.path.join(base_dir, 'Drone_Spatial_Metrics.svg')

fig1.savefig(img_path1, format='svg', dpi=300, bbox_inches='tight')
fig2.savefig(img_path2, format='svg', dpi=300, bbox_inches='tight')
fig3.savefig(img_path3, format='svg', dpi=300, bbox_inches='tight')
fig4.savefig(img_path4, format='svg', dpi=300, bbox_inches='tight')

print("All figures successfully saved as SVG.")
plt.show()
