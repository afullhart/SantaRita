import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load the raw dataset
df = pd.read_csv(r"C:\Users\andre\ScatterPlots\SRER_NRI_Plots_110m.csv")

# 2. Define the metrics mapping to calculate the stats dynamically
metrics_info = {
    'BGR': {'exact': 'Exact_BGR_Pct', 'nri_pattern': 'NRI_BGR_{scale}_Pct'},
    'Herb': {'exact': 'Exact_Herb_Pct', 'nri_pattern': 'NRI_Herb_{scale}_Pct'},
    'Woody': {'exact': 'Exact_Woody_Pct', 'nri_pattern': 'NRI_Woody_{scale}_Pct'},
    'Herb_Woody_Ratio': {'exact': 'Exact_Herb_Woody_Ratio', 'nri_pattern': 'NRI_HW_Ratio_{scale}'},
    'Fetch': {'exact': 'Exact_Fetch_m', 'nri_pattern': 'NRI_Fetch_{scale}'},
    'Gap_0_24': {'exact': 'Exact_Gap_0_24', 'nri_pattern': 'NRI_Gap_0_24', 'is_gap': True},
    'Gap_25_50': {'exact': 'Exact_Gap_25_50', 'nri_pattern': 'NRI_Gap_25_50', 'is_gap': True},
    'Gap_51_100': {'exact': 'Exact_Gap_51_100', 'nri_pattern': 'NRI_Gap_51_100', 'is_gap': True},
    'Gap_101_200': {'exact': 'Exact_Gap_101_200', 'nri_pattern': 'NRI_Gap_101_200', 'is_gap': True},
    'Gap_gt_200': {'exact': 'Exact_Gap_gt_200', 'nri_pattern': 'NRI_Gap_gt_200', 'is_gap': True},
}

scales = ['0cm', '25cm', '50cm', '100cm', '200cm']
results = []

# 3. Calculate Pearson r and Sen's Slope for all valid pairs
for metric, info in metrics_info.items():
    exact_col = info['exact']
    
    if info.get('is_gap'):
        # Gaps use only the 0cm baseline (No specific scale suffix)
        nri_col = info['nri_pattern']
        if exact_col in df.columns and nri_col in df.columns:
            # Filter out NaNs for the calculation
            valid_idx = df[exact_col].notna() & df[nri_col].notna()
            x = df.loc[valid_idx, exact_col]
            y = df.loc[valid_idx, nri_col]
            
            if len(x) > 1: # Need at least 2 points to calculate correlation
                r, _ = stats.pearsonr(x, y)
                sens, _, _, _ = stats.theilslopes(y, x)
                results.append({'Metric': metric, 'Scale': '0cm', 'r': r, 'Sens_Slope': sens})
    else:
        # Check across all specified scales (0cm, 25cm, 50cm, etc.)
        for scale in scales:
            nri_col = info['nri_pattern'].format(scale=scale)
            if exact_col in df.columns and nri_col in df.columns:
                valid_idx = df[exact_col].notna() & df[nri_col].notna()
                x = df.loc[valid_idx, exact_col]
                y = df.loc[valid_idx, nri_col]
                
                if len(x) > 1:
                    r, _ = stats.pearsonr(x, y)
                    sens, _, _, _ = stats.theilslopes(y, x)
                    results.append({'Metric': metric, 'Scale': scale, 'r': r, 'Sens_Slope': sens})

# Create a DataFrame from the calculated results
res_df = pd.DataFrame(results)

# 4. Pivot the data
r_pivot = res_df.pivot(index='Metric', columns='Scale', values='r')
sens_pivot = res_df.pivot(index='Metric', columns='Scale', values='Sens_Slope')

# 5. Reorder indices to match your desired plotting order
metric_order = ['BGR', 'Herb', 'Woody', 'Herb_Woody_Ratio', 'Fetch', 'Gap_0_24', 'Gap_25_50', 'Gap_51_100', 'Gap_101_200', 'Gap_gt_200']
scale_order = ['0cm', '25cm', '50cm', '100cm', '200cm']

r_pivot = r_pivot.reindex(index=[m for m in metric_order if m in r_pivot.index], columns=scale_order)
sens_pivot = sens_pivot.reindex(index=[m for m in metric_order if m in sens_pivot.index], columns=scale_order)


# ====================================================================
# HELPER TO DRAW STRIKES THROUGH BLANK CELLS
# ====================================================================
def add_strikes(ax, data_pivot):
    """Draws a diagonal strike through cells that contain NaN values."""
    for i in range(data_pivot.shape[0]):
        for j in range(data_pivot.shape[1]):
            if pd.isna(data_pivot.iloc[i, j]):
                # Draw a diagonal line from bottom-left to top-right of the cell
                ax.plot([j, j+1], [i+1, i], color='gray', lw=1.5)


# ====================================================================
# PLOTTING
# ====================================================================
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Left Plot: Pearson r
sns.heatmap(r_pivot, annot=True, fmt=".2f", cmap="Blues", ax=axes[0], 
            cbar_kws={'label': 'Pearson $r$'}, mask=r_pivot.isnull()) 
add_strikes(axes[0], r_pivot)  
axes[0].set_title('Pearson Correlation ($r$)\nExact vs Sampled', pad=15, fontweight='bold')
axes[0].set_ylabel('Ground Cover Metric', fontweight='bold')
axes[0].set_xlabel('NRI Transect Sampling Scale', fontweight='bold')

# Right Plot: Sen's Slope
sns.heatmap(sens_pivot, annot=True, fmt=".2f", cmap="vlag", center=1.0, ax=axes[1], 
            cbar_kws={'label': "Sen's Slope"}, mask=sens_pivot.isnull()) 
add_strikes(axes[1], sens_pivot)  
axes[1].set_title("Sen's Slope\nExact vs Sampled (1.0 = Perfect 1:1)", pad=15, fontweight='bold')
axes[1].set_ylabel('')
axes[1].set_xlabel('NRI Transect Sampling Scale', fontweight='bold')

plt.tight_layout()
plt.savefig(r'C:\Users\andre\ScatterPlots\Exact_vs_Sampled_Heatmaps.png', dpi=300, bbox_inches='tight')
plt.show()






import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# Suppress runtime warnings for division by zero (handled intentionally with NaNs)
warnings.filterwarnings('ignore', category=RuntimeWarning)

# 1. Load the raw dataset
df = pd.read_csv(r"C:\Users\andre\ScatterPlots\SRER_NRI_Plots_110m.csv")

# 2. Define the metrics mapping to calculate the stats dynamically
metrics_info = {
    'BGR': {'exact': 'Exact_BGR_Pct', 'nri_pattern': 'NRI_BGR_{scale}_Pct'},
    'Herb': {'exact': 'Exact_Herb_Pct', 'nri_pattern': 'NRI_Herb_{scale}_Pct'},
    'Woody': {'exact': 'Exact_Woody_Pct', 'nri_pattern': 'NRI_Woody_{scale}_Pct'},
    'Herb_Woody_Ratio': {'exact': 'Exact_Herb_Woody_Ratio', 'nri_pattern': 'NRI_HW_Ratio_{scale}'},
    'Fetch': {'exact': 'Exact_Fetch_m', 'nri_pattern': 'NRI_Fetch_{scale}'},
    'Gap_0_24': {'exact': 'Exact_Gap_0_24', 'nri_pattern': 'NRI_Gap_0_24', 'is_gap': True},
    'Gap_25_50': {'exact': 'Exact_Gap_25_50', 'nri_pattern': 'NRI_Gap_25_50', 'is_gap': True},
    'Gap_51_100': {'exact': 'Exact_Gap_51_100', 'nri_pattern': 'NRI_Gap_51_100', 'is_gap': True},
    'Gap_101_200': {'exact': 'Exact_Gap_101_200', 'nri_pattern': 'NRI_Gap_101_200', 'is_gap': True},
    'Gap_gt_200': {'exact': 'Exact_Gap_gt_200', 'nri_pattern': 'NRI_Gap_gt_200', 'is_gap': True},
}

scales = ['0cm', '25cm', '50cm', '100cm', '200cm']
results = []

# 3. Calculate MAE and MRE for all valid pairs
for metric, info in metrics_info.items():
    exact_col = info['exact']
    
    if info.get('is_gap'):
        # Gaps use only the 0cm baseline (No specific scale suffix)
        nri_col = info['nri_pattern']
        if exact_col in df.columns and nri_col in df.columns:
            valid_idx = df[exact_col].notna() & df[nri_col].notna()
            x = df.loc[valid_idx, exact_col]
            y = df.loc[valid_idx, nri_col]
            
            if len(x) > 0: 
                # Mean Absolute Error
                mae = np.mean(np.abs(y - x))
                # Mean Relative Error (multiplied by 100 for percentage)
                mre = np.nanmean(np.abs(y - x) / np.where(x == 0, np.nan, np.abs(x))) * 100
                
                results.append({'Metric': metric, 'Scale': '0cm', 'MAE': mae, 'MRE': mre})
    else:
        # Check across all specified scales (0cm, 25cm, 50cm, etc.)
        for scale in scales:
            nri_col = info['nri_pattern'].format(scale=scale)
            if exact_col in df.columns and nri_col in df.columns:
                valid_idx = df[exact_col].notna() & df[nri_col].notna()
                x = df.loc[valid_idx, exact_col]
                y = df.loc[valid_idx, nri_col]
                
                if len(x) > 0:
                    mae = np.mean(np.abs(y - x))
                    # Mean Relative Error (multiplied by 100 for percentage)
                    mre = np.nanmean(np.abs(y - x) / np.where(x == 0, np.nan, np.abs(x))) * 100
                    
                    results.append({'Metric': metric, 'Scale': scale, 'MAE': mae, 'MRE': mre})

# Create a DataFrame from the calculated results
res_df = pd.DataFrame(results)

# 4. Pivot the data for the heatmaps
mae_pivot = res_df.pivot(index='Metric', columns='Scale', values='MAE')
mre_pivot = res_df.pivot(index='Metric', columns='Scale', values='MRE')

# 5. Reorder indices to match desired plotting order
metric_order = ['BGR', 'Herb', 'Woody', 'Herb_Woody_Ratio', 'Fetch', 'Gap_0_24', 'Gap_25_50', 'Gap_51_100', 'Gap_101_200', 'Gap_gt_200']
scale_order = ['0cm', '25cm', '50cm', '100cm', '200cm']

mae_pivot = mae_pivot.reindex(index=[m for m in metric_order if m in mae_pivot.index], columns=scale_order)
mre_pivot = mre_pivot.reindex(index=[m for m in metric_order if m in mre_pivot.index], columns=scale_order)

# ====================================================================
# HELPER TO DRAW STRIKES THROUGH BLANK CELLS
# ====================================================================
def add_strikes(ax, data_pivot):
    """Draws a diagonal strike through cells that contain NaN values."""
    for i in range(data_pivot.shape[0]):
        for j in range(data_pivot.shape[1]):
            if pd.isna(data_pivot.iloc[i, j]):
                # Draw a diagonal line from bottom-left to top-right of the cell
                ax.plot([j, j+1], [i+1, i], color='gray', lw=1.5)

# ====================================================================
# PLOTTING
# ====================================================================
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Left Plot: Mean Absolute Error (MAE)
sns.heatmap(mae_pivot, annot=True, fmt=".2f", cmap="Reds", ax=axes[0], 
            cbar_kws={'label': 'Mean Absolute Error'}, mask=mae_pivot.isnull()) 
add_strikes(axes[0], mae_pivot)  
axes[0].set_title('Mean Absolute Error (MAE)\nExact vs Sampled', pad=15, fontweight='bold')
axes[0].set_ylabel('Ground Cover Metric', fontweight='bold')
axes[0].set_xlabel('NRI Transect Sampling Scale', fontweight='bold')

# Right Plot: Mean Relative Error (MRE) - Now as a Percentage
sns.heatmap(mre_pivot, annot=True, fmt=".2f", cmap="Oranges", ax=axes[1], 
            cbar_kws={'label': "Mean Relative Error (%)"}, mask=mre_pivot.isnull()) 
add_strikes(axes[1], mre_pivot)  
axes[1].set_title("Mean Relative Error (MRE)\nExact vs Sampled (% of Exact)", pad=15, fontweight='bold')
axes[1].set_ylabel('')
axes[1].set_xlabel('NRI Transect Sampling Scale', fontweight='bold')

plt.tight_layout()
plt.savefig(r'C:\Users\andre\ScatterPlots\Error_Metrics_Heatmaps.png', dpi=300, bbox_inches='tight')
plt.show()
