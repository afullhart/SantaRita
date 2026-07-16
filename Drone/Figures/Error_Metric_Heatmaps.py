import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# Suppress runtime warnings for division by zero (handled intentionally with NaNs)
warnings.filterwarnings('ignore', category=RuntimeWarning)

# ====================================================================
# 1. LOAD DATA & DEFINE MAPPINGS
# ====================================================================
df = pd.read_csv(r"C:\Users\andre\ScatterPlots\SRER_NRI_Plots_110m.csv")

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

# ====================================================================
# 2. CALCULATE METRICS IN A SINGLE PASS
# ====================================================================
def calculate_stats(scl, x_vals, y_vals):
    """Calculates Pearson r, Sen's Slope, MAE, and MRE"""
    if len(x_vals) > 1:
        r, _ = stats.pearsonr(x_vals, y_vals)
        sens, _, _, _ = stats.theilslopes(y_vals, x_vals)
        mae = np.mean(np.abs(y_vals - x_vals))
        mre = np.nanmean(np.abs(y_vals - x_vals) / np.where(x_vals == 0, np.nan, np.abs(x_vals))) * 100
        return {'Metric': metric, 'Scale': scl, 'r': r, 'Sens_Slope': sens, 'MAE': mae, 'MRE': mre}
    return None

for metric, info in metrics_info.items():
    exact_col = info['exact']
    
    if info.get('is_gap'):
        # Gaps use only the 0cm baseline (No specific scale suffix)
        nri_col = info['nri_pattern']
        if exact_col in df.columns and nri_col in df.columns:
            valid_idx = df[exact_col].notna() & df[nri_col].notna()
            x = df.loc[valid_idx, exact_col]
            y = df.loc[valid_idx, nri_col]
            
            res = calculate_stats('0cm', x, y)
            if res: results.append(res)
    else:
        # Check across all specified scales (0cm, 25cm, 50cm, etc.)
        for scale in scales:
            nri_col = info['nri_pattern'].format(scale=scale)
            if exact_col in df.columns and nri_col in df.columns:
                valid_idx = df[exact_col].notna() & df[nri_col].notna()
                x = df.loc[valid_idx, exact_col]
                y = df.loc[valid_idx, nri_col]
                
                res = calculate_stats(scale, x, y)
                if res: results.append(res)

res_df = pd.DataFrame(results)

# ====================================================================
# 3. PIVOT & REORDER DATA
# ====================================================================
r_pivot = res_df.pivot(index='Metric', columns='Scale', values='r')
sens_pivot = res_df.pivot(index='Metric', columns='Scale', values='Sens_Slope')
mae_pivot = res_df.pivot(index='Metric', columns='Scale', values='MAE')
mre_pivot = res_df.pivot(index='Metric', columns='Scale', values='MRE')

metric_order = ['BGR', 'Herb', 'Woody', 'Herb_Woody_Ratio', 'Fetch', 'Gap_0_24', 'Gap_25_50', 'Gap_51_100', 'Gap_101_200', 'Gap_gt_200']
scale_order = ['0cm', '25cm', '50cm', '100cm', '200cm']

# Safe reindexing to maintain desired visual plotting order
pivots = [r_pivot, sens_pivot, mae_pivot, mre_pivot]
for i in range(len(pivots)):
    pivots[i] = pivots[i].reindex(index=[m for m in metric_order if m in pivots[i].index], columns=scale_order)
r_pivot, sens_pivot, mae_pivot, mre_pivot = pivots

# ====================================================================
# 4. HELPER TO DRAW STRIKES THROUGH BLANK CELLS
# ====================================================================
def add_strikes(ax, data_pivot):
    """Draws a diagonal strike through cells that contain NaN values."""
    for i in range(data_pivot.shape[0]):
        for j in range(data_pivot.shape[1]):
            if pd.isna(data_pivot.iloc[i, j]):
                # Draw a diagonal line from bottom-left to top-right of the cell
                ax.plot([j, j+1], [i+1, i], color='gray', lw=1.5)

# ====================================================================
# 5. PLOTTING - FIGURE 1: Pearson r & Sen's Slope
# ====================================================================
fig1, axes1 = plt.subplots(1, 2, figsize=(16, 6))

# Left Plot: Pearson r
sns.heatmap(r_pivot, annot=True, fmt=".2f", cmap="Blues", ax=axes1[0], 
            cbar_kws={'label': 'Pearson $r$'}, mask=r_pivot.isnull()) 
add_strikes(axes1[0], r_pivot)  
axes1[0].set_title('Pearson Correlation ($r$)\nExact vs Sampled', pad=15, fontweight='bold')
axes1[0].set_ylabel('Ground Cover Metric', fontweight='bold')
axes1[0].set_xlabel('NRI Transect Sampling Scale', fontweight='bold')

# Right Plot: Sen's Slope
sns.heatmap(sens_pivot, annot=True, fmt=".2f", cmap="vlag", center=1.0, ax=axes1[1], 
            cbar_kws={'label': "Sen's Slope"}, mask=sens_pivot.isnull()) 
add_strikes(axes1[1], sens_pivot)  
axes1[1].set_title("Sen's Slope\nExact vs Sampled (1.0 = Perfect 1:1)", pad=15, fontweight='bold')
axes1[1].set_ylabel('')
axes1[1].set_xlabel('NRI Transect Sampling Scale', fontweight='bold')

plt.tight_layout()
plt.savefig(r'C:\Users\andre\ScatterPlots\Exact_vs_Sampled_Heatmaps.png', dpi=300, bbox_inches='tight')
plt.show()

# ====================================================================
# 6. PLOTTING - FIGURE 2: MAE & MRE
# ====================================================================
fig2, axes2 = plt.subplots(1, 2, figsize=(16, 6))

# Left Plot: Mean Absolute Error (MAE)
sns.heatmap(mae_pivot, annot=True, fmt=".2f", cmap="Reds", ax=axes2[0], 
            cbar_kws={'label': 'Mean Absolute Error'}, mask=mae_pivot.isnull()) 
add_strikes(axes2[0], mae_pivot)  
axes2[0].set_title('Mean Absolute Error (MAE)\nExact vs Sampled', pad=15, fontweight='bold')
axes2[0].set_ylabel('Ground Cover Metric', fontweight='bold')
axes2[0].set_xlabel('NRI Transect Sampling Scale', fontweight='bold')

# Right Plot: Mean Relative Error (MRE)
sns.heatmap(mre_pivot, annot=True, fmt=".2f", cmap="Oranges", ax=axes2[1], 
            cbar_kws={'label': "Mean Relative Error (%)"}, mask=mre_pivot.isnull()) 
add_strikes(axes2[1], mre_pivot)  
axes2[1].set_title("Mean Relative Error (MRE)\nExact vs Sampled (% of Exact)", pad=15, fontweight='bold')
axes2[1].set_ylabel('')
axes2[1].set_xlabel('NRI Transect Sampling Scale', fontweight='bold')

plt.tight_layout()
plt.savefig(r'C:\Users\andre\ScatterPlots\Error_Metrics_Heatmaps.png', dpi=300, bbox_inches='tight')
plt.show()
