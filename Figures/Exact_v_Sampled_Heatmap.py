import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load your pre-calculated exact vs sampled metrics
df = pd.read_csv(r"C:\Users\andre\ScatterPlots\Exact_vs_Sampled_Metrics.csv")

# 2. Corrected Helper functions
def get_scale(y_col):
    if 'Gap' in y_col: 
        return '0cm'  # All gaps strictly use the 0cm baseline sampling scale
    elif '25cm' in y_col: 
        return '25cm'
    elif '50cm' in y_col: 
        return '50cm'
    elif '100cm' in y_col: 
        return '100cm'
    elif '200cm' in y_col: 
        return '200cm'
    elif '0cm' in y_col: 
        return '0cm'
    return 'Unknown'

def get_base_metric(x_col):
    return x_col.replace('Exact_', '').replace('_Pct', '').replace('_m', '')

df['Scale'] = df['Y_Column'].apply(get_scale)
df['Metric'] = df['X_Column'].apply(get_base_metric)

# 3. Calculate Pearson r from R2
df['r'] = np.sqrt(df['R2']) * np.sign(df['OLS_Slope'])

# 4. Pivot the data
r_pivot = df.pivot(index='Metric', columns='Scale', values='r')
sens_pivot = df.pivot(index='Metric', columns='Scale', values='Sens_Slope')

# 5. Reorder indices (Includes the new 25cm and 200cm scales)
metric_order = ['BGR', 'Herb', 'Woody', 'Herb_Woody_Ratio', 'Fetch', 'Gap_0_24', 'Gap_25_50', 'Gap_51_100', 'Gap_101_200', 'Gap_gt_200']
scale_order = ['0cm', '25cm', '50cm', '100cm', '200cm']

r_pivot = r_pivot.reindex(index=[m for m in metric_order if m in r_pivot.index], 
                          columns=[s for s in scale_order if s in r_pivot.columns])
sens_pivot = sens_pivot.reindex(index=[m for m in metric_order if m in sens_pivot.index], 
                                columns=[s for s in scale_order if s in sens_pivot.columns])


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
            cbar_kws={'label': 'Pearson $r$'},
            mask=r_pivot.isnull()) 
add_strikes(axes[0], r_pivot)  # Apply strikes to empty cells
axes[0].set_title('Pearson Correlation ($r$)\nExact vs Sampled', pad=15, fontweight='bold')
axes[0].set_ylabel('Ground Cover Metric', fontweight='bold')
axes[0].set_xlabel('NRI Transect Sampling Scale', fontweight='bold') # Updated label

# Right Plot: Sen's Slope
sns.heatmap(sens_pivot, annot=True, fmt=".2f", cmap="vlag", center=1.0, ax=axes[1], 
            cbar_kws={'label': "Sen's Slope"},
            mask=sens_pivot.isnull()) 
add_strikes(axes[1], sens_pivot)  # Apply strikes to empty cells
axes[1].set_title("Sen's Slope\nExact vs Sampled (1.0 = Perfect 1:1)", pad=15, fontweight='bold')
axes[1].set_ylabel('')
axes[1].set_xlabel('NRI Transect Sampling Scale', fontweight='bold') # Updated label

plt.tight_layout()
plt.savefig(r'C:\Users\andre\ScatterPlots\Exact_vs_Sampled_Heatmaps.png', dpi=300, bbox_inches='tight')
plt.show()
