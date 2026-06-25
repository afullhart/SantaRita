import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ====================================================================
# CONFIGURATION
# ====================================================================
output_dir = r"C:\Users\andre\ScatterPlots"
os.makedirs(output_dir, exist_ok=True)

# Pointing to your raw 110m measurement data
csv_path = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\SRER_NRI_Plots_110m.csv"
df = pd.read_csv(csv_path)

# ====================================================================
# VARIABLE SELECTION & LABELING
# ====================================================================
# The exact metric columns from your raw dataset (Herb/Woody moved to end)
target_cols = [
    'Exact_BGR_Pct', 
    'Exact_Fetch_m', 
    'Exact_LPI_Pct',
    'Exact_Gap_0_24', 
    'Exact_Gap_25_50', 
    'Exact_Gap_51_100',
    'Exact_Gap_101_200', 
    'Exact_Gap_gt_200',
    'Exact_Herb_Pct', 
    'Exact_Woody_Pct', 
    'Exact_Herb_Woody_Ratio'
]

# Clean labels for the plot axes (Herb/Woody moved to end)
clean_labels = [
    'Bare Ground %', 
    'Mean Fetch (m)', 
    'LPI %',
    'Gap 0-24 cm', 
    'Gap 25-50 cm', 
    'Gap 51-100 cm',
    'Gap 101-200 cm', 
    'Gap > 200 cm',
    'Herb Cover %', 
    'Woody Cover %', 
    'Herb:Woody Ratio'
]

# ====================================================================
# CORRELATION CALCULATION & PLOTTING
# ====================================================================
print("Calculating split-triangle Pearson & Spearman correlation matrix...")

# Subset the dataframe to only the target columns
df_subset = df[target_cols]

# Calculate BOTH correlation matrices
corr_pearson = df_subset.corr(method='pearson')
corr_spearman = df_subset.corr(method='spearman')

# Create an empty matrix to hold the combined values (explicitly float to handle NaNs)
corr_combined = np.zeros_like(corr_pearson.values, dtype=float)

# Get the indices for the lower and upper triangles (offset by 1 to exclude diagonal)
lower_indices = np.tril_indices(len(target_cols), k=-1)
upper_indices = np.triu_indices(len(target_cols), k=1)

# Populate the empty matrix: Pearson in lower left, Spearman in upper right
corr_combined[lower_indices] = corr_pearson.values[lower_indices]
corr_combined[upper_indices] = corr_spearman.values[upper_indices]

# FIX: Set the diagonal to NaN BEFORE creating the DataFrame.
# Seaborn will leave NaNs transparent, allowing our black background to show through naturally.
np.fill_diagonal(corr_combined, np.nan)

# Convert back to a DataFrame for Seaborn
corr_combined_df = pd.DataFrame(corr_combined, index=corr_pearson.index, columns=corr_pearson.columns)

# Set up the matplotlib figure
fig, ax = plt.subplots(figsize=(14, 12))

# Paint the background of the axes black to highlight the diagonal separation
ax.set_facecolor('black')

# Generate a custom diverging colormap (blue to red)
cmap = sns.color_palette("vlag", as_cmap=True)

# Draw the heatmap using seaborn
sns.heatmap(
    corr_combined_df, 
    annot=True,          # Show the correlation values in the cells
    fmt=".2f",           # Format to 2 decimal places
    cmap=cmap,           # Color palette
    vmin=-1, vmax=1,     # Anchor the color scale from -1 to 1
    center=0,            # Center the colormap at 0 (white)
    square=True,         # Force cells to be square
    linewidths=.5,       # Add gridlines between cells
    cbar_kws={"shrink": .8, "label": "Correlation Coefficient"},
    xticklabels=clean_labels, 
    yticklabels=clean_labels,
    ax=ax                # Tell seaborn to plot on our customized black axes
)

# Rotate the x-axis labels so they don't overlap
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)

# Multi-line title to explicitly define the split matrix for reviewers
plt.title('Split Correlation Matrix of Exact Sampling Metrics\n'
          '(Lower Left: Pearson Linear $r$  |  Upper Right: Spearman Rank $\\rho$)', 
          fontsize=16, fontweight='bold', pad=20)

# Final formatting and export
plt.tight_layout()
filename = "Split_Metric_Correlation_Matrix_110m.png"
filepath = os.path.join(output_dir, filename)

plt.savefig(filepath, dpi=300, bbox_inches='tight')
plt.close()

print(f"Success! Split correlation matrix saved to: {filepath}")
