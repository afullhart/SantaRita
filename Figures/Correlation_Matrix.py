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

csv_path = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\SRER_NRI_Plots.csv"
df = pd.read_csv(csv_path)

# ====================================================================
# VARIABLE SELECTION & LABELING
# ====================================================================
# The exact metric columns from your dataset
target_cols = [
  'Exact_BGR_Pct', 
  'Exact_Fetch_m', 
  'Exact_LPI_Pct',
  'Exact_Gap_0_24', 
  'Exact_Gap_25_50', 
  'Exact_Gap_51_100',
  'Exact_Gap_101_200', 
  'Exact_Gap_gt_200', 
  'Exact_Herb_Woody_Ratio'
]

# Clean labels for the plot axes
clean_labels = [
  'Bare Ground %', 
  'Mean Fetch (m)', 
  'LPI %',
  'Gap 0-24 cm', 
  'Gap 25-50 cm', 
  'Gap 51-100 cm',
  'Gap 101-200 cm', 
  'Gap > 200 cm', 
  'Herb:Woody Ratio'
]

# ====================================================================
# CORRELATION CALCULATION & PLOTTING
# ====================================================================
print("Extracting variables and calculating Pearson correlation matrix...")

# Subset the dataframe to only the target columns
df_subset = df[target_cols]

# Calculate the correlation matrix
corr_matrix = df_subset.corr(method='pearson')

# Create a boolean mask to hide the main diagonal
diagonal_mask = np.eye(len(corr_matrix), dtype=bool)

# Set up the matplotlib figure and explicitly define the axes (ax)
fig, ax = plt.subplots(figsize=(12, 10))

# Paint the background of the axes black. 
# Because the diagonal is masked, the black background will show through.
ax.set_facecolor('black')

# Generate a custom diverging colormap (blue to red)
cmap = sns.color_palette("vlag", as_cmap=True)

# Draw the heatmap using seaborn
sns.heatmap(
  corr_matrix, 
  mask=diagonal_mask,  # Apply the mask to hide the diagonal
  annot=True,          # Show the correlation values in the cells
  fmt=".2f",           # Format to 2 decimal places
  cmap=cmap,           # Color palette
  vmin=-1, vmax=1,     # Anchor the color scale from -1 to 1
  center=0,            # Center the colormap at 0 (white)
  square=True,         # Force cells to be square
  linewidths=.5,       # Add gridlines between cells
  cbar_kws={"shrink": .8, "label": "Pearson Correlation Coefficient"},
  xticklabels=clean_labels, 
  yticklabels=clean_labels,
  ax=ax                # Tell seaborn to plot on our customized black axes
)

# Rotate the x-axis labels so they don't overlap
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)

plt.title('Correlation Matrix of Exact Sampling Metrics', fontsize=16, fontweight='bold', pad=20)

# Final formatting and export
plt.tight_layout()
filename = "Metric_Correlation_Matrix.png"
filepath = os.path.join(output_dir, filename)

plt.savefig(filepath, dpi=300, bbox_inches='tight')
plt.close()

print(f"Success! Correlation matrix saved to: {filepath}")
