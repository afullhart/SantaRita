import os
import pandas as pd
from scipy.stats import pearsonr, kendalltau, theilslopes
import matplotlib.pyplot as plt

# --- Set publication-quality plot parameters ---
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'axes.linewidth': 1.2,
    'font.family': 'sans-serif'
})

# 1. Define the file path
nri_file_path = r'C:\Users\andre\ScatterPlots\SRER_NRI_Plots_110m.csv'

# Get the directory of the NRI file
nri_directory = os.path.dirname(os.path.abspath(nri_file_path))
output_filepath = os.path.join(nri_directory, 'Full_Cover_Estimation_Error.svg')

# 2. Load the dataset
df = pd.read_csv(nri_file_path)

# 3. Define the exact values and the sampled transect estimates
exact_wp = df['Exact_Woody_Pct']
sampled_wp = df['NRI_Woody_0cm_Pct']

# 4. Calculate the raw estimation error for each plot
df['WP_Error'] = sampled_wp - exact_wp

# --- Calculate statistics for annotations ---
# 1. Pearson R and R-squared
r_val, p_val_r = pearsonr(exact_wp, df['WP_Error'])
r_squared = r_val**2

# 2. Theil-Sen Slope of Residuals
res_error = theilslopes(df['WP_Error'], exact_wp, 0.95)
error_slope = res_error[0]

# 3. Kendall's Tau (Significance)
tau, p_val_tau = kendalltau(exact_wp, df['WP_Error'])
p_val_one_tailed = p_val_tau / 2

# Create bins from 0 to 60 in 5% increments
bins = list(range(0, 65, 5))

# 6. Visualize the error across the full range
fig, ax = plt.subplots(figsize=(8, 6))

for b in bins:
    ax.axvline(b, color='gray', linestyle=':', alpha=0.4, zorder=0)

ax.axhline(0, color='black', linestyle='--', linewidth=1.5, zorder=1)

scatter = ax.scatter(df['Exact_Woody_Pct'], df['WP_Error'], 
                     color='#4C72B0', alpha=0.8, s=60, 
                     edgecolor='black', linewidth=0.8, zorder=2)

ax.set_xlabel('True Woody Cover (%)', fontweight='bold')
ax.set_ylabel('Estimation Error (Percentage Points)', fontweight='bold')
ax.set_xlim(0, 60)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.2, linestyle='-')

# --- Add the statistical annotation text box ---
stats_text = (
    f"$R^2$: {r_squared:.3f}\n"
    f"Sen's Slope: {error_slope:.3f}\n"
    f"Kendall's $\\tau$: {tau:.3f}\n"
    f"$p$-value: {p_val_one_tailed:.3f}"
)

# --- NEW: Place the text box outside the plot to the right ---
props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray')
ax.text(1.02, 0.95, stats_text, transform=ax.transAxes, fontsize=11,
        verticalalignment='top', horizontalalignment='left', bbox=props)

plt.tight_layout()

# Save the figure (bbox_inches='tight' ensures the outside legend isn't clipped)
plt.savefig(output_filepath, format='svg', bbox_inches='tight')
plt.close()

print(f"\nScatter plot with stats successfully saved to: {output_filepath}")
