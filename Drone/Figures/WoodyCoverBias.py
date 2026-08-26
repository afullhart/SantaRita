import os
import pandas as pd
import numpy as np
from scipy.stats import kendalltau, theilslopes
import matplotlib.pyplot as plt

# --- Set publication-quality plot parameters ---
plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
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
df['WP_Error'] = sampled_wp - exact_wp

# ==========================================
# --- CALCULATIONS FOR LEFT PANEL (FIT) ---
# ==========================================
print("Calculating non-parametric regression and 2000 bootstrap iterations. Please wait...")
np.random.seed(42) # For reproducible bootstrapping

# Main Theil-Sen fit
res_left = theilslopes(sampled_wp, exact_wp, 0.95)
slope, intercept, lo_slope, up_slope = res_left

# Bootstrapping for CI and PI bounds
n_boot = 2000
n_points = len(exact_wp)
x_grid = np.linspace(0, 60, 100) # Grid matches the 0-60 x-axis

boot_ci = np.zeros((n_boot, len(x_grid)))
boot_pi = np.zeros((n_boot, len(x_grid)))

main_preds = slope * exact_wp + intercept
residuals = sampled_wp - main_preds

for i in range(n_boot):
    # Resample pairs for CI
    idx = np.random.choice(np.arange(n_points), size=n_points, replace=True)
    x_b = exact_wp.iloc[idx]
    y_b = sampled_wp.iloc[idx]
    
    # Fit resampled data
    res_b = theilslopes(y_b, x_b)
    slope_b, int_b = res_b[0], res_b[1]
    
    # Generate line for CI
    boot_ci[i, :] = slope_b * x_grid + int_b
    
    # Resample residuals for PI
    res_idx = np.random.choice(np.arange(n_points), size=len(x_grid), replace=True)
    boot_pi[i, :] = boot_ci[i, :] + residuals.iloc[res_idx].values

# Extract 95% bounds
ci_lower = np.percentile(boot_ci, 2.5, axis=0)
ci_upper = np.percentile(boot_ci, 97.5, axis=0)
pi_lower = np.percentile(boot_pi, 2.5, axis=0)
pi_upper = np.percentile(boot_pi, 97.5, axis=0)

# Calculate PI Coverage and Widths
pi_lower_interp = np.interp(exact_wp, x_grid, pi_lower)
pi_upper_interp = np.interp(exact_wp, x_grid, pi_upper)
coverage = np.mean((sampled_wp >= pi_lower_interp) & (sampled_wp <= pi_upper_interp)) * 100

pi_widths = pi_upper - pi_lower
mean_pi_width = np.mean(pi_widths)

# Left panel specific stats
mae = np.mean(np.abs(sampled_wp - exact_wp))
mre = np.mean(np.abs(sampled_wp - exact_wp) / exact_wp) * 100
pbias = (np.sum(sampled_wp - exact_wp) / np.sum(exact_wp)) * 100

tau_left, p_val_left_two_sided = kendalltau(exact_wp, df['WP_Error'])
p_val_left = p_val_left_two_sided / 2  # Convert to one-tailed p-value

sig_diff = "Yes" if p_val_left < 0.05 else "No"
dir_diff = "< 1" if tau_left < 0 else "> 1"
sig_str = f"{sig_diff} ({dir_diff})" if sig_diff == "Yes" else "No"

stats_text_left = (
    f"Sen's Slope: {slope:.3f} ({lo_slope:.2f}, {up_slope:.2f})\n"
    f"Intercept (Med): {intercept:.3f}\n"
    f"MAE: {mae:.3f}\n"
    f"MRE%: {mre:.2f}\n"
    f"PBIAS: {pbias:.2f}%\n"
    f"--\n"
    f"MK Tau (y-x): {tau_left:.3f} (p: {p_val_left:.3f})\n"
    f"Sig diff from 1: {sig_str}\n"
    f"--\n"
    f"PI Coverage: {coverage:.1f}%\n"
    f"Mean PI Width: {mean_pi_width:.3f}"
)


# ================================================
# --- CALCULATIONS FOR RIGHT PANEL (RESIDUALS) ---
# ================================================
res_error = theilslopes(df['WP_Error'], exact_wp, 0.95)
error_slope = res_error[0]

tau_right, p_val_tau_right = kendalltau(exact_wp, df['WP_Error'])
p_val_one_tailed = p_val_tau_right / 2

stats_text_right = (
    f"Sen's Slope: {error_slope:.3f}\n"
    f"Kendall's $\\tau$: {tau_right:.3f}\n"
    f"$p$-value: {p_val_one_tailed:.3f}"
)


# ==========================
# --- PLOTTING ROUTINE ---
# ==========================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# ---------------------------------
# 1. LEFT PANEL: Regression & Bounds
# ---------------------------------
ax1.set_title(r'$\bf{a.}$ Non-Parametric Regression', fontsize=20, pad=15)

# Bootstrapped CI
ax1.fill_between(x_grid, ci_lower, ci_upper, color='red', alpha=0.3, zorder=1, label='95% Bootstrapped CI')

# Bootstrapped PI
ax1.plot(x_grid, pi_lower, color='orange', linestyle='--', linewidth=1.5, zorder=2, label='95% Bootstrapped PI')
ax1.plot(x_grid, pi_upper, color='orange', linestyle='--', linewidth=1.5, zorder=2)

# Scatter 
ax1.scatter(exact_wp, sampled_wp, color='#4C72B0', alpha=0.8, s=60, 
            edgecolor='black', linewidth=0.8, zorder=3, label='Data points')

# Regression Line
main_line = slope * x_grid + intercept
ax1.plot(x_grid, main_line, color='red', linewidth=2.0, zorder=4, label='Theil-Sen Regression Line')

# Aesthetics
ax1.set_xlabel('True Woody Cover (%)', fontweight='bold')
ax1.set_ylabel('Estimated Woody Cover (%)', fontweight='bold')
ax1.set_xlim(0, 60)
ax1.set_ylim(-5, 75)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.grid(axis='both', alpha=0.2, linestyle='-')
ax1.legend(loc='lower right', frameon=True, fontsize=14, framealpha=0.3)

# Stats text box (Inside top left)
props = dict(boxstyle='round', facecolor='white', alpha=0.3, edgecolor='gray')
ax1.text(0.04, 0.96, stats_text_left, transform=ax1.transAxes, fontsize=13,
         verticalalignment='top', horizontalalignment='left', bbox=props)


# ---------------------------------
# 2. RIGHT PANEL: Residuals
# ---------------------------------
ax2.set_title(r'$\bf{b.}$ Estimation Residuals', fontsize=20, pad=15)

ax2.axhline(0, color='black', linestyle='--', linewidth=1.5, zorder=1)

ax2.scatter(df['Exact_Woody_Pct'], df['WP_Error'], 
            color='#4C72B0', alpha=0.8, s=60, 
            edgecolor='black', linewidth=0.8, zorder=2)

# Aesthetics
ax2.set_xlabel('True Woody Cover (%)', fontweight='bold')
ax2.set_ylabel('Estimation Error (Percentage Points)', fontweight='bold')
ax2.set_xlim(0, 60)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.grid(axis='both', alpha=0.2, linestyle='-')

# Stats text box (Inside bottom right)
ax2.text(0.96, 0.04, stats_text_right, transform=ax2.transAxes, fontsize=13,
         verticalalignment='bottom', horizontalalignment='right', bbox=props)

# Layout adjustment and save
plt.tight_layout()
plt.savefig(output_filepath, format='svg', bbox_inches='tight')
plt.close()

print(f"\n1x2 subplot successfully saved to: {output_filepath}")
