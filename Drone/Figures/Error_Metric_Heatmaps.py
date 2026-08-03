import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# Suppress runtime warnings for division by zero (handled intentionally with NaNs)
warnings.filterwarnings('ignore', category=RuntimeWarning)

# ====================================================================
# GLOBAL PLOT SETTINGS
# ====================================================================
# Base font size for colorbars and ticks to ensure readability
plt.rcParams.update({'font.size': 14}) 

# ====================================================================
# 1. LOAD DATA & DEFINE MAPPINGS
# ====================================================================
df = pd.read_csv(r"C:\Users\andre\ScatterPlots\SRER_NRI_Plots_110m.csv")

# Metrics dictionary with abbreviated names, spaced CGF labels, and no ratio
metrics_info = {
    'BGR': {'exact': 'Exact_BGR_Pct', 'nri_pattern': 'NRI_BGR_{scale}_Pct'},
    'HP': {'exact': 'Exact_Herb_Pct', 'nri_pattern': 'NRI_Herb_{scale}_Pct'},
    'WP': {'exact': 'Exact_Woody_Pct', 'nri_pattern': 'NRI_Woody_{scale}_Pct'},
    'MF': {'exact': 'Exact_Fetch_m', 'nri_pattern': 'NRI_Fetch_{scale}'},
    'CGF 0-24 cm': {'exact': 'Exact_Gap_0_24', 'nri_pattern': 'NRI_Gap_0_24', 'is_gap': True},
    'CGF 25-50 cm': {'exact': 'Exact_Gap_25_50', 'nri_pattern': 'NRI_Gap_25_50', 'is_gap': True},
    'CGF 51-100 cm': {'exact': 'Exact_Gap_51_100', 'nri_pattern': 'NRI_Gap_51_100', 'is_gap': True},
    'CGF 101-200 cm': {'exact': 'Exact_Gap_101_200', 'nri_pattern': 'NRI_Gap_101_200', 'is_gap': True},
    'CGF +200 cm': {'exact': 'Exact_Gap_gt_200', 'nri_pattern': 'NRI_Gap_gt_200', 'is_gap': True},
}

scales = ['0cm', '25cm', '50cm', '100cm', '200cm']
results = []

# ====================================================================
# 2. NON-PARAMETRIC BOOTSTRAP FUNCTION
# ====================================================================
def compute_pi_width(X, Y, n_boot=2000, alpha=0.05):
    """Calculates the mean width of non-parametric prediction intervals."""
    n = len(X)
    pi_preds = np.zeros((n_boot, n))

    for i in range(n_boot):
        # Pairs Bootstrapping
        indices = np.random.choice(n, size=n, replace=True)
        x_boot = X[indices]
        y_boot = Y[indices]

        # Fit Bootstrap Model using Theil-Sen Estimator
        if np.std(x_boot) == 0:
            slope, intercept = 0, np.median(y_boot)
        else:
            res = stats.mstats.theilslopes(y_boot, x_boot)
            slope, intercept = res[0], res[1]

        # Calculate Empirical Residuals
        residuals = y_boot - (slope * x_boot + intercept)

        # Predict over the ORIGINAL X values
        y_hat = slope * X + intercept

        # Add Resampled Residuals
        noise = np.random.choice(residuals, size=n, replace=True)
        pi_preds[i, :] = y_hat + noise

    # Extract Percentiles for the Bounds
    lower_pct = (alpha / 2) * 100
    upper_pct = (1 - alpha / 2) * 100

    pi_lo = np.percentile(pi_preds, lower_pct, axis=0)
    pi_up = np.percentile(pi_preds, upper_pct, axis=0)

    # Return the Mean PI Width
    return np.mean(pi_up - pi_lo)


# ====================================================================
# 3. CALCULATE METRICS IN A SINGLE PASS
# ====================================================================
def calculate_stats(scl, x_vals, y_vals):
    """Calculates Pearson r, Spearman rho, Sen's Slope (with Sig), PBIAS, MAE, MRE, and MPIW"""
    if len(x_vals) > 1:
        # Correlations
        r, _ = stats.pearsonr(x_vals, y_vals)
        rho, _ = stats.spearmanr(x_vals, y_vals)
        
        # Regression & Errors
        sens, _, _, _ = stats.theilslopes(y_vals, x_vals)
        mae = np.mean(np.abs(y_vals - x_vals))
        mre = np.nanmean(np.abs(y_vals - x_vals) / np.where(x_vals == 0, np.nan, np.abs(x_vals))) * 100
        
        # Percent Bias
        sum_x = np.sum(x_vals)
        pbias = (np.sum(y_vals - x_vals) / sum_x) * 100 if sum_x != 0 else np.nan
        
        # Mean Prediction Interval Width (Bootstrapped)
        mpiw = compute_pi_width(x_vals.values, y_vals.values, n_boot=2000)
        
        # Significance of Sen's Slope != 1 (1-Tailed Mann-Kendall)
        y_trans = y_vals - x_vals
        tau, p_val_2tail = stats.kendalltau(x_vals, y_trans)
        if not pd.isna(tau):
            if sens > 1:
                p_1tail = p_val_2tail / 2 if tau > 0 else 1.0
            elif sens < 1:
                p_1tail = p_val_2tail / 2 if tau < 0 else 1.0
            else:
                p_1tail = 1.0
            is_sig = p_1tail < 0.05
        else:
            is_sig = False
            
        return {
            'Metric': metric, 'Scale': scl, 
            'r': r, 'rho': rho, 
            'Sens_Slope': sens, 'Sig_Slope': is_sig, 
            'PBIAS': pbias, 'MAE': mae, 'MRE': mre, 'MPIW': mpiw
        }
    return None

print("Calculating metrics and non-parametric prediction intervals. This may take a moment...")

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
# 4. PIVOT & REORDER DATA
# ====================================================================
r_pivot = res_df.pivot(index='Metric', columns='Scale', values='r')
rho_pivot = res_df.pivot(index='Metric', columns='Scale', values='rho')
sens_pivot = res_df.pivot(index='Metric', columns='Scale', values='Sens_Slope')
sig_pivot = res_df.pivot(index='Metric', columns='Scale', values='Sig_Slope')
pbias_pivot = res_df.pivot(index='Metric', columns='Scale', values='PBIAS')
mae_pivot = res_df.pivot(index='Metric', columns='Scale', values='MAE')
mre_pivot = res_df.pivot(index='Metric', columns='Scale', values='MRE')
mpiw_pivot = res_df.pivot(index='Metric', columns='Scale', values='MPIW')

# Applied ordering from top-to-bottom
metric_order = ['BGR', 'HP', 'WP', 'MF', 'CGF 0-24 cm', 'CGF 25-50 cm', 'CGF 51-100 cm', 'CGF 101-200 cm', 'CGF +200 cm']
scale_order = ['0cm', '25cm', '50cm', '100cm', '200cm']

# Safe reindexing to maintain desired visual plotting order
pivots = [r_pivot, rho_pivot, sens_pivot, sig_pivot, pbias_pivot, mae_pivot, mre_pivot, mpiw_pivot]
for i in range(len(pivots)):
    pivots[i] = pivots[i].reindex(index=[m for m in metric_order if m in pivots[i].index], columns=scale_order)
r_pivot, rho_pivot, sens_pivot, sig_pivot, pbias_pivot, mae_pivot, mre_pivot, mpiw_pivot = pivots

# ====================================================================
# 5. HELPERS FOR HEATMAP VISUALS
# ====================================================================
def add_strikes(ax, data_pivot):
    """Draws a diagonal strike through cells that contain NaN values."""
    for i in range(data_pivot.shape[0]):
        for j in range(data_pivot.shape[1]):
            if pd.isna(data_pivot.iloc[i, j]):
                # Draw a diagonal line from bottom-left to top-right of the cell
                ax.plot([j, j+1], [i+1, i], color='gray', lw=2.0)

# Create custom string annotations for Sen's Slope to append asterisks
sens_annot = np.empty_like(sens_pivot.values, dtype=object)
for i in range(sens_pivot.shape[0]):
    for j in range(sens_pivot.shape[1]):
        val = sens_pivot.iloc[i, j]
        is_sig = sig_pivot.iloc[i, j]
        if pd.isna(val):
            sens_annot[i, j] = ""
        else:
            sens_annot[i, j] = f"{val:.2f}*" if is_sig else f"{val:.2f}"

# ====================================================================
# 6. PLOTTING - FIGURE 1: Pearson r, Spearman rho, & Sen's Slope
# ====================================================================
# Widened to 21x8 to prevent horizontal text overlap while keeping fonts large
fig1, axes1 = plt.subplots(1, 3, figsize=(21, 8))

# Left Plot: Pearson r
sns.heatmap(r_pivot, annot=True, fmt=".2f", cmap="Blues", ax=axes1[0], 
            annot_kws={"size": 15}, cbar_kws={'label': 'Pearson $r$'}, mask=r_pivot.isnull()) 
add_strikes(axes1[0], r_pivot)  
axes1[0].set_title('Pearson Correlation ($r$)\nExact vs Sampled', pad=20, fontsize=18, fontweight='bold')
axes1[0].set_ylabel('Ground Cover Metric', fontsize=16, fontweight='bold')
axes1[0].set_xlabel('NRI Transect Sampling Scale', fontsize=16, fontweight='bold')
axes1[0].tick_params(axis='both', which='major', labelsize=14)

# Middle Plot: Spearman rho
sns.heatmap(rho_pivot, annot=True, fmt=".2f", cmap="Blues", ax=axes1[1], 
            annot_kws={"size": 15}, cbar_kws={'label': r'Spearman $\rho$'}, mask=rho_pivot.isnull()) 
add_strikes(axes1[1], rho_pivot)  
axes1[1].set_title(r'Spearman Rank Correlation ($\rho$)' + '\nExact vs Sampled', pad=20, fontsize=18, fontweight='bold')
axes1[1].set_ylabel('')
axes1[1].set_xlabel('NRI Transect Sampling Scale', fontsize=16, fontweight='bold')
axes1[1].tick_params(axis='both', which='major', labelsize=14)

# Right Plot: Sen's Slope (with Significance Asterisks)
sns.heatmap(sens_pivot, annot=sens_annot, fmt="", cmap="vlag", center=1.0, ax=axes1[2], 
            annot_kws={"size": 15}, cbar_kws={'label': "Sen's Slope"}, mask=sens_pivot.isnull()) 
add_strikes(axes1[2], sens_pivot)  
axes1[2].set_title("Sen's Slope\nExact vs Sampled (* Sig < 1 or > 1)", pad=20, fontsize=18, fontweight='bold')
axes1[2].set_ylabel('')
axes1[2].set_xlabel('NRI Transect Sampling Scale', fontsize=16, fontweight='bold')
axes1[2].tick_params(axis='both', which='major', labelsize=14)

# Use w_pad to ensure the labels between subplots don't overlap
plt.tight_layout(w_pad=2.0)
plt.savefig(r'C:\Users\andre\ScatterPlots\Exact_vs_Sampled_Heatmaps.svg', format='svg', dpi=300, bbox_inches='tight')
plt.show()

# ====================================================================
# 7. PLOTTING - FIGURE 2: PBIAS, MAE, & MRE
# ====================================================================
# Widened to 21x8 to prevent horizontal text overlap
fig2, axes2 = plt.subplots(1, 3, figsize=(21, 8))

# Left Plot: Percent Bias (PBIAS)
sns.heatmap(pbias_pivot, annot=True, fmt=".2f", cmap="vlag", center=0, ax=axes2[0], 
            annot_kws={"size": 15}, cbar_kws={'label': 'Percent Bias (%)'}, mask=pbias_pivot.isnull()) 
add_strikes(axes2[0], pbias_pivot)  
axes2[0].set_title('Percent Bias (PBIAS)\nExact vs Sampled (%)', pad=20, fontsize=18, fontweight='bold')
axes2[0].set_ylabel('Ground Cover Metric', fontsize=16, fontweight='bold')
axes2[0].set_xlabel('NRI Transect Sampling Scale', fontsize=16, fontweight='bold')
axes2[0].tick_params(axis='both', which='major', labelsize=14)

# Middle Plot: Mean Absolute Error (MAE)
sns.heatmap(mae_pivot, annot=True, fmt=".2f", cmap="Reds", ax=axes2[1], 
            annot_kws={"size": 15}, cbar_kws={'label': 'Mean Absolute Error'}, mask=mae_pivot.isnull()) 
add_strikes(axes2[1], mae_pivot)  
axes2[1].set_title('Mean Absolute Error (MAE)\nExact vs Sampled', pad=20, fontsize=18, fontweight='bold')
axes2[1].set_ylabel('')
axes2[1].set_xlabel('NRI Transect Sampling Scale', fontsize=16, fontweight='bold')
axes2[1].tick_params(axis='both', which='major', labelsize=14)

# Right Plot: Mean Relative Error (MRE)
sns.heatmap(mre_pivot, annot=True, fmt=".2f", cmap="Oranges", ax=axes2[2], 
            annot_kws={"size": 15}, cbar_kws={'label': "Mean Relative Error (%)"}, mask=mre_pivot.isnull()) 
add_strikes(axes2[2], mre_pivot)  
axes2[2].set_title("Mean Relative Error (MRE)\nExact vs Sampled (% of Exact)", pad=20, fontsize=18, fontweight='bold')
axes2[2].set_ylabel('')
axes2[2].set_xlabel('NRI Transect Sampling Scale', fontsize=16, fontweight='bold')
axes2[2].tick_params(axis='both', which='major', labelsize=14)

plt.tight_layout(w_pad=2.0)
plt.savefig(r'C:\Users\andre\ScatterPlots\Error_Metrics_Heatmaps.svg', format='svg', dpi=300, bbox_inches='tight')
plt.show()

# ====================================================================
# 8. PLOTTING - FIGURE 3: Mean Prediction Interval Width (MPIW)
# ====================================================================
# Adjusted to match the individual subplot proportions of the widened figures
fig3, ax3 = plt.subplots(figsize=(8, 8))

# Single Plot: Mean Prediction Interval Width
sns.heatmap(mpiw_pivot, annot=True, fmt=".2f", cmap="Purples", ax=ax3, 
            annot_kws={"size": 15}, cbar_kws={'label': 'Mean PI Width'}, mask=mpiw_pivot.isnull()) 
add_strikes(ax3, mpiw_pivot)  
ax3.set_title("Mean 95% Prediction Interval Width\nExact vs Sampled", pad=20, fontsize=18, fontweight='bold')
ax3.set_ylabel('Ground Cover Metric', fontsize=16, fontweight='bold')
ax3.set_xlabel('NRI Transect Sampling Scale', fontsize=16, fontweight='bold')
ax3.tick_params(axis='both', which='major', labelsize=14)

plt.tight_layout()
plt.savefig(r'C:\Users\andre\ScatterPlots\Mean_PI_Width_Heatmap.svg', format='svg', dpi=300, bbox_inches='tight')
plt.show()

print("Processing complete. Heatmaps saved to C:\\Users\\andre\\ScatterPlots\\")
