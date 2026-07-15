import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats

# ====================================================================
# NON-PARAMETRIC BOOTSTRAPPING FUNCTION
# ====================================================================
def compute_bootstrap_intervals(X, Y, n_boot=10000, alpha=0.05):
    """
    Calculates non-parametric confidence and prediction intervals
    using pairs bootstrapping and empirical residual resampling.
    """
    n = len(X)
    pi_preds = np.zeros((n_boot, n))
    ci_preds = np.zeros((n_boot, n))

    for i in range(n_boot):
        # 1. Pairs Bootstrapping: Resample indices with replacement
        indices = np.random.choice(n, size=n, replace=True)
        x_boot = X[indices]
        y_boot = Y[indices]

        # 2. Fit Bootstrap Model (using np.polyfit for speed inside the loop)
        if np.std(x_boot) == 0:
            slope, intercept = 0, np.mean(y_boot)
        else:
            slope, intercept = np.polyfit(x_boot, y_boot, 1)

        # 3. Calculate Empirical Residuals of this specific Bootstrap Model
        residuals = y_boot - (slope * x_boot + intercept)

        # 4. Predict over the ORIGINAL X values (to anchor our plot lines)
        y_hat = slope * X + intercept
        ci_preds[i, :] = y_hat

        # 5. Add Resampled Residuals to simulate Prediction Interval scatter
        noise = np.random.choice(residuals, size=n, replace=True)
        pi_preds[i, :] = y_hat + noise

    # 6. Extract Percentiles for the Bounds (Default 95% = 2.5th and 97.5th percentiles)
    lower_pct = (alpha / 2) * 100
    upper_pct = (1 - alpha / 2) * 100

    ci_lo = np.percentile(ci_preds, lower_pct, axis=0)
    ci_up = np.percentile(ci_preds, upper_pct, axis=0)
    pi_lo = np.percentile(pi_preds, lower_pct, axis=0)
    pi_up = np.percentile(pi_preds, upper_pct, axis=0)

    return ci_lo, ci_up, pi_lo, pi_up


# ====================================================================
# FOLDER SETUP & ROUTING
# ====================================================================
# Define base output directory
base_dir = r'C:\Users\andre\ScatterPlots'

# Define and create the new subfolders
exact_dir = os.path.join(base_dir, 'Exact_vs_Exact_NonPara')
sampled_dir = os.path.join(base_dir, 'Exact_vs_Sampled_NonPara')

os.makedirs(exact_dir, exist_ok=True)
os.makedirs(sampled_dir, exist_ok=True)

# Read the uploaded CSV
df = pd.read_csv(r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\SRER_NRI_Plots_110m.csv')

# ====================================================================
# PAIR DEFINITIONS
# ====================================================================
pairs = [
    # --- Exact vs Sampled (NRI) Pairs ---
    ('Exact_BGR_Pct', 'NRI_BGR_0cm_Pct'),
    ('Exact_BGR_Pct', 'NRI_BGR_25cm_Pct'),
    ('Exact_BGR_Pct', 'NRI_BGR_50cm_Pct'),
    ('Exact_BGR_Pct', 'NRI_BGR_100cm_Pct'),
    ('Exact_BGR_Pct', 'NRI_BGR_200cm_Pct'),
    
    ('Exact_Herb_Pct', 'NRI_Herb_0cm_Pct'),
    ('Exact_Herb_Pct', 'NRI_Herb_25cm_Pct'),
    ('Exact_Herb_Pct', 'NRI_Herb_50cm_Pct'),
    ('Exact_Herb_Pct', 'NRI_Herb_100cm_Pct'),
    ('Exact_Herb_Pct', 'NRI_Herb_200cm_Pct'),
    
    ('Exact_Woody_Pct', 'NRI_Woody_0cm_Pct'),
    ('Exact_Woody_Pct', 'NRI_Woody_25cm_Pct'),
    ('Exact_Woody_Pct', 'NRI_Woody_50cm_Pct'),
    ('Exact_Woody_Pct', 'NRI_Woody_100cm_Pct'),
    ('Exact_Woody_Pct', 'NRI_Woody_200cm_Pct'),

    ('Exact_Fetch_m', 'NRI_Fetch_25cm'),
    ('Exact_Fetch_m', 'NRI_Fetch_50cm'),
    ('Exact_Fetch_m', 'NRI_Fetch_100cm'),
    ('Exact_Fetch_m', 'NRI_Fetch_200cm'),
    
    ('Exact_Herb_Woody_Ratio', 'NRI_HW_Ratio_0cm'),
    ('Exact_Herb_Woody_Ratio', 'NRI_HW_Ratio_25cm'),
    ('Exact_Herb_Woody_Ratio', 'NRI_HW_Ratio_50cm'),
    ('Exact_Herb_Woody_Ratio', 'NRI_HW_Ratio_100cm'),
    ('Exact_Herb_Woody_Ratio', 'NRI_HW_Ratio_200cm'),
    
    ('Exact_Gap_0_24', 'NRI_Gap_0_24'),
    ('Exact_Gap_25_50', 'NRI_Gap_25_50'),
    ('Exact_Gap_51_100', 'NRI_Gap_51_100'),
    ('Exact_Gap_101_200', 'NRI_Gap_101_200'),
    ('Exact_Gap_gt_200', 'NRI_Gap_gt_200'),
    
    # --- Exact vs Exact Pairs ---
    ('Exact_Gap_0_24', 'Exact_LPI_Pct'),
    ('Exact_Gap_25_50', 'Exact_LPI_Pct'),
    ('Exact_Gap_51_100', 'Exact_LPI_Pct'),
    ('Exact_Gap_101_200', 'Exact_LPI_Pct'),
    ('Exact_Gap_gt_200', 'Exact_LPI_Pct'),

    ('Exact_BGR_Pct', 'Exact_Fetch_m'),
    ('Exact_BGR_Pct', 'Exact_LPI_Pct'),
    ('Exact_BGR_Pct', 'Exact_Gap_0_24'),
    ('Exact_BGR_Pct', 'Exact_Gap_25_50'),
    ('Exact_BGR_Pct', 'Exact_Gap_51_100'),
    ('Exact_BGR_Pct', 'Exact_Gap_101_200'),
    ('Exact_BGR_Pct', 'Exact_Gap_gt_200'),
    ('Exact_BGR_Pct', 'Exact_Herb_Woody_Ratio'),

    ('Exact_Fetch_m', 'Exact_LPI_Pct'),
    ('Exact_Fetch_m', 'Exact_Gap_0_24'),
    ('Exact_Fetch_m', 'Exact_Gap_25_50'),
    ('Exact_Fetch_m', 'Exact_Gap_51_100'),
    ('Exact_Fetch_m', 'Exact_Gap_101_200'),
    ('Exact_Fetch_m', 'Exact_Gap_gt_200'),
    ('Exact_Fetch_m', 'Exact_Herb_Woody_Ratio'),

    ('Exact_LPI_Pct', 'Exact_Herb_Woody_Ratio'),

    ('Exact_Gap_0_24', 'Exact_Gap_25_50'),
    ('Exact_Gap_0_24', 'Exact_Gap_51_100'),
    ('Exact_Gap_0_24', 'Exact_Gap_101_200'),
    ('Exact_Gap_0_24', 'Exact_Gap_gt_200'),
    ('Exact_Gap_0_24', 'Exact_Herb_Woody_Ratio'),

    ('Exact_Gap_25_50', 'Exact_Gap_51_100'),
    ('Exact_Gap_25_50', 'Exact_Gap_101_200'),
    ('Exact_Gap_25_50', 'Exact_Gap_gt_200'),
    ('Exact_Gap_25_50', 'Exact_Herb_Woody_Ratio'),

    ('Exact_Gap_51_100', 'Exact_Gap_101_200'),
    ('Exact_Gap_51_100', 'Exact_Gap_gt_200'),
    ('Exact_Gap_51_100', 'Exact_Herb_Woody_Ratio'),

    ('Exact_Gap_101_200', 'Exact_Gap_gt_200'),
    ('Exact_Gap_101_200', 'Exact_Herb_Woody_Ratio'),

    ('Exact_Gap_gt_200', 'Exact_Herb_Woody_Ratio'),

    # --- Herb & Woody Percentages ---
    ('Exact_Herb_Pct', 'Exact_Woody_Pct'),
    
    ('Exact_BGR_Pct', 'Exact_Herb_Pct'),
    ('Exact_BGR_Pct', 'Exact_Woody_Pct'),
    
    ('Exact_Fetch_m', 'Exact_Herb_Pct'),
    ('Exact_Fetch_m', 'Exact_Woody_Pct'),
    
    ('Exact_LPI_Pct', 'Exact_Herb_Pct'),
    ('Exact_LPI_Pct', 'Exact_Woody_Pct'),
    
    ('Exact_Herb_Pct', 'Exact_Herb_Woody_Ratio'),
    ('Exact_Woody_Pct', 'Exact_Herb_Woody_Ratio'),
    
    ('Exact_Gap_0_24', 'Exact_Herb_Pct'),
    ('Exact_Gap_25_50', 'Exact_Herb_Pct'),
    ('Exact_Gap_51_100', 'Exact_Herb_Pct'),
    ('Exact_Gap_101_200', 'Exact_Herb_Pct'),
    ('Exact_Gap_gt_200', 'Exact_Herb_Pct'),
    
    ('Exact_Gap_0_24', 'Exact_Woody_Pct'),
    ('Exact_Gap_25_50', 'Exact_Woody_Pct'),
    ('Exact_Gap_51_100', 'Exact_Woody_Pct'),
    ('Exact_Gap_101_200', 'Exact_Woody_Pct'),
    ('Exact_Gap_gt_200', 'Exact_Woody_Pct')
]

# Lists to store metrics for the CSV outputs
sampled_metrics = []
exact_metrics = []

# ====================================================================
# PLOTTING LOOP
# ====================================================================
for i, (x_col, y_col) in enumerate(pairs):
    plt.figure(figsize=(8, 6))
    
    # Check if columns exist before plotting
    if x_col not in df.columns or y_col not in df.columns:
        print(f"Skipping plot: Columns {x_col} or {y_col} not found in the dataset.")
        plt.close()
        continue
    
    # Drop missing values and sort by X to ensure continuous fill boundaries
    subset = df[[x_col, y_col]].dropna().sort_values(by=x_col)
    
    # Skip plotting if dropping NaNs leaves us with no data
    if subset.empty:
        print(f"Skipping plot for {y_col} vs {x_col} due to insufficient valid data.")
        plt.close()
        continue

    X = subset[x_col].values
    Y = subset[y_col].values
    
    # Fit Ordinary Least Squares (OLS) model for the regression calculations
    X_sm = sm.add_constant(X)
    model = sm.OLS(Y, X_sm).fit()
    
    # Extract Slope, Intercept, and R-squared
    intercept = model.params[0]
    slope = model.params[1]
    r_squared = model.rsquared
    
    # Base Predictions from OLS
    Y_pred = model.predict(X_sm)

    # ------------------------------------------------------------------
    # Calculate Non-Parametric Bootstrapped Intervals
    # ------------------------------------------------------------------
    ci_lo, ci_up, pi_lo, pi_up = compute_bootstrap_intervals(X, Y, n_boot=10000)

    # ------------------------------------------------------------------
    # PREDICTION INTERVAL METRICS (MPIW & PICP)
    # ------------------------------------------------------------------
    pi_widths = pi_up - pi_lo
    
    min_pi_width = np.min(pi_widths)
    max_pi_width = np.max(pi_widths)
    mean_pi_width = np.mean(pi_widths)
    
    # Calculate Coverage (PICP): How many actual Y points fall inside the bounds?
    points_inside = np.sum((Y >= pi_lo) & (Y <= pi_up))
    picp = (points_inside / len(Y)) * 100

    # ==================================================================
    # CORE METRICS CALCULATION
    # ==================================================================
    rmse = np.sqrt(np.mean((Y - Y_pred)**2))
    mae = np.mean(np.abs(Y - Y_pred))
    
    # Calculate MRE%, ignoring true zeros to prevent division by zero errors
    with np.errstate(divide='ignore', invalid='ignore'):
        mre_array = np.abs(Y - Y_pred) / Y
        mre_array[np.isinf(mre_array)] = np.nan
        mre_pct = np.nanmean(mre_array) * 100

        # Calculate Percent Bias (PBIAS) between Evaluated (Y) and Reference (X)
        sum_X = np.sum(X)
        if sum_X != 0:
            pbias = (np.sum(Y - X) / sum_X) * 100
        else:
            pbias = np.nan

    # Calculate Sen's Slope for ALL pairs
    sens_slope, sens_intercept, sens_lo, sens_up = stats.mstats.theilslopes(Y, X, alpha=0.95)

    # Check which category this pair belongs to
    is_sampled = 'NRI_' in x_col or 'NRI_' in y_col

    if is_sampled:
        # ==================================================================
        # MANN-KENDALL (DYNAMIC 1-TAILED TEST) - SAMPLED ONLY
        # ==================================================================
        Y_transformed = Y - (1.0 * X)
        tau, p_value_two_tailed = stats.kendalltau(X, Y_transformed)
        
        # Determine direction based on Sen's Slope
        if not pd.isna(tau):
            if sens_slope > 1:
                test_direction = "> 1"
                p_value_one_tailed = p_value_two_tailed / 2 if tau > 0 else 1.0
            elif sens_slope < 1:
                test_direction = "< 1"
                p_value_one_tailed = p_value_two_tailed / 2 if tau < 0 else 1.0
            else:
                test_direction = "== 1"
                p_value_one_tailed = 1.0
            
            is_sig = p_value_one_tailed < 0.05
            sig_text = f"Yes ({test_direction})" if is_sig else "No"
        else:
            p_value_one_tailed = np.nan
            test_direction = "N/A"
            is_sig = False
            sig_text = "N/A"

        # Base text string with MK Test, PBIAS, and PI Metrics included
        text_str = (f'OLS Slope: {slope:.3f}\nIntercept: {intercept:.3f}\n$R^2$: {r_squared:.3f}\n'
                    f'RMSE: {rmse:.3f}\nMAE: {mae:.3f}\nMRE%: {mre_pct:.2f}\n'
                    f'PBIAS: {pbias:.2f}%\n'
                    f'---\n'
                    f"Sen's Slope: {sens_slope:.3f} ({sens_lo:.2f}, {sens_up:.2f})\n"
                    f"MK Tau (y-x): {tau:.3f} (p: {p_value_one_tailed:.3f})\n"
                    f"Sig diff from 1: {sig_text}\n"
                    f'---\n'
                    f'PI Coverage: {picp:.1f}%\n'
                    f'Mean PI Width: {mean_pi_width:.3f}\n'
                    f'Min/Max PI Width: {min_pi_width:.3f} / {max_pi_width:.3f}')

        # Variables to pass into the CSV dictionary
        mk_tau_out = tau
        mk_pval_out = p_value_one_tailed
        test_dir_out = test_direction
        is_sig_out = is_sig
        target_dir = sampled_dir

    else:
        # ==================================================================
        # EXACT VS EXACT - SKIP MANN KENDALL
        # ==================================================================
        text_str = (f'OLS Slope: {slope:.3f}\nIntercept: {intercept:.3f}\n$R^2$: {r_squared:.3f}\n'
                    f'RMSE: {rmse:.3f}\nMAE: {mae:.3f}\nMRE%: {mre_pct:.2f}\n'
                    f'PBIAS: {pbias:.2f}%\n'
                    f'---\n'
                    f"Sen's Slope: {sens_slope:.3f} ({sens_lo:.2f}, {sens_up:.2f})\n"
                    f'---\n'
                    f'PI Coverage: {picp:.1f}%\n'
                    f'Mean PI Width: {mean_pi_width:.3f}\n'
                    f'Min/Max PI Width: {min_pi_width:.3f} / {max_pi_width:.3f}')

        mk_tau_out = np.nan
        mk_pval_out = np.nan
        test_dir_out = None
        is_sig_out = None
        target_dir = exact_dir

    # Package metrics for the CSVs (Now including PI Metrics)
    current_metrics = {
        'X_Column': x_col,
        'Y_Column': y_col,
        'OLS_Slope': slope,
        'OLS_Intercept': intercept,
        'R2': r_squared,
        'RMSE': rmse,
        'MAE': mae,
        'MRE_Pct': mre_pct,
        'PBIAS_Pct': pbias,
        'Sens_Slope': sens_slope,
        'Sens_CI_Lower': sens_lo,
        'Sens_CI_Upper': sens_up,
        'MK_Tau_Transformed': mk_tau_out,
        'MK_p_value_1tailed': mk_pval_out,
        'MK_Test_Direction': test_dir_out,
        'Slope_Sig_Diff_From_1': is_sig_out,
        'PICP_Pct': picp,
        'Mean_PI_Width': mean_pi_width,
        'Min_PI_Width': min_pi_width,
        'Max_PI_Width': max_pi_width
    }

    if is_sampled:
        sampled_metrics.append(current_metrics)
    else:
        exact_metrics.append(current_metrics)

    # ==================================================================
    # GRAPHING
    # ==================================================================
    plt.scatter(X, Y, alpha=0.6, label='Data points', color='blue')
    plt.plot(X, Y_pred, color='red', label='Regression Line')
    
    # Using the new bootstrapped arrays for Confidence and Prediction Intervals
    plt.fill_between(X, ci_lo, ci_up, 
                     color='red', alpha=0.3, label='95% Bootstrapped CI')

    plt.plot(X, pi_lo, color='orange', linestyle='--', label='95% Bootstrapped PI')
    plt.plot(X, pi_up, color='orange', linestyle='--')
    
    plt.gca().text(0.05, 0.95, text_str, transform=plt.gca().transAxes, fontsize=9,
                   verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))

    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(f'{y_col} vs {x_col}')
    
    plt.legend(loc='lower right')
    plt.tight_layout()

    filename = f'{x_col}_vs_{y_col}.png'
    filepath = os.path.join(target_dir, filename)
    plt.savefig(filepath)
    plt.close()
    
    print(f'Saved: {filename} -> {os.path.basename(target_dir)}')

# ====================================================================
# EXPORT METRICS TO CSV
# ====================================================================
if sampled_metrics:
    sampled_df = pd.DataFrame(sampled_metrics)
    sampled_csv_path = os.path.join(base_dir, 'Exact_vs_Sampled_Metrics.csv')
    sampled_df.to_csv(sampled_csv_path, index=False)
    print(f'\nSuccessfully exported Exact_vs_Sampled metrics to: {sampled_csv_path}')

if exact_metrics:
    exact_df = pd.DataFrame(exact_metrics)
    exact_csv_path = os.path.join(base_dir, 'Exact_vs_Exact_Metrics.csv')
    exact_df.to_csv(exact_csv_path, index=False)
    print(f'Successfully exported Exact_vs_Exact metrics to: {exact_csv_path}')

print(f'\nAll plots and tables successfully generated inside: {base_dir}')
