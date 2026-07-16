import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats

# ====================================================================
# FOLDER SETUP & ROUTING
# ====================================================================
# Define base output directory
base_dir = r'C:\Users\andre\ScatterPlots'

# Define and create the new subfolders
exact_dir = os.path.join(base_dir, 'Exact_vs_Exact')
sampled_dir = os.path.join(base_dir, 'Exact_vs_Sampled')

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

    ('Exact_Fetch_m', 'NRI_Fetch_0cm'),
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
    
    # Retrieve predictions along with confidence and prediction intervals (alpha=0.05 -> 95%)
    pred = model.get_prediction(X_sm)
    pred_df = pred.summary_frame(alpha=0.05)
    Y_pred = pred_df['mean'].values

    # ==================================================================
    # METRICS CALCULATION
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

        # Base text string with MK Test and PBIAS included
        text_str = (f'OLS Slope: {slope:.3f}\nIntercept: {intercept:.3f}\n$R^2$: {r_squared:.3f}\n'
                    f'RMSE: {rmse:.3f}\nMAE: {mae:.3f}\nMRE%: {mre_pct:.2f}\n'
                    f'PBIAS: {pbias:.2f}%\n'
                    f'---\n'
                    f"Sen's Slope: {sens_slope:.3f} ({sens_lo:.2f}, {sens_up:.2f})\n"
                    f"MK Tau (y-x): {tau:.3f} (p: {p_value_one_tailed:.3f})\n"
                    f"Sig diff from 1: {sig_text}")

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
                    f"Sen's Slope: {sens_slope:.3f} ({sens_lo:.2f}, {sens_up:.2f})")

        mk_tau_out = np.nan
        mk_pval_out = np.nan
        test_dir_out = None
        is_sig_out = None
        target_dir = exact_dir

    # Package metrics for the CSVs
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
        'Slope_Sig_Diff_From_1': is_sig_out
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
    
    plt.fill_between(X, pred_df['mean_ci_lower'], pred_df['mean_ci_upper'], 
                     color='red', alpha=0.3, label='95% Confidence Interval')

    plt.plot(X, pred_df['obs_ci_lower'], color='orange', linestyle='--', label='95% Prediction Interval')
    plt.plot(X, pred_df['obs_ci_upper'], color='orange', linestyle='--')
    
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
