import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

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
df = pd.read_csv(r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\SRER_NRI_Plots.csv')

# ====================================================================
# PAIR DEFINITIONS
# ====================================================================
pairs = [
  # --- Exact vs Sampled (NRI) Pairs ---
  ('Exact_BGR_Pct', 'NRI_BGR_0cm_Pct'),
  ('Exact_BGR_Pct', 'NRI_BGR_50cm_Pct'),
  ('Exact_BGR_Pct', 'NRI_BGR_100cm_Pct'),
  ('Exact_Fetch_m', 'NRI_Fetch_50cm'),
  ('Exact_Fetch_m', 'NRI_Fetch_100cm'),
  ('Exact_Herb_Woody_Ratio', 'NRI_HW_Ratio_0cm'),
  ('Exact_Herb_Woody_Ratio', 'NRI_HW_Ratio_50cm'),
  ('Exact_Herb_Woody_Ratio', 'NRI_HW_Ratio_100cm'),
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

  ('Exact_Gap_gt_200', 'Exact_Herb_Woody_Ratio')
]

# Lists to store metrics for the CSV outputs
sampled_metrics = []
exact_metrics = []

# ====================================================================
# PLOTTING LOOP
# ====================================================================
for i, (x_col, y_col) in enumerate(pairs):
  plt.figure(figsize=(8, 6))
  
  # Drop missing values and sort by X to ensure continuous fill boundaries
  # This perfectly handles our undefined Herb-to-Woody ratios (0 woody pixels)
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
  # METRICS CALCULATION (FOR ALL PLOTS)
  # ==================================================================
  rmse = np.sqrt(np.mean((Y - Y_pred)**2))
  mae = np.mean(np.abs(Y - Y_pred))
  
  # Calculate MRE%, ignoring true zeros to prevent division by zero errors
  with np.errstate(divide='ignore', invalid='ignore'):
    mre_array = np.abs(Y - Y_pred) / Y
    mre_array[np.isinf(mre_array)] = np.nan
    mre_pct = np.nanmean(mre_array) * 100

  # Package metrics for the CSVs
  current_metrics = {
    'X_Column': x_col,
    'Y_Column': y_col,
    'Slope': slope,
    'Intercept': intercept,
    'R2': r_squared,
    'RMSE': rmse,
    'MAE': mae,
    'MRE_Pct': mre_pct
  }

  # Extended text string applied to ALL plots
  text_str = (f'Slope: {slope:.3f}\nIntercept: {intercept:.3f}\n$R^2$: {r_squared:.3f}\n'
              f'RMSE: {rmse:.3f}\nMAE: {mae:.3f}\nMRE%: {mre_pct:.2f}')

  # ==================================================================
  # ROUTING LOGIC
  # ==================================================================
  is_sampled = 'NRI_' in x_col or 'NRI_' in y_col
  
  if is_sampled:
    target_dir = sampled_dir
    sampled_metrics.append(current_metrics)
  else:
    target_dir = exact_dir
    exact_metrics.append(current_metrics)

  # ==================================================================
  # GRAPHING
  # ==================================================================
  # Plot Datapoints
  plt.scatter(X, Y, alpha=0.6, label='Data points', color='blue')
  
  # Plot Regression Line
  plt.plot(X, Y_pred, color='red', label='Regression Line')
  
  # 95% Confidence Interval for the regression line
  plt.fill_between(X, pred_df['mean_ci_lower'], pred_df['mean_ci_upper'], 
    color='red', alpha=0.3, label='95% Confidence Interval')

  # 95% Prediction Interval for individual data points
  plt.plot(X, pred_df['obs_ci_lower'], color='orange', linestyle='--', label='95% Prediction Interval')
  plt.plot(X, pred_df['obs_ci_upper'], color='orange', linestyle='--')
  
  # Add the generated text box to the plot (top left)
  plt.gca().text(0.05, 0.95, text_str, transform=plt.gca().transAxes, fontsize=10,
    verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))

  # Formatting the plot
  plt.xlabel(x_col)
  plt.ylabel(y_col)
  plt.title(f'{y_col} vs {x_col}')
  
  # Move legend to lower right to avoid overlapping the text box
  plt.legend(loc='lower right')
  plt.tight_layout()

  # Save the file to the chosen folder
  filename = f'{x_col}_vs_{y_col}.png'
  filepath = os.path.join(target_dir, filename)
  plt.savefig(filepath)
  plt.close()
  
  print(f'Saved: {filename} -> {os.path.basename(target_dir)}')

# ====================================================================
# EXPORT METRICS TO CSV
# ====================================================================
# Export Sampled Metrics
if sampled_metrics:
  sampled_df = pd.DataFrame(sampled_metrics)
  sampled_csv_path = os.path.join(base_dir, 'Exact_vs_Sampled_Metrics.csv')
  sampled_df.to_csv(sampled_csv_path, index=False)
  print(f'\nSuccessfully exported Exact_vs_Sampled metrics to: {sampled_csv_path}')

# Export Exact Metrics
if exact_metrics:
  exact_df = pd.DataFrame(exact_metrics)
  exact_csv_path = os.path.join(base_dir, 'Exact_vs_Exact_Metrics.csv')
  exact_df.to_csv(exact_csv_path, index=False)
  print(f'Successfully exported Exact_vs_Exact metrics to: {exact_csv_path}')

print(f'\nAll plots and tables successfully generated inside: {base_dir}')
