import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

# Define output directory and ensure it exists
output_dir = r"C:\Users\andre\ScatterPlots"
os.makedirs(output_dir, exist_ok=True)

# Read the uploaded CSV
df = pd.read_csv(r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\SRER_NRI_Plots.csv")

# Define the pairs of columns to plot
pairs = [
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
  ('Exact_Gap_0_24', 'Exact_LPI_Pct'),
  ('Exact_Gap_25_50', 'Exact_LPI_Pct'),
  ('Exact_Gap_51_100', 'Exact_LPI_Pct'),
  ('Exact_Gap_101_200', 'Exact_LPI_Pct'),
  ('Exact_Gap_gt_200', 'Exact_LPI_Pct')
]

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
  # (index 0 is the constant/intercept, index 1 is the slope for X)
  intercept = model.params[0]
  slope = model.params[1]
  r_squared = model.rsquared
  
  # Retrieve predictions along with confidence and prediction intervals (alpha=0.05 -> 95%)
  pred = model.get_prediction(X_sm)
  pred_df = pred.summary_frame(alpha=0.05)
  
  # Plot Datapoints
  plt.scatter(X, Y, alpha=0.6, label='Data points', color='blue')
  
  # Plot Regression Line
  plt.plot(X, pred_df['mean'], color='red', label='Regression Line')
  
  # 95% Confidence Interval for the regression line
  plt.fill_between(X, pred_df['mean_ci_lower'], pred_df['mean_ci_upper'], 
                   color='red', alpha=0.3, label='95% Confidence Interval')
  
  # 95% Prediction Interval for individual data points
  plt.plot(X, pred_df['obs_ci_lower'], color='orange', linestyle='--', label='95% Prediction Interval')
  plt.plot(X, pred_df['obs_ci_upper'], color='orange', linestyle='--')
  
  # Add Slope, Intercept, and R-squared text box to the plot (top left)
  text_str = f'Slope: {slope:.3f}\nIntercept: {intercept:.3f}\n$R^2$: {r_squared:.3f}'
  plt.gca().text(0.05, 0.95, text_str, transform=plt.gca().transAxes, fontsize=10,
                 verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
  
  # Formatting the plot
  plt.xlabel(x_col)
  plt.ylabel(y_col)
  plt.title(f'{y_col} vs {x_col}')
  
  # Move legend to lower right to avoid overlapping the text box
  plt.legend(loc='lower right')
  plt.tight_layout()
  
  # Save the file to the new folder
  filename = f'{x_col}_vs_{y_col}.png'
  filepath = os.path.join(output_dir, filename)
  plt.savefig(filepath)
  plt.close()

print(f"All plots successfully generated and saved to: {output_dir}")
