import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

# ====================================================================
# CONFIGURATION
# ====================================================================
output_dir = r"C:\Users\andre\ScatterPlots"
os.makedirs(output_dir, exist_ok=True)

csv_path = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\SRER_NRI_Plots.csv"
df = pd.read_csv(csv_path)

# ====================================================================
# PIECEWISE FUNCTION & PLOTTING LOGIC
# ====================================================================
def piecewise_linear(x, x0, y0, k1, k2):
  """
  x0 = breakpoint (x-coordinate)
  y0 = value at breakpoint (y-coordinate)
  k1 = slope of the first segment
  k2 = slope of the second segment
  """
  return np.piecewise(x, 
                      [x < x0], 
                      [lambda x: k1*x + y0 - k1*x0, 
                       lambda x: k2*x + y0 - k2*x0])

def generate_piecewise_plot(data, x_col, y_col, initial_guess):
  print(f"Fitting piecewise regression for {y_col} vs {x_col}...")
  
  # Drop NaNs specifically for the two columns being plotted and sort by X
  subset = data[[x_col, y_col]].dropna().sort_values(by=x_col)
  
  if subset.empty:
    print(f"  -> Skipping: No valid data found for {y_col}")
    return
    
  x = subset[x_col].values
  y = subset[y_col].values

  # Fit the model using scipy's curve_fit
  try:
    # Maxfev increased to give the optimizer more attempts to converge on tricky distributions
    params, covariance = curve_fit(piecewise_linear, x, y, p0=initial_guess, maxfev=10000)
  except Exception as e:
    print(f"  -> Optimizer failed to find a fit: {e}")
    return

  x0, y0, k1, k2 = params

  # Generate predictions and calculate R-squared
  y_pred = piecewise_linear(x, *params)
  r2 = r2_score(y, y_pred)

  # Create the plot
  plt.figure(figsize=(10, 6))

  # Plot original data points
  plt.scatter(x, y, alpha=0.6, color='blue', label='Data points')

  # Plot the fitted piecewise line
  plt.plot(x, y_pred, color='red', linewidth=2, label=f'Piecewise Fit ($R^2$={r2:.3f})')

  # Add a vertical dashed line at the breakpoint
  plt.axvline(x=x0, color='gray', linestyle='--', label=f'Breakpoint: {x0:.1f}')

  # Add a text box with the slope and breakpoint information
  stats_text = f"Slope 1: {k1:.3f}\nSlope 2: {k2:.3f}\nBreakpoint: {x0:.1f}"
  plt.text(0.05, 0.90, stats_text, 
           transform=plt.gca().transAxes,
           bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))

  # Labels, Title, and formatting
  plt.xlabel(x_col)
  plt.ylabel(y_col)
  plt.title(f'Piecewise Regression: {y_col} vs {x_col}')
  plt.legend(loc='lower right')
  plt.grid(True, alpha=0.3)

  # Save the plot directly to the target folder
  plt.tight_layout()
  filename = f"Piecewise_{x_col}_vs_{y_col}.png"
  filepath = os.path.join(output_dir, filename)
  plt.savefig(filepath, dpi=300)
  plt.close()
  
  print(f"  -> Saved successfully: {filename}")

# ====================================================================
# EXECUTION
# ====================================================================
print("\n--- Generating Piecewise Regressions ---")

# Plot 1: Exact BGR vs Exact LPI
guess_lpi = [25.0, 5.0, 0.2, 1.2]
generate_piecewise_plot(df, 'Exact_BGR_Pct', 'Exact_LPI_Pct', initial_guess=guess_lpi)

# Plot 2: Exact BGR vs Exact Herb-to-Woody Ratio
hw_mean = df['Exact_Herb_Woody_Ratio'].mean()
guess_hw = [25.0, hw_mean, 0.0, 0.0]
generate_piecewise_plot(df, 'Exact_BGR_Pct', 'Exact_Herb_Woody_Ratio', initial_guess=guess_hw)

# Plot 3: Exact LPI vs Exact Herb-to-Woody Ratio
guess_lpi_hw = [25.0, hw_mean, 0.0, 0.0]
generate_piecewise_plot(df, 'Exact_LPI_Pct', 'Exact_Herb_Woody_Ratio', initial_guess=guess_lpi_hw)

# Plots 4-8: Exact Canopy Gap Categories vs Exact LPI
gap_cols = [
  'Exact_Gap_0_24',
  'Exact_Gap_25_50',
  'Exact_Gap_51_100',
  'Exact_Gap_101_200',
  'Exact_Gap_gt_200'
]

lpi_mean = df['Exact_LPI_Pct'].mean()

for gap_col in gap_cols:
  # Dynamically anchor the breakpoint to the specific mean of the current gap category
  gap_mean = df[gap_col].mean()
  guess_gap = [gap_mean, lpi_mean, 0.0, 0.0]
  
  # Plotting Gap fraction on X-axis, LPI on Y-axis
  generate_piecewise_plot(df, gap_col, 'Exact_LPI_Pct', initial_guess=guess_gap)

print(f"\nAll operations complete! Plots are located in: {output_dir}")
