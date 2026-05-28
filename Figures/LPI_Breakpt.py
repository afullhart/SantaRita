import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

# 1. Load and sort the data
# Make sure 'ok.csv' is in your current working directory
df = pd.read_csv(r'C:\Users\andre\Desktop\data.csv')
df = df.sort_values(by='Exact_BGR_Pct')

x = df['Exact_BGR_Pct'].values
y = df['Exact_LPI_Pct'].values

# 2. Define the piecewise linear function
# x0 = breakpoint (x-coordinate)
# y0 = value at breakpoint (y-coordinate)
# k1 = slope of the first segment
# k2 = slope of the second segment
def piecewise_linear(x, x0, y0, k1, k2):
  return np.piecewise(x, 
                      [x < x0], 
                      [lambda x: k1*x + y0 - k1*x0, 
                       lambda x: k2*x + y0 - k2*x0])

# 3. Fit the model
# Provide an initial guess [x0, y0, k1, k2] to help the optimizer
initial_guess = [25.0, 5.0, 0.2, 1.2]

# curve_fit finds the best parameters
params, covariance = curve_fit(piecewise_linear, x, y, p0=initial_guess)
x0, y0, k1, k2 = params

# 4. Generate predictions and calculate R-squared
y_pred = piecewise_linear(x, *params)
r2 = r2_score(y, y_pred)

# 5. Create the plot
plt.figure(figsize=(10, 6))

# Plot original data points
plt.scatter(x, y, alpha=0.6, color='blue', label='Data points')

# Plot the fitted piecewise line
plt.plot(x, y_pred, color='red', linewidth=2, label=f'Piecewise Fit ($R^2$={r2:.3f})')

# Add a vertical dashed line at the breakpoint
plt.axvline(x=x0, color='gray', linestyle='--', label=f'Breakpoint: {x0:.1f}%')

# Add a text box with the slope and breakpoint information
stats_text = f"Slope 1: {k1:.3f}\nSlope 2: {k2:.3f}\nBreakpoint: {x0:.1f}%"
plt.text(0.05, 0.90, stats_text, 
         transform=plt.gca().transAxes,
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))

# Labels, Title, and formatting
plt.xlabel('Exact_BGR_Pct')
plt.ylabel('Exact_LPI_Pct')
plt.title('Piecewise Regression: BGR vs LPI')
plt.legend(loc='lower right')
plt.grid(True, alpha=0.3)

# Show or save the plot
plt.tight_layout()
plt.show()
# plt.savefig('piecewise_regression.png', dpi=300) # Use this to save to file
