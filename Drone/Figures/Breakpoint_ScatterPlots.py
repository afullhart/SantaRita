import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

# ====================================================================
# CONFIGURATION & PUBLICATION STYLE
# ====================================================================
output_dir = r'C:\Users\andre\ScatterPlots\Piecewise_Regression'
os.makedirs(output_dir, exist_ok=True)

csv_path = r'C:\Users\andre\ScatterPlots\SRER_NRI_Plots_110m.csv'
df = pd.read_csv(csv_path)

# Apply publication-ready global styling (Matched to Figure 14)
plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'axes.linewidth': 1.2,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans']
})

# Define styling constants for consistency across all plots
SCATTER_COLOR = '#4C72B0'  # Muted seaborn-style blue
SCATTER_ALPHA = 0.8        # Increased opacity
SCATTER_SIZE = 60
LINE_COLOR = '#c44e52'     # Muted crimson
LINE_WIDTH = 2.5
VLINE_COLOR = '#555555'

# ====================================================================
# PIECEWISE FUNCTION & FITTING LOGIC
# ====================================================================
def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, 
                        [x < x0], 
                        [lambda x: k1*x + y0 - k1*x0, 
                         lambda x: k2*x + y0 - k2*x0])

def get_piecewise_fit(data, x_col, y_col, initial_guess):
    subset = data[[x_col, y_col]].dropna().sort_values(by=x_col)
    
    if subset.empty:
        raise ValueError(f"No valid data found for {y_col} vs {x_col}")
        
    x = subset[x_col].values
    y = subset[y_col].values

    params, _ = curve_fit(piecewise_linear, x, y, p0=initial_guess, maxfev=10000)
    y_pred = piecewise_linear(x, *params)
    r2 = r2_score(y, y_pred)
    
    return x, y, y_pred, params, r2

# ====================================================================
# PLOTTING FUNCTIONS
# ====================================================================
def style_axes(ax):
    """Applies standard despining and grid styling to a given axis."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='both', linestyle='-', alpha=0.2)

def generate_individual_plot(data, x_col, y_col, initial_guess):
    print(f"Generating single plot for {y_col} vs {x_col}...")
    
    try:
        x, y, y_pred, params, r2 = get_piecewise_fit(data, x_col, y_col, initial_guess)
    except Exception as e:
        print(f"  -> Skipping. Error during fitting: {e}")
        return

    x0, y0, k1, k2 = params

    fig, ax = plt.subplots(figsize=(8, 7))

    ax.scatter(x, y, alpha=SCATTER_ALPHA, color=SCATTER_COLOR, edgecolor='black', linewidth=0.8, s=SCATTER_SIZE, label='Plot Data')
    ax.plot(x, y_pred, color=LINE_COLOR, linewidth=LINE_WIDTH, label=f'Piecewise Fit ($R^2$={r2:.2f})')
    
    # Label simply as "Breakpoint"
    ax.axvline(x=x0, color=VLINE_COLOR, linestyle='--', linewidth=1.5, zorder=0, label='Breakpoint')

    stats_text = f"Slope 1: {k1:.2f}\nSlope 2: {k2:.2f}\nBreakpoint: {x0:.1f}"
    props = dict(boxstyle='round', facecolor='white', alpha=0.3, edgecolor='gray')
    ax.text(0.04, 0.96, stats_text, transform=ax.transAxes, fontsize=17,
            verticalalignment='top', horizontalalignment='left', bbox=props)

    # Format labels (Replace underscores for readability)
    ax.set_xlabel(x_col.replace('_', ' '), weight='bold')
    ax.set_ylabel(y_col.replace('_', ' '), weight='bold')
    ax.set_title(f'{y_col.replace("_", " ")} vs {x_col.replace("_", " ")}', fontsize=20, pad=15)
    ax.legend(loc='lower right', frameon=True, fontsize=16, framealpha=0.3)
    
    style_axes(ax)
    plt.tight_layout()
    
    filename = f"Piecewise_{x_col}_vs_{y_col}.svg"
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, format='svg', bbox_inches='tight')
    plt.close()
    print(f"  -> Saved successfully: {filename}")

def generate_1x2_publication_subplot(data):
    print("\nGenerating Publication-Ready 1x2 LPI Subplot...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Shared text box styling
    props = dict(boxstyle='round', facecolor='white', alpha=0.3, edgecolor='gray')

    # --- Subplot 1 (Left): BGR vs LPI ---
    guess_bgr = [25.0, 5.0, 0.2, 1.2]
    x1, y1, y_pred1, params1, r2_1 = get_piecewise_fit(data, 'Exact_BGR_Pct', 'Exact_LPI_Pct', guess_bgr)
    x0_1, y0_1, k1_1, k2_1 = params1

    axes[0].scatter(x1, y1, alpha=SCATTER_ALPHA, color=SCATTER_COLOR, edgecolor='black', linewidth=0.8, s=SCATTER_SIZE, label='Plot Data')
    axes[0].plot(x1, y_pred1, color=LINE_COLOR, linewidth=LINE_WIDTH, label=f'Piecewise Fit ($R^2$={r2_1:.2f})')
    
    # Label simply as "Breakpoint"
    axes[0].axvline(x=x0_1, color=VLINE_COLOR, linestyle='--', linewidth=1.5, zorder=0, label='Breakpoint')

    stats_text1 = f"Slope 1: {k1_1:.2f}\nSlope 2: {k2_1:.2f}\nBreakpoint: {x0_1:.1f}%"
    axes[0].text(0.04, 0.96, stats_text1, transform=axes[0].transAxes, fontsize=17,
                 verticalalignment='top', horizontalalignment='left', bbox=props)

    axes[0].set_xlabel('Total Bare Ground (%)', weight='bold')
    axes[0].set_ylabel('Largest Patch Index (%)', weight='bold')
    
    # Split Title Alignment
    axes[0].set_title(r'$\bf{a.}$', loc='left', fontsize=20, pad=15)
    axes[0].set_title('Bare Ground Coalescence Threshold', loc='center', fontsize=20, pad=15)
    
    axes[0].legend(loc='lower right', frameon=True, fontsize=16, framealpha=0.3)
    style_axes(axes[0])

    # --- Subplot 2 (Right): Gap 101-200 vs LPI ---
    gap_col = 'Exact_Gap_101_200'
    gap_mean = data[gap_col].mean()
    lpi_mean = data['Exact_LPI_Pct'].mean()
    guess_gap = [gap_mean, lpi_mean, 0.0, 0.0]

    x2, y2, y_pred2, params2, r2_2 = get_piecewise_fit(data, gap_col, 'Exact_LPI_Pct', guess_gap)
    x0_2, y0_2, k1_2, k2_2 = params2

    axes[1].scatter(x2, y2, alpha=SCATTER_ALPHA, color=SCATTER_COLOR, edgecolor='black', linewidth=0.8, s=SCATTER_SIZE, label='Plot Data')
    axes[1].plot(x2, y_pred2, color=LINE_COLOR, linewidth=LINE_WIDTH, label=f'Piecewise Fit ($R^2$={r2_2:.2f})')
    
    # Label simply as "Breakpoint"
    axes[1].axvline(x=x0_2, color=VLINE_COLOR, linestyle='--', linewidth=1.5, zorder=0, label='Breakpoint')

    stats_text2 = f"Slope 1: {k1_2:.2f}\nSlope 2: {k2_2:.2f}\nBreakpoint: {x0_2:.1f}%"
    axes[1].text(0.04, 0.96, stats_text2, transform=axes[1].transAxes, fontsize=17,
                 verticalalignment='top', horizontalalignment='left', bbox=props)

    axes[1].set_xlabel('Canopy Gap 101-200 cm (%)', weight='bold')
    axes[1].set_ylabel('Largest Patch Index (%)', weight='bold')
    
    # Split Title Alignment
    axes[1].set_title(r'$\bf{b.}$', loc='left', fontsize=20, pad=15)
    axes[1].set_title('Canopy Gap Coalescence Threshold', loc='center', fontsize=20, pad=15)
    
    axes[1].legend(loc='lower right', frameon=True, fontsize=16, framealpha=0.3)
    style_axes(axes[1])

    plt.tight_layout()
    filename = "Publication_Piecewise_Subplots_LPI.svg"
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, format='svg', bbox_inches='tight')
    plt.close()
    print(f"  -> Saved successfully: {filename}")

# ====================================================================
# EXECUTION
# ====================================================================
print("\n--- Generating Individual Piecewise Regressions ---")

# Setup common means for initial guesses
hw_mean = df['Exact_Herb_Woody_Ratio'].mean()
lpi_mean = df['Exact_LPI_Pct'].mean()

# 1. BGR vs LPI
generate_individual_plot(df, 'Exact_BGR_Pct', 'Exact_LPI_Pct', initial_guess=[25.0, 5.0, 0.2, 1.2])

# 2. BGR vs HW Ratio
generate_individual_plot(df, 'Exact_BGR_Pct', 'Exact_Herb_Woody_Ratio', initial_guess=[25.0, hw_mean, 0.0, 0.0])

# 3. LPI vs HW Ratio
generate_individual_plot(df, 'Exact_LPI_Pct', 'Exact_Herb_Woody_Ratio', initial_guess=[25.0, hw_mean, 0.0, 0.0])

# 4-8. Gap Categories vs LPI
gap_cols = [
  'Exact_Gap_0_24',
  'Exact_Gap_25_50',
  'Exact_Gap_51_100',
  'Exact_Gap_101_200',
  'Exact_Gap_gt_200'
]

for gap_col in gap_cols:
    gap_mean = df[gap_col].mean()
    generate_individual_plot(df, gap_col, 'Exact_LPI_Pct', initial_guess=[gap_mean, lpi_mean, 0.0, 0.0])

# 9. Mean Fetch vs LPI
if 'Exact_Fetch_m' in df.columns:
    fetch_mean = df['Exact_Fetch_m'].mean()
    generate_individual_plot(df, 'Exact_Fetch_m', 'Exact_LPI_Pct', initial_guess=[fetch_mean, lpi_mean, 0.0, 0.0])
else:
    print("  -> Column 'Exact_Fetch_m' not found in CSV. Skipping Mean Fetch vs LPI plot.")

# 10. Generate the final 1x2 publication-ready subplot
generate_1x2_publication_subplot(df)

print(f'\nAll operations complete! All SVG plots are located in: {output_dir}')
