import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import os

# ====================================================================
# CONFIGURATION
# ====================================================================
data_folder = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1'
out_folder = r'C:\Users\andre\ScatterPlots\Spoke_Convergence'
os.makedirs(out_folder, exist_ok=True)

csv_files = {
    '110m_NRI': 'SRER_NRI_Plots_110m_Spoke_Convergence.csv',
    '30m_Grid': 'SRER_Grid_30m_Spoke_Convergence.csv',
    '10m_Grid': 'SRER_Grid_10m_Spoke_Convergence.csv'
}

NRI_INTERVALS_CM = [500, 250, 200, 100, 50, 25, 10, 5]

# ====================================================================
# HELPER FUNCTIONS
# ====================================================================
def get_ols_metrics(df, sampled_col, exact_col):
    """Calculates the Diagnostic Slope (Predicting Sample from Truth) and R-Squared"""
    subset = df[[sampled_col, exact_col]].dropna()
    if len(subset) < 2: 
        return np.nan, np.nan
    
    # X is the Exact baseline (assumed truth)
    X = sm.add_constant(subset[exact_col].values)
    # Y is the Sampled data
    model = sm.OLS(subset[sampled_col].values, X).fit()
    
    return model.params[1], model.rsquared

def plot_slope_r2(ax, x_vals, slopes_diag, r2_vals, title):
    """Plots the Diagnostic Slope and R-Squared on dual axes with standardized limits"""
    color_d = 'tab:green'
    color_r2 = 'tab:orange'
    
    ax.set_title(title, fontsize=12, fontweight='bold')
    
    # Left Axis: Diagnostic Slope
    ax.set_ylabel('Diagnostic Slope\n(Sampled ~ Exact)', color='black', fontsize=10)
    ax.plot(x_vals, slopes_diag, marker='^', linestyle='-', color=color_d, label='Diagnostic Slope')
    ax.tick_params(axis='y', labelcolor='black')
    ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, label='Perfect 1.0')

    ax.set_ylim(0.5, 1.5)
    
    # X-Axis Formatting
    ax.set_xlabel('Number of Sampled Pins (Log Scale)')
    ax.set_xscale('log')
    ax.set_xticks(x_vals)
    ax.set_xticklabels([str(int(x)) for x in x_vals], rotation=45)
    ax.grid(True, ls="--", alpha=0.5)
    
    # Right Axis: R-Squared
    ax2 = ax.twinx()
    ax2.set_ylabel('R-Squared ($R^2$)', color=color_r2, fontsize=10)
    ax2.plot(x_vals, r2_vals, marker='s', linestyle='-', color=color_r2, alpha=0.8, label='$R^2$')
    ax2.tick_params(axis='y', labelcolor=color_r2)
    
    # --- APPLY FIXED R2 Y-LIMITS ---
    ax2.set_ylim(0.5, 1.0)

    # Unified Legend
    lines_1, labels_1 = ax.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax2.legend(lines_1 + lines_2, labels_1 + labels_2, loc='center right', fontsize=9)

# ====================================================================
# MAIN PLOTTING LOOP
# ====================================================================
if __name__ == '__main__':
    for scale_name, csv_filename in csv_files.items():
        csv_path = os.path.join(data_folder, csv_filename)
        
        if not os.path.exists(csv_path): 
            print(f"Skipping {scale_name}: Could not find {csv_path}")
            continue
            
        df = pd.read_csv(csv_path)
        
        # Extract the average number of pins for each interval to use as X-axis
        pin_counts = [df[f'NRI_TotalPins_{inv}cm'].mean() for inv in NRI_INTERVALS_CM]
        
        metrics = {
            'Bare Ground Percentage': ('NRI_BGR_{}cm_Pct', 'Exact_BGR_Pct'),
            'Herb Cover Percentage': ('NRI_Herb_{}cm_Pct', 'Exact_Herb_Pct'),
            'Woody Cover Percentage': ('NRI_Woody_{}cm_Pct', 'Exact_Woody_Pct'),
            'Herb-to-Woody Ratio': ('NRI_HW_Ratio_{}cm', 'Exact_HW_Ratio'),
            'Mean Fetch': ('NRI_Fetch_{}cm', 'Exact_Fetch_m')
        }

        fig, axes = plt.subplots(5, 1, figsize=(10, 25))
        
        for ax, (title, (x_template, y_col)) in zip(axes, metrics.items()):
            slopes_d, r2s = [], []
            
            for inv in NRI_INTERVALS_CM:
                s_d, r = get_ols_metrics(df, x_template.format(inv), y_col)
                slopes_d.append(s_d)
                r2s.append(r)
            
            plot_slope_r2(ax, pin_counts, slopes_d, r2s, title)

        fig.suptitle(f'Diagnostic Spoke Convergence ({scale_name})', fontsize=16, fontweight='bold', y=0.99)
        fig.tight_layout(pad=3.0)  
        
        # Save the plot
        output_img_path = os.path.join(out_folder, f'{scale_name}_Diagnostic_Slope_Convergence.png')
        plt.savefig(output_img_path, dpi=300)
        plt.close() 
        
        print(f"Saved: {os.path.basename(output_img_path)}")

    print("\nAll standardized diagnostic plotting complete!")
