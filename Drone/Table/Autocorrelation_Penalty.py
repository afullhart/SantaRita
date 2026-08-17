import pandas as pd
import numpy as np
import os
from scipy.optimize import curve_fit

# --- POWER LAW FUNCTION ---
def power_law(x, a, b, c):
    return a * np.power(x, -b) + c

# --- HELPER TO EXTRACT MAE AND FIT CURVE ---
def get_random_fit(df, exact_col, pt_col_prefix, pts):
    mae_vals = []
    for p in pts:
        mae = np.abs(df[f'{pt_col_prefix}{p}'] - df[exact_col]).mean()
        mae_vals.append(mae)
    
    x = np.array(pts, dtype=float)
    y = np.array(mae_vals, dtype=float)
    
    # Fit only on points >= 5 to match your established methodology
    mask = x >= 5
    x_fit = x[mask]
    y_fit = y[mask]
    
    p0 = [np.max(y_fit), 0.5, np.min(y_fit)]
    bounds = ([-np.inf, 0.01, -np.inf], [np.inf, 5.0, np.inf])
    
    params, _ = curve_fit(power_law, x_fit, y_fit, p0=p0, bounds=bounds, maxfev=10000)
    return params

def main():
    data_folder = r'C:\Users\andre\ScatterPlots'
    
    line_csv = os.path.join(data_folder, 'SRER_NRI_Plots_110m.csv')
    rand_csv = os.path.join(data_folder, 'SRER_NRI_Plots_110m_Random.csv')
    
    if not os.path.exists(line_csv) or not os.path.exists(rand_csv):
        print("Data files not found. Please verify the path.")
        return

    df_line = pd.read_csv(line_csv)
    df_rand = pd.read_csv(rand_csv)
    
    # Random sampling densities to fit the curve
    pts = [1, 3, 5, 7, 10, 15, 20, 30, 50, 70, 100, 200, 300, 500, 700, 1000, 10000]
    
    # Fit the convergence curves for the Random point sampling
    params_bgr = get_random_fit(df_rand, 'BGR_Exact', 'BGR_pt_', pts)
    params_herb = get_random_fit(df_rand, 'Herb_Pct_Exact', 'Herb_Pct_pt_', pts)
    params_woody = get_random_fit(df_rand, 'Woody_Pct_Exact', 'Woody_Pct_pt_', pts)
    
    # Print the derived equations to the console for confirmation
    print("\n--- DERIVED POWER LAW EQUATIONS (110m NRI Base) ---")
    print(f"BGR:   y = {params_bgr[0]:.2f}x^-{params_bgr[1]:.2f} + {params_bgr[2]:.2f}")
    print(f"Herb:  y = {params_herb[0]:.2f}x^-{params_herb[1]:.2f} + {params_herb[2]:.2f}")
    print(f"Woody: y = {params_woody[0]:.2f}x^-{params_woody[1]:.2f} + {params_woody[2]:.2f}")
    
    # Define the line-intercept spacings and their equivalent random point count (N)
    # Based on a 150m total transect length per 110m plot (150m / spacing)
    intervals = {
        '200cm': 75, 
        '100cm': 150, 
        '50cm': 300, 
        '25cm': 600
    }
    
    results = []
    
    for label, n in intervals.items():
        # 1. Calculate actual Transect MAE directly from the line file
        t_bgr = np.abs(df_line[f'NRI_BGR_{label}_Pct'] - df_line['Exact_BGR_Pct']).mean()
        t_herb = np.abs(df_line[f'NRI_Herb_{label}_Pct'] - df_line['Exact_Herb_Pct']).mean()
        t_woody = np.abs(df_line[f'NRI_Woody_{label}_Pct'] - df_line['Exact_Woody_Pct']).mean()
        
        # 2. Calculate theoretical Random MAE from the fitted equations
        r_bgr = power_law(n, *params_bgr)
        r_herb = power_law(n, *params_herb)
        r_woody = power_law(n, *params_woody)
        
        # 3. Compile the row data
        results.append({
            'Spacing': label,
            'Equiv. N': n,
            'BGR (Rand)': r_bgr,
            'BGR (Line)': t_bgr,
            'BGR Penalty': t_bgr - r_bgr,
            'Herb (Rand)': r_herb,
            'Herb (Line)': t_herb,
            'Herb Penalty': t_herb - r_herb,
            'Woody (Rand)': r_woody,
            'Woody (Line)': t_woody,
            'Woody Penalty': t_woody - r_woody
        })
        
    # Build and display DataFrame
    df_results = pd.DataFrame(results)
    
    # Format the table for neat console printing
    pd.set_option('display.float_format', '{:.2f}'.format)
    print("\n--- SPATIAL AUTOCORRELATION PENALTY TABLE (110m NRI Base) ---\n")
    print(df_results.to_string(index=False))

    # Save to CSV
    out_folder = r'C:\Users\andre\ScatterPlots'
    out_path = os.path.join(out_folder, 'Autocorrelation_Penalty_Table.csv')
    df_results.to_csv(out_path, index=False)
    print(f"\nSaved table to: {out_path}")

if __name__ == '__main__':
    main()
