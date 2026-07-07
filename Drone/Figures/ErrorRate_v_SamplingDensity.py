import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

data_folder = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data'
out_folder = r'C:\Users\andre\ScatterPlots\Random_Sampling'

csv_files = {
  '10m_Grid': 'SRER_Grid_10m_Random.csv',
  '30m_Grid': 'SRER_Grid_30m_Random.csv',
  '110m_NRI': 'SRER_NRI_Plots_110m_Random.csv'
}

# --- HELPER FUNCTION FOR PLOTTING ---
def plot_dual_axis(ax, x_vals, mae_errs, rel_errs, title, x_label, y1_label):
  color_mae = 'tab:red'
  color_rel = 'tab:blue'
  
  ax.set_title(title, fontsize=14, fontweight='bold')
  ax.set_xlabel(x_label, fontsize=12)
  
  # Left Axis: Absolute Error
  ax.set_ylabel(y1_label, color=color_mae, fontsize=12)
  ax.plot(x_vals, mae_errs, marker='s', linestyle='-', color=color_mae, label='Absolute Error')
  ax.tick_params(axis='y', labelcolor=color_mae)
  
  # Format X-Axis
  ax.set_xscale('log')
  ax.set_xticks(x_vals)
  ax.set_xticklabels([str(x) for x in x_vals], rotation=45)
  ax.grid(True, which="major", ls="--", alpha=0.5)

  # Right Axis: Relative Error
  ax2 = ax.twinx()
  ax2.set_ylabel('Mean Relative % Error', color=color_rel, fontsize=12)
  ax2.plot(x_vals, rel_errs, marker='o', linestyle='-', color=color_rel, label='Relative % Error')
  ax2.tick_params(axis='y', labelcolor=color_rel)

for scale_name, csv_filename in csv_files.items():
  csv_path = os.path.join(data_folder, csv_filename)
  
  if not os.path.exists(csv_path):
    print(f"Skipping {scale_name}: Could not find {csv_path}")
    continue
    
  print(f"\n--- Generating Plots for {scale_name} ---")
  df = pd.read_csv(csv_path)

  # 1. Dynamically find sampling increments
  pt_cols = [c for c in df.columns if c.startswith('BGR_pt_')]
  if not pt_cols:
    print("No point sampling columns found. Skipping.")
    continue
  points = sorted([int(c.split('_')[-1]) for c in pt_cols])

  ln_cols = [c for c in df.columns if c.startswith('Gap_0_24_L_')]
  lines = sorted([int(c.split('_')[-1]) for c in ln_cols])

  # 2. Calculate errors for BGR, Herb, Woody, Fetch, and Herb-to-Woody Ratio
  bgr_mae, bgr_rel = [], []
  herb_mae, herb_rel = [], []
  woody_mae, woody_rel = [], []
  fetch_mae, fetch_rel = [], []
  hw_mae, hw_rel = [], []

  for pt in points:
    # BGR
    b_col = f'BGR_pt_{pt}'
    bgr_mae.append(np.abs(df[b_col] - df['BGR_Exact']).mean())
    bgr_rel.append((np.abs(df[b_col] - df['BGR_Exact']) / (df['BGR_Exact'] + 1e-6) * 100).mean())
    
    # Herb Pct
    h_col = f'Herb_Pct_pt_{pt}'
    if h_col in df.columns and 'Herb_Pct_Exact' in df.columns:
      herb_mae.append(np.abs(df[h_col] - df['Herb_Pct_Exact']).mean())
      herb_rel.append((np.abs(df[h_col] - df['Herb_Pct_Exact']) / (df['Herb_Pct_Exact'] + 1e-6) * 100).mean())
    else:
      herb_mae.append(np.nan)
      herb_rel.append(np.nan)

    # Woody Pct
    w_col = f'Woody_Pct_pt_{pt}'
    if w_col in df.columns and 'Woody_Pct_Exact' in df.columns:
      woody_mae.append(np.abs(df[w_col] - df['Woody_Pct_Exact']).mean())
      woody_rel.append((np.abs(df[w_col] - df['Woody_Pct_Exact']) / (df['Woody_Pct_Exact'] + 1e-6) * 100).mean())
    else:
      woody_mae.append(np.nan)
      woody_rel.append(np.nan)

    # Fetch
    f_col = f'Fetch_pt_{pt}'
    fetch_mae.append(np.abs(df[f_col] - df['Fetch_Exact']).mean())
    fetch_rel.append((np.abs(df[f_col] - df['Fetch_Exact']) / (df['Fetch_Exact'] + 1e-6) * 100).mean())

    # Herb-to-Woody Ratio
    hw_col = f'HW_Ratio_pt_{pt}'
    if hw_col in df.columns and 'HW_Ratio_Exact' in df.columns:
      hw_mae.append(np.abs(df[hw_col] - df['HW_Ratio_Exact']).mean())
      hw_rel.append((np.abs(df[hw_col] - df['HW_Ratio_Exact']) / (df['HW_Ratio_Exact'] + 1e-6) * 100).mean())
    else:
      hw_mae.append(np.nan)
      hw_rel.append(np.nan)

  # 3. Calculate errors for the 5 Canopy Gap Categories
  gap_cats = {
    'Gap_0_24': 'Canopy Gap: 0-24 cm',
    'Gap_25_50': 'Canopy Gap: 25-50 cm',
    'Gap_51_100': 'Canopy Gap: 51-100 cm',
    'Gap_101_200': 'Canopy Gap: 101-200 cm',
    'Gap_gt_200': 'Canopy Gap: >200 cm'
  }
  
  gap_data = {}
  for cat, title in gap_cats.items():
    g_mae, g_rel = [], []
    exact_col = f'{cat}_Exact'
    
    for ln in lines:
      l_col = f'{cat}_L_{ln}'
      g_mae.append(np.abs(df[l_col] - df[exact_col]).mean())
      g_rel.append((np.abs(df[l_col] - df[exact_col]) / (df[exact_col] + 1e-6) * 100).mean())
      
    gap_data[cat] = {'title': title, 'mae': g_mae, 'rel': g_rel}

  # 4. Generate the 10-Subplot Figure
  fig, axes = plt.subplots(nrows=10, ncols=1, figsize=(12, 50))

  # Plot 1: BGR
  plot_dual_axis(axes[0], points, bgr_mae, bgr_rel, 
                 'Bare Ground Percentage', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (pp)')
  
  # Plot 2: Herb Cover Percentage
  plot_dual_axis(axes[1], points, herb_mae, herb_rel, 
                 'Herb Cover Percentage', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (pp)')

  # Plot 3: Woody Cover Percentage
  plot_dual_axis(axes[2], points, woody_mae, woody_rel, 
                 'Woody Cover Percentage', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (pp)')
  
  # Plot 4: Fetch
  plot_dual_axis(axes[3], points, fetch_mae, fetch_rel, 
                 'Mean Fetch', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (m)')

  # Plot 5: Herb-to-Woody Ratio
  plot_dual_axis(axes[4], points, hw_mae, hw_rel, 
                 'Herb-to-Woody Ratio', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (Ratio)')

  # Plots 6-10: Canopy Gaps
  ax_idx = 5
  for cat in gap_cats.keys():
    plot_dual_axis(axes[ax_idx], lines, gap_data[cat]['mae'], gap_data[cat]['rel'], 
                   gap_data[cat]['title'], 'Virtual Transect Length in Meters (Log Scale)', 'Absolute Error (pp)')
    ax_idx += 1

  # 5. Final Formatting and Export
  fig.suptitle(f'Convergence of Sampled Metrics to Exact Values ({scale_name.replace("_", " ")})', 
               fontsize=20, fontweight='bold', y=0.99)
               
  # Generous padding to prevent overlapping labels
  fig.tight_layout(pad=4.0)  

  output_img_path = os.path.join(out_folder, f'{scale_name}_Convergence.png')
  plt.savefig(output_img_path, dpi=300)
  plt.close() 
  
  print(f"Saved comprehensive 10-plot figure to: {output_img_path}")

print("\nAll comprehensive plotting complete!")
