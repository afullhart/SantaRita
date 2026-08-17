import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

data_folder = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data'
out_folder = r'C:\Users\andre\ScatterPlots\Random_Sampling'

csv_files = {
  '10m_Grid': 'SRER_Grid_10m_Random.csv',
  '30m_Grid': 'SRER_Grid_30m_Random.csv',
  '110m_NRI': 'SRER_NRI_Plots_110m_Random.csv'
}

fit_results = {'10m_Grid': {}, '30m_Grid': {}, '110m_NRI': {}}
fit_params = {'10m_Grid': {}, '30m_Grid': {}, '110m_NRI': {}}
data_cache = {}

# --- POWER LAW FUNCTION ---
def power_law(x, a, b, c):
    return a * np.power(x, -b) + c

# --- HELPER FUNCTION FOR STANDARD PLOTTING ---
def plot_error_convergence(ax, x_vals, mae_errs, rel_errs, title, x_label, y_label, combined_fits=None, plot_secondary=True, extra_curves=None, mae_label='Absolute Error', zoom_to_fit=False):
  color_mae = 'tab:red'
  color_rel = 'tab:blue'
  color_fit_main = 'black'
  
  ax.set_title(title, fontsize=18, fontweight='bold')
  ax.set_xlabel(x_label, fontsize=16)
  
  # Left Axis: Absolute Error
  ax.set_ylabel(y_label, color=color_mae, fontsize=16)
  
  l1 = ax.plot(x_vals, mae_errs, marker='s', linestyle='-', color=color_mae, alpha=0.4, label=mae_label)
  ax.tick_params(axis='y', labelcolor=color_mae, labelsize=14)
  ax.tick_params(axis='x', labelsize=14)
  
  return_data = {'eq_text': "Fit Failed", 'params': None}
  
  # --- POWER LAW FITTING ($x \ge 5$) ---
  x_arr = np.array(x_vals, dtype=float)
  y_arr = np.array(mae_errs, dtype=float)
  
  valid_mask = ~np.isnan(y_arr) & ~np.isnan(x_arr)
  x_valid = x_arr[valid_mask]
  y_valid = y_arr[valid_mask]
  
  fit_mask = x_valid >= 5
  x_fit = x_valid[fit_mask]
  y_fit = y_valid[fit_mask]
  
  l2, l_extra = [], []
  y_min_zoom, y_max_zoom = np.inf, -np.inf
  
  if len(y_fit) > 0:
      y_min_zoom = min(y_min_zoom, np.min(y_fit))
      y_max_zoom = max(y_max_zoom, np.max(y_fit))
  
  if len(x_fit) >= 3:
      try:
          p0 = [np.max(y_fit), 0.5, np.min(y_fit)]
          bounds = ([-np.inf, 0.01, -np.inf], [np.inf, 5.0, np.inf])
          
          params, _ = curve_fit(power_law, x_fit, y_fit, p0=p0, bounds=bounds, maxfev=10000)
          a, b, c = params
          
          x_smooth = np.geomspace(np.min(x_fit), np.max(x_fit), 200)
          y_smooth = power_law(x_smooth, a, b, c)
          
          main_fit_label = r'Power Law Fit ($x \geq 5$)' if combined_fits is None else r'Fit: 110m NRI'
          l2 = ax.plot(x_smooth, y_smooth, linestyle='--', color=color_fit_main, linewidth=2, label=main_fit_label)
          
          y_pred = power_law(x_fit, a, b, c)
          ss_res = np.sum((y_fit - y_pred) ** 2)
          ss_tot = np.sum((y_fit - np.mean(y_fit)) ** 2)
          r2 = 1 - (ss_res / ss_tot)
          
          # Dynamic formatting to prevent rounding to 1.000
          r2_str = ">0.999" if r2 > 0.999 else f"{r2:.3f}"
          
          eq_text_return = fr'$y = {a:.2f}x^{{-{b:.2f}}} + {c:.2f} \ \ \ (R^2 = {r2_str})$'
          return_data = {'eq_text': eq_text_return, 'params': (a, b, c, np.min(x_fit), np.max(x_fit))}
          
          if extra_curves is not None:
              extra_colors = {'10m_Grid': 'tab:green', '30m_Grid': 'tab:orange'}
              extra_styles = {'10m_Grid': ':', '30m_Grid': '-.'}
              
              for scale_key, p_tuple in extra_curves.items():
                  if p_tuple is not None:
                      ea, eb, ec, ex_min, ex_max = p_tuple
                      ex_smooth = np.geomspace(ex_min, ex_max, 200)
                      ey_smooth = power_law(ex_smooth, ea, eb, ec)
                      
                      mask_zoom = ex_smooth >= 5
                      if np.any(mask_zoom):
                          y_min_zoom = min(y_min_zoom, np.min(ey_smooth[mask_zoom]))
                          y_max_zoom = max(y_max_zoom, np.max(ey_smooth[mask_zoom]))
                          
                      lbl = scale_key.replace('_', ' ')
                      l_extra += ax.plot(ex_smooth, ey_smooth, 
                                         linestyle=extra_styles.get(scale_key, '-'), 
                                         color=extra_colors.get(scale_key, 'gray'), 
                                         linewidth=2, label=f'Fit: {lbl}')

          if combined_fits is not None:
              eq_10 = combined_fits['10m_Grid'].get(title, "N/A")
              eq_30 = combined_fits['30m_Grid'].get(title, "N/A")
              
              display_text = (
                  r'$\bf{10m:}$ ' + eq_10 + '\n' +
                  r'$\bf{30m:}$ ' + eq_30 + '\n' +
                  r'$\bf{110m:}$ ' + eq_text_return
              )
          else:
              display_text = fr'$y = {a:.2f}x^{{-{b:.2f}}} + {c:.2f}$' + '\n' + fr'$R^2 = {r2_str}$'
              
          ax.text(0.95, 0.95, display_text, transform=ax.transAxes, fontsize=14,
                  verticalalignment='top', horizontalalignment='right',
                  bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
          
      except Exception as e:
          print(f"Warning: Could not fit curve for {title} - {e}")

  ax.set_xscale('log')
  
  if zoom_to_fit:
      filtered_x_vals = [x for x in x_vals if x >= 5]
      ax.set_xticks(filtered_x_vals)
      ax.set_xticklabels([str(x) for x in filtered_x_vals], rotation=45)
      if filtered_x_vals:
          ax.set_xlim(left=min(filtered_x_vals) * 0.9, right=max(filtered_x_vals) * 1.1)
  else:
      ax.set_xticks(x_vals)
      ax.set_xticklabels([str(x) for x in x_vals], rotation=45)
      
  ax.grid(True, which="major", ls="--", alpha=0.5)

  if zoom_to_fit and y_min_zoom != np.inf and y_max_zoom != -np.inf:
      y_margin = (y_max_zoom - y_min_zoom) * 0.1
      if y_margin == 0:
          y_margin = y_max_zoom * 0.1
      ax.set_ylim(max(0, y_min_zoom - y_margin), y_max_zoom + y_margin)

  # Right Axis: Relative Error
  l3 = []
  if plot_secondary and rel_errs is not None:
      ax2 = ax.twinx()
      ax2.set_ylabel('Mean Relative % Error', color=color_rel, fontsize=16)
      l3 = ax2.plot(x_vals, rel_errs, marker='o', linestyle='-', color=color_rel, alpha=0.4, label='Relative % Error')
      ax2.tick_params(axis='y', labelcolor=color_rel, labelsize=14)

  lines = l1 + l3 + l_extra + l2 
  labels = [l.get_label() for l in lines]
  
  return return_data, lines, labels


# --- HELPER FUNCTION FOR BGR SUBSET PLOTTING (5TH FIGURE) ---
def plot_subset_convergence(ax, x_vals, mae_dict, title, x_label, y_label):
    color_map = {'0-17.5% BGR': 'tab:blue', '17.5-35% BGR': 'tab:orange', '>35% BGR': 'tab:green'}
    style_map = {'0-17.5% BGR': '-', '17.5-35% BGR': '--', '>35% BGR': '-.'}
    
    ax.set_title(title, fontsize=18, fontweight='bold')
    ax.set_xlabel(x_label, fontsize=16)
    ax.set_ylabel(y_label, color='black', fontsize=16) 
    ax.tick_params(axis='y', labelsize=14, labelcolor='black')
    ax.tick_params(axis='x', labelsize=14)
    
    eq_texts, all_lines = [], []
    x_arr = np.array(x_vals, dtype=float)
    
    y_min_zoom = np.inf
    y_max_zoom = -np.inf
    
    for subset_name, mae_errs in mae_dict.items():
        c = color_map[subset_name]
        ls = style_map[subset_name]
        
        ax.plot(x_vals, mae_errs, marker='s', linestyle='', color=c, alpha=0.3)
        
        y_arr = np.array(mae_errs, dtype=float)
        valid_mask = ~np.isnan(y_arr) & ~np.isnan(x_arr)
        x_valid = x_arr[valid_mask]
        y_valid = y_arr[valid_mask]
        
        fit_mask = x_valid >= 5
        x_fit = x_valid[fit_mask]
        y_fit = y_valid[fit_mask]
        
        if len(y_fit) > 0:
            y_min_zoom = min(y_min_zoom, np.min(y_fit))
            y_max_zoom = max(y_max_zoom, np.max(y_fit))
        
        safe_name = subset_name.replace('%', r'\%')
        
        if len(x_fit) >= 3:
            try:
                p0 = [np.max(y_fit), 0.5, np.min(y_fit)]
                bounds = ([-np.inf, 0.01, -np.inf], [np.inf, 5.0, np.inf])
                params, _ = curve_fit(power_law, x_fit, y_fit, p0=p0, bounds=bounds, maxfev=10000)
                a, b, c_param = params
                
                x_smooth = np.geomspace(np.min(x_fit), np.max(x_fit), 200)
                y_smooth = power_law(x_smooth, a, b, c_param)
                
                l = ax.plot(x_smooth, y_smooth, linestyle=ls, color=c, linewidth=2, label=f'Fit: {subset_name}')
                all_lines.extend(l)
                
                y_pred = power_law(x_fit, a, b, c_param)
                ss_res = np.sum((y_fit - y_pred) ** 2)
                ss_tot = np.sum((y_fit - np.mean(y_fit)) ** 2)
                r2 = 1 - (ss_res / ss_tot)
                
                # Dynamic formatting to prevent rounding to 1.000
                r2_str = ">0.999" if r2 > 0.999 else f"{r2:.3f}"
                
                eq_texts.append(fr'$\bf{{{safe_name}:}}$ $y = {a:.2f}x^{{-{b:.2f}}} + {c_param:.2f} \ (R^2 = {r2_str})$')
            except Exception as e:
                eq_texts.append(fr'$\bf{{{safe_name}:}}$ Fit Failed')
        else:
            eq_texts.append(fr'$\bf{{{safe_name}:}}$ Not enough data')

    ax.set_xscale('log')
    
    filtered_x_vals = [x for x in x_vals if x >= 5]
    ax.set_xticks(filtered_x_vals)
    ax.set_xticklabels([str(x) for x in filtered_x_vals], rotation=45)
    
    if filtered_x_vals:
        ax.set_xlim(left=min(filtered_x_vals) * 0.9, right=max(filtered_x_vals) * 1.1)

    ax.grid(True, which="major", ls="--", alpha=0.5)
    
    if y_min_zoom != np.inf and y_max_zoom != -np.inf:
        y_margin = (y_max_zoom - y_min_zoom) * 0.1
        if y_margin == 0:
            y_margin = y_max_zoom * 0.1 
        ax.set_ylim(max(0, y_min_zoom - y_margin), y_max_zoom + y_margin)
    
    display_text = '\n'.join(eq_texts)
    ax.text(0.95, 0.95, display_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
            
    return all_lines


# --- HELPER TO BUILD MASTER LEGEND ---
def build_master_legend(legend_ax, all_lines, all_labels):
    handles_dict = {}
    for l, lbl in zip(all_lines, all_labels):
        if lbl not in handles_dict:
            handles_dict[lbl] = l
            
    legend_ax.axis('off')
    legend_ax.legend(handles_dict.values(), handles_dict.keys(), loc='center right', fontsize=16, frameon=True, borderpad=1.5, title="Metric Legend", title_fontsize=18)

# --- GENERATE INITIAL 3 FIGURES ---
gap_cats = {
    'Gap_0_24': 'Canopy Gap: 0-24 cm',
    'Gap_25_50': 'Canopy Gap: 25-50 cm',
    'Gap_51_100': 'Canopy Gap: 51-100 cm',
    'Gap_101_200': 'Canopy Gap: 101-200 cm',
    'Gap_gt_200': 'Canopy Gap: >200 cm'
}
gap_keys = list(gap_cats.keys())
gap_coords = [(2, 1), (3, 0), (3, 1), (4, 0), (4, 1)]

for scale_name, csv_filename in csv_files.items():
  csv_path = os.path.join(data_folder, csv_filename)
  
  if not os.path.exists(csv_path):
    print(f"Skipping {scale_name}: Could not find {csv_path}")
    continue
    
  print(f"\n--- Generating Plots for {scale_name} ---")
  df = pd.read_csv(csv_path)

  pt_cols = [c for c in df.columns if c.startswith('BGR_pt_')]
  if not pt_cols:
    continue
  points = sorted([int(c.split('_')[-1]) for c in pt_cols])

  ln_cols = [c for c in df.columns if c.startswith('Gap_0_24_L_')]
  lines = sorted([int(c.split('_')[-1]) for c in ln_cols])

  bgr_mae, bgr_rel = [], []
  herb_mae, herb_rel = [], []
  woody_mae, woody_rel = [], []
  fetch_mae, fetch_rel = [], []

  for pt in points:
    b_col = f'BGR_pt_{pt}'
    bgr_mae.append(np.abs(df[b_col] - df['BGR_Exact']).mean())
    bgr_rel.append((np.abs(df[b_col] - df['BGR_Exact']) / (df['BGR_Exact'] + 1e-6) * 100).mean())
    
    h_col = f'Herb_Pct_pt_{pt}'
    if h_col in df.columns and 'Herb_Pct_Exact' in df.columns:
      herb_mae.append(np.abs(df[h_col] - df['Herb_Pct_Exact']).mean())
      herb_rel.append((np.abs(df[h_col] - df['Herb_Pct_Exact']) / (df['Herb_Pct_Exact'] + 1e-6) * 100).mean())
    else:
      herb_mae.append(np.nan), herb_rel.append(np.nan)

    w_col = f'Woody_Pct_pt_{pt}'
    if w_col in df.columns and 'Woody_Pct_Exact' in df.columns:
      woody_mae.append(np.abs(df[w_col] - df['Woody_Pct_Exact']).mean())
      woody_rel.append((np.abs(df[w_col] - df['Woody_Pct_Exact']) / (df['Woody_Pct_Exact'] + 1e-6) * 100).mean())
    else:
      woody_mae.append(np.nan), woody_rel.append(np.nan)

    f_col = f'Fetch_pt_{pt}'
    fetch_mae.append(np.abs(df[f_col] - df['Fetch_Exact']).mean())
    fetch_rel.append((np.abs(df[f_col] - df['Fetch_Exact']) / (df['Fetch_Exact'] + 1e-6) * 100).mean())
  
  gap_data = {}
  for cat, title in gap_cats.items():
    g_mae, g_rel = [], []
    exact_col = f'{cat}_Exact'
    for ln in lines:
      l_col = f'{cat}_L_{ln}'
      g_mae.append(np.abs(df[l_col] - df[exact_col]).mean())
      g_rel.append((np.abs(df[l_col] - df[exact_col]) / (df[exact_col] + 1e-6) * 100).mean())
    gap_data[cat] = {'title': title, 'mae': g_mae, 'rel': g_rel}

  data_cache[scale_name] = {
      'points': points, 'lines': lines,
      'bgr_mae': bgr_mae, 'bgr_rel': bgr_rel,
      'herb_mae': herb_mae, 'herb_rel': herb_rel,
      'woody_mae': woody_mae, 'woody_rel': woody_rel,
      'fetch_mae': fetch_mae, 'fetch_rel': fetch_rel,
      'gap_data': gap_data
  }

  fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(18, 18))
  master_lines, master_labels = [], []

  res_bgr, l, lbl = plot_error_convergence(axes[0, 0], points, bgr_mae, bgr_rel, 'Bare Ground Percentage', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (pp)')
  fit_results[scale_name]['Bare Ground Percentage'] = res_bgr['eq_text']
  fit_params[scale_name]['Bare Ground Percentage'] = res_bgr['params']
  master_lines.extend(l); master_labels.extend(lbl)

  res_herb, l, lbl = plot_error_convergence(axes[1, 0], points, herb_mae, herb_rel, 'Herb Cover Percentage', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (pp)')
  fit_results[scale_name]['Herb Cover Percentage'] = res_herb['eq_text']
  fit_params[scale_name]['Herb Cover Percentage'] = res_herb['params']
  master_lines.extend(l); master_labels.extend(lbl)
  
  res_woody, l, lbl = plot_error_convergence(axes[1, 1], points, woody_mae, woody_rel, 'Woody Cover Percentage', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (pp)')
  fit_results[scale_name]['Woody Cover Percentage'] = res_woody['eq_text']
  fit_params[scale_name]['Woody Cover Percentage'] = res_woody['params']
  master_lines.extend(l); master_labels.extend(lbl)

  res_fetch, l, lbl = plot_error_convergence(axes[2, 0], points, fetch_mae, fetch_rel, 'Mean Fetch', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (m)')
  fit_results[scale_name]['Mean Fetch'] = res_fetch['eq_text']
  fit_params[scale_name]['Mean Fetch'] = res_fetch['params']
  master_lines.extend(l); master_labels.extend(lbl)
  
  for i, cat in enumerate(gap_keys):
      r, c = gap_coords[i]
      cat_title = gap_data[cat]['title']
      res_gap, l, lbl = plot_error_convergence(axes[r, c], lines, gap_data[cat]['mae'], gap_data[cat]['rel'], cat_title, 'Virtual Transect Length in Meters (Log Scale)', 'Absolute Error (pp)')
      fit_results[scale_name][cat_title] = res_gap['eq_text']
      fit_params[scale_name][cat_title] = res_gap['params']
      master_lines.extend(l); master_labels.extend(lbl)

  build_master_legend(axes[0, 1], master_lines, master_labels)

  fig.suptitle(f'Convergence of Sampled Metrics to Exact Values ({scale_name.replace("_", " ")})', fontsize=28, fontweight='bold', y=0.99)
  fig.tight_layout(pad=3.0)  

  output_img_path = os.path.join(out_folder, f'{scale_name}_Convergence.svg')
  plt.savefig(output_img_path, format='svg', dpi=300)
  plt.close() 
  print(f"Saved figure to: {output_img_path}")


# --- GENERATE 4TH SUMMARY FIGURE ---
if '110m_NRI' in data_cache:
    print(f"\n--- Generating 4th Cross-Scale Summary Plot (110m NRI Base) ---")
    c_data = data_cache['110m_NRI']
    
    fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(18, 18))
    master_lines, master_labels = [], []

    def get_extra_curves(metric_title):
        return {
            '10m_Grid': fit_params['10m_Grid'].get(metric_title),
            '30m_Grid': fit_params['30m_Grid'].get(metric_title)
        }

    summary_mae_label = 'Mean Absolute Error (110m NRI; 10m & 30m not shown)'

    _, l, lbl = plot_error_convergence(axes[0, 0], c_data['points'], c_data['bgr_mae'], None, 
                   'Bare Ground Percentage', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (pp)', 
                   combined_fits=fit_results, plot_secondary=False, extra_curves=get_extra_curves('Bare Ground Percentage'), mae_label=summary_mae_label, zoom_to_fit=True)
    master_lines.extend(l); master_labels.extend(lbl)
    
    _, l, lbl = plot_error_convergence(axes[1, 0], c_data['points'], c_data['herb_mae'], None,  
                   'Herb Cover Percentage', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (pp)', 
                   combined_fits=fit_results, plot_secondary=False, extra_curves=get_extra_curves('Herb Cover Percentage'), mae_label=summary_mae_label, zoom_to_fit=True)
    master_lines.extend(l); master_labels.extend(lbl)

    _, l, lbl = plot_error_convergence(axes[1, 1], c_data['points'], c_data['woody_mae'], None, 
                   'Woody Cover Percentage', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (pp)', 
                   combined_fits=fit_results, plot_secondary=False, extra_curves=get_extra_curves('Woody Cover Percentage'), mae_label=summary_mae_label, zoom_to_fit=True)
    master_lines.extend(l); master_labels.extend(lbl)
    
    _, l, lbl = plot_error_convergence(axes[2, 0], c_data['points'], c_data['fetch_mae'], None, 
                   'Mean Fetch', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (m)', 
                   combined_fits=fit_results, plot_secondary=False, extra_curves=get_extra_curves('Mean Fetch'), mae_label=summary_mae_label, zoom_to_fit=True)
    master_lines.extend(l); master_labels.extend(lbl)

    for i, cat in enumerate(gap_keys):
        r, c = gap_coords[i]
        cat_title = c_data['gap_data'][cat]['title']
        _, l, lbl = plot_error_convergence(axes[r, c], c_data['lines'], c_data['gap_data'][cat]['mae'], None, 
                     cat_title, 'Virtual Transect Length in Meters (Log Scale)', 'Absolute Error (pp)', 
                     combined_fits=fit_results, plot_secondary=False, extra_curves=get_extra_curves(cat_title), mae_label=summary_mae_label, zoom_to_fit=True)
        master_lines.extend(l); master_labels.extend(lbl)

    build_master_legend(axes[0, 1], master_lines, master_labels)

    fig.suptitle('Convergence Summary Across Spatial Scales (Mapped to 110m NRI Base)', fontsize=28, fontweight='bold', y=0.99)
    fig.tight_layout(pad=3.0)  

    summary_img_path = os.path.join(out_folder, '110m_NRI_CrossScale_Summary.svg')
    plt.savefig(summary_img_path, format='svg', dpi=300)
    plt.close() 
    print(f"Saved summary figure to: {summary_img_path}")


# --- GENERATE 5TH BGR SUBSET SUMMARY FIGURE ---
if '110m_NRI' in csv_files:
    print(f"\n--- Generating 5th BGR Subset Summary Plot (110m NRI Base) ---")
    
    df_110 = pd.read_csv(os.path.join(data_folder, csv_files['110m_NRI']))
    
    subset_masks = {
        '0-17.5% BGR': df_110['BGR_Exact'] <= 17.5,
        '17.5-35% BGR': (df_110['BGR_Exact'] > 17.5) & (df_110['BGR_Exact'] <= 35.0),
        '>35% BGR': df_110['BGR_Exact'] > 35.0
    }
    
    s_data = {s: {'bgr_mae': [], 'herb_mae': [], 'woody_mae': [], 'fetch_mae': [], 'gap_data': {cat: {'mae': []} for cat in gap_keys}} for s in subset_masks}
    
    c_points = data_cache['110m_NRI']['points']
    c_lines = data_cache['110m_NRI']['lines']
    
    for s_name, mask in subset_masks.items():
        sub_df = df_110[mask]
        
        for pt in c_points:
            s_data[s_name]['bgr_mae'].append(np.abs(sub_df[f'BGR_pt_{pt}'] - sub_df['BGR_Exact']).mean())
            
            h_col = f'Herb_Pct_pt_{pt}'
            if h_col in sub_df.columns and 'Herb_Pct_Exact' in sub_df.columns:
                s_data[s_name]['herb_mae'].append(np.abs(sub_df[h_col] - sub_df['Herb_Pct_Exact']).mean())
            else:
                s_data[s_name]['herb_mae'].append(np.nan)
                
            w_col = f'Woody_Pct_pt_{pt}'
            if w_col in sub_df.columns and 'Woody_Pct_Exact' in sub_df.columns:
                s_data[s_name]['woody_mae'].append(np.abs(sub_df[w_col] - sub_df['Woody_Pct_Exact']).mean())
            else:
                s_data[s_name]['woody_mae'].append(np.nan)
                
            f_col = f'Fetch_pt_{pt}'
            s_data[s_name]['fetch_mae'].append(np.abs(sub_df[f_col] - sub_df['Fetch_Exact']).mean())

        for cat in gap_keys:
            exact_col = f'{cat}_Exact'
            for ln in c_lines:
                l_col = f'{cat}_L_{ln}'
                s_data[s_name]['gap_data'][cat]['mae'].append(np.abs(sub_df[l_col] - sub_df[exact_col]).mean())

    fig_sub, axes_sub = plt.subplots(nrows=5, ncols=2, figsize=(18, 18))
    master_lines_sub = []
    
    def get_subset_dict(metric_type, is_gap=False, gap_cat=None):
        out = {}
        for s in subset_masks.keys():
            if is_gap:
                out[s] = s_data[s]['gap_data'][gap_cat]['mae']
            else:
                out[s] = s_data[s][metric_type]
        return out

    l = plot_subset_convergence(axes_sub[0, 0], c_points, get_subset_dict('bgr_mae'), 'Bare Ground Percentage', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (pp)')
    master_lines_sub.extend(l)
    
    l = plot_subset_convergence(axes_sub[1, 0], c_points, get_subset_dict('herb_mae'), 'Herb Cover Percentage', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (pp)')
    master_lines_sub.extend(l)
    
    l = plot_subset_convergence(axes_sub[1, 1], c_points, get_subset_dict('woody_mae'), 'Woody Cover Percentage', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (pp)')
    master_lines_sub.extend(l)
    
    l = plot_subset_convergence(axes_sub[2, 0], c_points, get_subset_dict('fetch_mae'), 'Mean Fetch', 'Number of Randomly Sampled Points (Log Scale)', 'Absolute Error (m)')
    master_lines_sub.extend(l)
    
    for i, cat in enumerate(gap_keys):
        r, c = gap_coords[i]
        cat_title = gap_cats[cat]
        l = plot_subset_convergence(axes_sub[r, c], c_lines, get_subset_dict('gap_data', is_gap=True, gap_cat=cat), cat_title, 'Virtual Transect Length in Meters (Log Scale)', 'Absolute Error (pp)')
        master_lines_sub.extend(l)

    master_labels_sub = [line.get_label() for line in master_lines_sub]
    build_master_legend(axes_sub[0, 1], master_lines_sub, master_labels_sub)
    
    fig_sub.suptitle('Convergence Summary Across BGR Subsets (110m NRI Base)', fontsize=28, fontweight='bold', y=0.99)
    fig_sub.tight_layout(pad=3.0)  

    subset_img_path = os.path.join(out_folder, '110m_NRI_Subset_Summary.svg')
    plt.savefig(subset_img_path, format='svg', dpi=300)
    plt.close() 
    print(f"Saved subset summary figure to: {subset_img_path}")

print("\nAll comprehensive plotting complete!")
