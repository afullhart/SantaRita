import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DEFAULT_OUTPUT_DIR = r"C:\Users\andre"

CSV_SCALE_MAP = {
    "5cm": "0cm",
    "25cm": "25cm",
    "50cm": "50cm",
    "100cm": "100cm",
    "200cm": "200cm",
}

def load_metrics(csv_path: str) -> pd.DataFrame:
    """Load aggregated simulation metrics and sort along the bare ground gradient."""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Simulation CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    if "True_BG_Mean" not in df.columns:
        raise ValueError("Expected column 'True_BG_Mean' in CSV.")

    return df.sort_values("True_BG_Mean").reset_index(drop=True)


def plot_metric_group(
    axes_row_or_col,
    mae_df: pd.DataFrame,
    prefix: str,
    scales: list[str],
    col_color: str | None,
    scale_colors: dict[str, str],
    exact_col_map: dict[str, str],
    is_vertical_stack: bool,
    show_x_vertical_50: bool = False,
    tick_labelsize: int | None = None,
) -> tuple[plt.Axes, plt.Axes, plt.Axes, plt.Axes]:
    """Plot Val, MAE, MRE, and Bias curves for one metric across sampling scales."""
    if is_vertical_stack:
        ax_val, ax_mae, ax_mre, ax_bias = axes_row_or_col[:, 0]
    else:
        ax_val, ax_mae, ax_mre, ax_bias = axes_row_or_col

    ax_val.plot(
        mae_df["True_BG_Mean"],
        mae_df[exact_col_map[prefix]],
        color="gray",
        linestyle="--",
        linewidth=2,
        label="Exact True Value",
        zorder=1,
    )

    for scale in scales:
        csv_scale = CSV_SCALE_MAP[scale]
        var_base = f"{prefix}_{csv_scale}"
        color = scale_colors[scale] if len(scales) > 1 else col_color
        label_name = f"{scale} Point"
        line_kws = {
            "marker": "o",
            "color": color,
            "linestyle": "-",
            "linewidth": 2.5,
            "markersize": 7,
        }

        ax_val.plot(mae_df["True_BG_Mean"], mae_df[f"{var_base}_Val"], label=label_name, **line_kws)
        ax_mae.plot(mae_df["True_BG_Mean"], mae_df[f"{var_base}_MAE"], **line_kws)
        ax_mre.plot(mae_df["True_BG_Mean"], mae_df[f"{var_base}_MRE"], **line_kws)
        ax_bias.plot(mae_df["True_BG_Mean"], mae_df[f"{var_base}_Bias"], **line_kws)

    for ax in (ax_val, ax_mae, ax_mre, ax_bias):
        ax.set_xlim(0, 100)
        ax.set_xticks(np.arange(0, 101, 20))
        ax.grid(True, alpha=0.3)
        if show_x_vertical_50:
            ax.axvline(50, color="gray", linestyle=":", linewidth=2)
        if tick_labelsize is not None:
            ax.tick_params(axis="both", labelsize=tick_labelsize)

    ax_bias.axhline(0, color="gray", linestyle="--", linewidth=1.5)
    return ax_val, ax_mae, ax_mre, ax_bias


def generate_plots(mae_df: pd.DataFrame, output_dir: str, show: bool = True) -> None:
    """Generate and save the multi-scale, cover, and spatial metric figures."""
    os.makedirs(output_dir, exist_ok=True)

    scale_colors = {
        "5cm": "black",
        "25cm": "forestgreen",
        "50cm": "dodgerblue",
        "100cm": "darkorange",
        "200cm": "crimson",
    }

    exact_col_map = {
        "BG": "True_BG_Mean",
        "Fetch": "Exact_Fetch_Mean",
        "Herb": "True_Herb_Mean",
        "Woody": "True_Woody_Mean",
        "Gap_0_24": "Exact_Gap_0_24_Mean",
        "Gap_25_50": "Exact_Gap_25_50_Mean",
        "Gap_51_100": "Exact_Gap_51_100_Mean",
        "Gap_101_200": "Exact_Gap_101_200_Mean",
        "Gap_gt_200": "Exact_Gap_gt_200_Mean",
    }

    # ====================================================================
    # FIGURE 1: Original 4x9 Configuration
    # ====================================================================
    plot_config = [
        ("BG", "Total Bare Ground (%)", ["5cm", "25cm", "50cm", "100cm", "200cm"], None),
        ("Fetch", "Mean Fetch (m)", ["5cm", "25cm", "50cm", "100cm", "200cm"], None),
        ("Gap_0_24", "Canopy Gap 0-24cm (%)", ["5cm"], "steelblue"),
        ("Gap_25_50", "Canopy Gap 25-50cm (%)", ["5cm"], "cadetblue"),
        ("Gap_51_100", "Canopy Gap 51-100cm (%)", ["5cm"], "mediumseagreen"),
        ("Gap_101_200", "Canopy Gap 101-200cm (%)", ["5cm"], "darkorange"),
        ("Gap_gt_200", "Canopy Gap >200cm (%)", ["5cm"], "firebrick"),
        ("Herb", "Herbaceous Cover (%)", ["5cm", "25cm", "50cm", "100cm", "200cm"], None),
        ("Woody", "Woody Cover (%)", ["5cm", "25cm", "50cm", "100cm", "200cm"], None),
    ]

    fig1, axes1 = plt.subplots(4, 9, figsize=(36, 16), constrained_layout=True)
    fig1.suptitle(
        "Simulation Metrics (Values, MAE, MRE, Bias) Across Bare Ground Gradient",
        fontsize=20,
        weight="bold",
    )

    for col, (prefix, title, scales, col_color) in enumerate(plot_config):
        ax_val, ax_mae, ax_mre, ax_bias = plot_metric_group(
            axes1[:, col],
            mae_df,
            prefix,
            scales,
            col_color,
            scale_colors,
            exact_col_map,
            is_vertical_stack=False,
        )
        ax_val.set_title(title, fontsize=15, pad=12)
        ax_val.legend(loc="upper left", fontsize=10)
        ax_bias.set_xlabel("True Bare Ground (%)", fontsize=12)

        if col == 0:
            ax_val.set_ylabel("Sampled Value", fontsize=13)
            ax_mae.set_ylabel("Mean Absolute Error", fontsize=13)
            ax_mre.set_ylabel("Mean Relative Error (%)", fontsize=13)
            ax_bias.set_ylabel("Mean Absolute Bias", fontsize=13)

    # ====================================================================
    # FIGURE 2: COVER METRICS (BG, HERB, WOODY) - 3 Rows x 4 Cols
    # ====================================================================
    fig2_config = [
        ("BG", "Total Bare Ground (%)", ["5cm", "25cm", "50cm", "100cm", "200cm"], None),
        ("Herb", "Herbaceous Cover (%)", ["5cm", "25cm", "50cm", "100cm", "200cm"], None),
        ("Woody", "Woody Cover (%)", ["5cm", "25cm", "50cm", "100cm", "200cm"], None),
    ]

    fig2, axes2 = plt.subplots(3, 4, figsize=(16, 12), constrained_layout=True)
    fig2.suptitle(
        "Simulation Cover Metrics Across Bare Ground Gradient",
        fontsize=20,
        weight="bold",
    )

    for row, (prefix, title, scales, col_color) in enumerate(fig2_config):
        ax_val, ax_mae, ax_mre, ax_bias = plot_metric_group(
            axes2[row, :],
            mae_df,
            prefix,
            scales,
            col_color,
            scale_colors,
            exact_col_map,
            is_vertical_stack=False,
        )
        ax_val.set_ylabel(title, fontsize=13, weight="bold")
        if row == 0:
            ax_val.legend(loc="upper left", fontsize=10)
            ax_val.set_title("Sampled Value", fontsize=15, pad=12)
            ax_mae.set_title("Mean Absolute Error", fontsize=15, pad=12)
            ax_mre.set_title("Mean Relative Error (%)", fontsize=15, pad=12)
            ax_bias.set_title("Mean Absolute Bias", fontsize=15, pad=12)

        if row == len(fig2_config) - 1:
            for ax in (ax_val, ax_mae, ax_mre, ax_bias):
                ax.set_xlabel("True Bare Ground (%)", fontsize=12)

    # ====================================================================
    # FIGURE 3: SPATIAL METRICS (FETCH, GAPS) - 6 Rows x 4 Cols
    # ====================================================================
    fig3_config = [
        ("Fetch", "Mean Fetch (m)", ["5cm", "25cm", "50cm", "100cm", "200cm"], None),
        ("Gap_0_24", "Canopy Gap\n0-24cm (%)", ["5cm"], "steelblue"),
        ("Gap_25_50", "Canopy Gap\n25-50cm (%)", ["5cm"], "cadetblue"),
        ("Gap_51_100", "Canopy Gap\n51-100cm (%)", ["5cm"], "mediumseagreen"),
        ("Gap_101_200", "Canopy Gap\n101-200cm (%)", ["5cm"], "darkorange"),
        ("Gap_gt_200", "Canopy Gap\n>200cm (%)", ["5cm"], "firebrick"),
    ]

    fig3, axes3 = plt.subplots(6, 4, figsize=(16, 24), constrained_layout=True)
    fig3.suptitle(
        "Simulation Spatial Structure Metrics Across Bare Ground Gradient",
        fontsize=24,
        weight="bold",
    )

    for row, (prefix, title, scales, col_color) in enumerate(fig3_config):
        ax_val, ax_mae, ax_mre, ax_bias = plot_metric_group(
            axes3[row, :],
            mae_df,
            prefix,
            scales,
            col_color,
            scale_colors,
            exact_col_map,
            is_vertical_stack=False,
            show_x_vertical_50=True,
            tick_labelsize=12,
        )
        ax_val.set_ylabel(title, fontsize=16, weight="bold")

        last_var_base = f"{prefix}_{CSV_SCALE_MAP[scales[-1]]}"
        start_val = mae_df[f"{last_var_base}_Val"].iloc[0]
        end_val = mae_df[f"{last_var_base}_Val"].iloc[-1]
        leg_loc = "upper left" if end_val >= start_val else "upper right"
        ax_val.legend(loc=leg_loc, fontsize=12)

        if row == len(fig3_config) - 1:
            for ax in (ax_val, ax_mae, ax_mre, ax_bias):
                ax.set_xlabel("True Bare Ground (%)", fontsize=15)

    # ====================================================================
    # SAVE AND DISPLAY PLOTS
    # ====================================================================
    img_path1 = os.path.join(output_dir, "Simulation_Metrics_MultiScale.svg")
    img_path2 = os.path.join(output_dir, "Simulation_Cover_Metrics.svg")
    img_path3 = os.path.join(output_dir, "Simulation_Spatial_Metrics_Rows.svg")

    fig1.savefig(img_path1, format="svg", dpi=300, bbox_inches="tight")
    fig2.savefig(img_path2, format="svg", dpi=300, bbox_inches="tight")
    fig3.savefig(img_path3, format="svg", dpi=300, bbox_inches="tight")

    print(f"Saved:\n  - {img_path1}\n  - {img_path2}\n  - {img_path3}")

    if show:
        plt.show()
    else:
        plt.close(fig1)
        plt.close(fig2)
        plt.close(fig3)


def parse_args() -> argparse.Namespace:
    default_csv = os.path.join(DEFAULT_OUTPUT_DIR, "Simulation_Metrics_MultiScale.csv")

    parser = argparse.ArgumentParser(
        description="Plot simulation metrics from Simulation_Metrics_MultiScale.csv."
    )
    parser.add_argument(
        "--csv",
        default=default_csv,
        help=f"Path to Simulation_Metrics_MultiScale.csv (default: {default_csv}).",
    )
    parser.add_argument(
        "--out-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory to save SVG plots (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Save figures without calling plt.show().",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    metrics_df = load_metrics(args.csv)
    generate_plots(metrics_df, output_dir=args.out_dir, show=not args.no_show)
