import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Load the exported CSV data from the Landsat workflow. 
# Added "r" before the string to safely handle Windows file paths.
df = pd.read_csv(r'G:\My Drive\GEE_Downloads\SRER_Regional_Average_TimeSeries_Landsat.csv')

# Ensure Date column is recognized as datetime objects and sorted
df['Date'] = pd.to_datetime(df['Date'])
df = df.sort_values('Date')

print("Data loaded successfully. Generating 40-year charts...")

# Adjust marker size and line width for higher data density
msize = 3
lwidth = 1.5

# =========================================================
# CHART 1: Core Metrics (BGR, LPI, Fetch)
# =========================================================
fig1, ax1 = plt.subplots(figsize=(12, 6))

# Plot Percentages on the primary (left) y-axis
l1, = ax1.plot(df['Date'], df['Mean_BGR_pct'], marker='.', markersize=msize, color='#d73027', linewidth=lwidth, label='Mean BGR (%)')
l2, = ax1.plot(df['Date'], df['Mean_LPI_pct'], marker='.', markersize=msize, color='#fc8d59', linewidth=lwidth, label='Mean LPI (%)')
ax1.set_xlabel('Date')
ax1.set_ylabel('Percentage (%)', color='black')
ax1.set_ylim(0, 60)
ax1.grid(True, linestyle='--', alpha=0.5)

# Plot Fetch on the secondary (right) y-axis
ax2 = ax1.twinx()
l3, = ax2.plot(df['Date'], df['Mean_Fetch_m'], marker='.', markersize=msize, color='#1a9850', linewidth=lwidth, label='Mean Fetch (m)')
ax2.set_ylabel('Mean Fetch Distance (m)', color='black')
ax2.set_ylim(0, 0.30)

# Configure Title and Unified Legend
plt.title('SRER Regional Average: Predicted Core Metrics (30m Landsat Resolution)', pad=45, fontsize=14, fontweight='bold')
lines = [l1, l2, l3]
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper left', bbox_to_anchor=(0, 1.10), ncol=3, frameon=False)

# Format x-axis for 40 years (Tick every 5 years, display Year only)
ax1.xaxis.set_major_locator(mdates.YearLocator(5))
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')

plt.tight_layout()
fig1.savefig('core_metrics_30m.png', dpi=300) 

# =========================================================
# CHART 2: Log Herb-to-Woody Ratio (HWR)
# =========================================================
fig2, ax3 = plt.subplots(figsize=(12, 5))

# Plot Log HWR on a single axis
ax3.plot(df['Date'], df['Mean_Log_HWR'], marker='.', markersize=msize, color='#4575b4', linewidth=lwidth, label='Mean Log HWR')
ax3.set_xlabel('Date')
ax3.set_ylabel('Mean Log Ratio', color='black')
ax3.set_ylim(0, 3.5)
ax3.grid(True, linestyle='--', alpha=0.5)

# Configure Title and Legend
plt.title('SRER Regional Average: Predicted Herb-to-Woody Ratio (30m Landsat Resolution)', pad=45, fontsize=14, fontweight='bold')
ax3.legend(loc='upper left', bbox_to_anchor=(0, 1.10), frameon=False)

# Format x-axis dates
ax3.xaxis.set_major_locator(mdates.YearLocator(5))
ax3.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.setp(ax3.xaxis.get_majorticklabels(), rotation=45, ha='right')

plt.tight_layout()
fig2.savefig('log_hwr_30m.png', dpi=300) 

# =========================================================
# CHART 3: Predicted Herbaceous & Woody Cover
# =========================================================
fig3, ax4 = plt.subplots(figsize=(12, 5))

# Plot actual predicted Herb and Woody columns
ax4.plot(df['Date'], df['Mean_Herb_pct'], marker='.', markersize=msize, color='#91cf60', linewidth=lwidth, label='Mean Herbaceous (%)')
ax4.plot(df['Date'], df['Mean_Woody_pct'], marker='.', markersize=msize, color='#8c510a', linewidth=lwidth, label='Mean Woody (%)')

ax4.set_xlabel('Date')
ax4.set_ylabel('Cover Percentage (%)', color='black')
ax4.set_ylim(0, 100) # Lock to 0-100% scale for easy visual comparison
ax4.grid(True, linestyle='--', alpha=0.5)

# Configure Title and Legend
plt.title('SRER Regional Average: Predicted Absolute Cover (30m Landsat Resolution)', pad=45, fontsize=14, fontweight='bold')
ax4.legend(loc='upper left', bbox_to_anchor=(0, 1.10), ncol=2, frameon=False)

# Format x-axis dates
ax4.xaxis.set_major_locator(mdates.YearLocator(5))
ax4.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.setp(ax4.xaxis.get_majorticklabels(), rotation=45, ha='right')

plt.tight_layout()
fig3.savefig('herb_woody_cover_30m.png', dpi=300) 

# Show all charts
plt.show()

print("Charts successfully saved to your directory!")
