import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Load the exported CSV data. This is the output from Model_RegionalTrends.js
df = pd.read_csv('G:\My Drive\GEE_Downloads\SRER_Regional_Average_TimeSeries_10m.csv')

# Ensure Date column is recognized as datetime objects and sorted
df['Date'] = pd.to_datetime(df['Date'])
df = df.sort_values('Date')

print("Data loaded successfully. Generating charts...")

# =========================================================
# CHART 1: Core Metrics (BGR, LPI, Fetch)
# =========================================================
fig1, ax1 = plt.subplots(figsize=(12, 6))

# Plot Percentages on the primary (left) y-axis
l1, = ax1.plot(df['Date'], df['Mean_BGR_pct'], marker='.', color='#d73027', linewidth=2, label='Mean BGR (%)')
l2, = ax1.plot(df['Date'], df['Mean_LPI_pct'], marker='.', color='#fc8d59', linewidth=2, label='Mean LPI (%)')
ax1.set_xlabel('Date')
ax1.set_ylabel('Percentage (%)', color='black')
ax1.set_ylim(0, 60)
ax1.grid(True, linestyle='--', alpha=0.5)

# Plot Fetch on the secondary (right) y-axis
ax2 = ax1.twinx()
l3, = ax2.plot(df['Date'], df['Mean_Fetch_m'], marker='.', color='#1a9850', linewidth=2, label='Mean Fetch (m)')
ax2.set_ylabel('Mean Fetch Distance (m)', color='black')
ax2.set_ylim(0, 0.30)

# Configure Title and Unified Legend
# FIX: Increased pad to 45 to make room for the legend underneath
plt.title('SRER Regional Average: Predicted Core Metrics (10m Native Resolution)', pad=45, fontsize=14, fontweight='bold')
lines = [l1, l2, l3]
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper left', bbox_to_anchor=(0, 1.10), ncol=3, frameon=False)

# Format x-axis dates beautifully
ax1.xaxis.set_major_locator(mdates.AutoDateLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')

plt.tight_layout()
fig1.savefig('core_metrics_10m.png', dpi=300) # Saves high-res image
plt.show()

# =========================================================
# CHART 2: Log Herb-to-Woody Ratio (HWR)
# =========================================================
fig2, ax3 = plt.subplots(figsize=(12, 5))

# Plot Log HWR on a single axis
ax3.plot(df['Date'], df['Mean_Log_HWR'], marker='.', color='#4575b4', linewidth=2, label='Mean Log HWR')
ax3.set_xlabel('Date')
ax3.set_ylabel('Mean Log Ratio', color='black')
ax3.set_ylim(0, 3.5)
ax3.grid(True, linestyle='--', alpha=0.5)

# Configure Title and Legend
# FIX: Increased pad to 45 to make room for the legend underneath
plt.title('SRER Regional Average: Predicted Herb-to-Woody Ratio (10m Native Resolution)', pad=45, fontsize=14, fontweight='bold')
ax3.legend(loc='upper left', bbox_to_anchor=(0, 1.10), frameon=False)

# Format x-axis dates beautifully
ax3.xaxis.set_major_locator(mdates.AutoDateLocator())
ax3.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
plt.setp(ax3.xaxis.get_majorticklabels(), rotation=45, ha='right')

plt.tight_layout()
fig2.savefig('log_hwr_10m.png', dpi=300) # Saves high-res image
plt.show()

print("Charts successfully saved to your directory!")
