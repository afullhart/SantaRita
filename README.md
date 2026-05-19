```mermaid
---
title: GEE Workflow
config:
  themeCSS: ".title { font-size: 36px !important; font-weight: bold; }"
---
graph TD
  %% Styling Definitions
  classDef input fill:#e1f5fe,stroke:#039be5,stroke-width:2px,stroke-dasharray: 3 3;
  classDef intermediate fill:#fff3e0,stroke:#ffb74d,stroke-width:2px;
  classDef script fill:#e8f5e9,stroke:#4caf50,stroke-width:2px;
  classDef output fill:#ede7f6,stroke:#7e57c2,stroke-width:2px;
  classDef junction fill:#333,stroke:#333,circle;

  %% --- LEVEL 1: INPUT DATA & RAW TEXTURES ---
  subgraph Inputs [Data Sources & Footprints]
    I1["5cm Classified Images"]
    I2["Drone Image Footprints"]
    I3["Ecological Survey Footprints"]
    I4["10m Sentinel-2 Images"]
    I5["SRER Footprint"]
    
    %% Script-based pre-processing inputs
    S1["Export_10m_DEM.js"]
    S2["Export_Cloud_FeatureClass.js"]
  end

  %% --- LEVEL 2: PRE-PROCESSING & INTERMEDIATE PRODUCTS ---
  subgraph PreProcessing [Pre-Processing & Derived Features]
    P1["10m USGS Resampled DEM"]
    P2["Cloud-Free Feature Class"]
  end

  %% --- JOINTS ---
  J3(( )):::junction
  J1(( )):::junction
  J2(( )):::junction

  %% --- LEVEL 3: FEATURE ENGINEERING & TRAINING ---
  subgraph CoreModel [Model Feature Engineering & Training]
    S3["Model_FeatureClass.js"]
    T1["Training Feature Class"]
  end

  %% --- LEVEL 4: ANALYSIS & VISUALIZATION SCRIPTS ---
  subgraph Evaluation [Analysis, Regression & Visualization]
    O1["Model_Regressions.js"]
    O2["Model_Visualization.js"]
    O3["RAP_Comparison.js"]
    O4["RAP_Visualization.js"]
  end

  %% --- FLOW CONNECTIONS ---
  %% Left-side footprint pipeline gathered to a joint to narrow the subgraph
  I1 --> J3
  I2 --> J3
  I3 --> J3
  J3 --> S3

  %% Pre-processing script connections
  S1 --> P1
  S2 --> P2

  %% --- INPUTS TO THE FIRST JOINT ---
  I4 --> J1
  P1 --> J1
  P2 --> J1
  I5 --> J1

  %% JOINT 1 TO FEATURE CLASS SCRIPT
  J1 --> S3

  %% Feature generation to training data
  S3 --> T1

  %% --- Top Joint directly to Bottom Joint ---
  J1 --> J2

  %% --- TRAINING DATA TO SECOND JOINT ---
  T1 --> J2

  %% --- JOINT 2 TO DOWNSTREAM SCRIPTS ---
  J2 --> O1
  J2 --> O2
  J2 --> O3
  J2 --> O4

  %% Assigning Classes for Visual Distinction
  class I1,I2,I3,I4,I5 input;
  class P1,P2,T1 intermediate;
  class S1,S2,S3 script;
  class O1,O2,O3,O4 output;
```

### Scripts
* **`Export_10m_DEM.js`**: Resamples the high-resolution 1m USGS 3DEP Digital Elevation Model to a 10m spatial resolution using a true spatial mean reduction to avoid aliasing errors, and exports it as a stripe-free asset.
* **`Export_Cloud_FeatureClass.js`**: Scans the 2018-2025 Sentinel-2 Harmonized record to locate the optimal 7-day imagery windows for each month based on a custom "True Obscuration" metric, exporting a feature class asset of the clearest dates.
* **`Model_FeatureClass.js`**: Generates 10m Sentinel-2 grid feature class asset that is filtered by spatial overlap with footprints. Adds BGR, LPI, and MFT. Adds Sentinel-2 band values. Adds NDVI, MCARI, BSI, and NBR2. Adds slope, aspect, and illumination based on resampled 10m USGS DEM.
* **`Model_Regressions.js`**: Executes K-folds cross-validation and final models, visualizes prediction map, and exports predicted vs. true values.
* **`Model_Visualization.js`**: Loads the tuned Gradient Tree Boost models, resampled DEM, and optimal cloud windows to dynamically generate interactive, topography-aware predictive maps and local time-series charts.
* **`RAP_Comparison.js`**: Isolates drone-based ground truth data and compares it against the 10m Rangeland Analysis Platform (RAP) by computing Cumulative Distribution Functions (CDFs) to visualize fractional cover distributions.
* **`RAP_Visualization.js`**: Maps the spatial difference (residuals) between the custom SRER predictive model and the 10m RAP product, aggregating monthly predictions into an annual mean and providing an interactive pixel inspector for historical trends.

### Data sources
* **`5cm Classified Images`**: https://gee-community-catalog.org/projects/srer_drone/
* **`Drone Image Footprints`**: https://www.arcgis.com/home/item.html?id=50b30d505bd2491e9412217139b7df83
* **`Ecological Survey Footprints`**: https://www.arcgis.com/home/item.html?id=50b30d505bd2491e9412217139b7df83
* **`10m Sentinel-2 Images`**: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR_HARMONIZED
* **`Export_10m_DEM.js`**: https://developers.google.com/earth-engine/datasets/catalog/USGS_3DEP_1m
* **`Export_Cloud_FeatureClass.js`**: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR_HARMONIZED
* **`SRER Footprint`**: https://www.arcgis.com/home/item.html?id=1fa8b9b97e844aaf8170a95c6e2a3b76
