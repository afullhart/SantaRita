```mermaid
---
config:
  layout: elk
  theme: default
  look: handDrawn
  edges: straight
---
flowchart TD
    classDef input fill:#f9f9f9,stroke:#333,stroke-width:2px;
    classDef script fill:#e1f5fe,stroke:#3182bd,stroke-width:2px;
    classDef asset fill:#fff3e0,stroke:#ff9800,stroke-width:2px;

    %% Data Sources (Top Row Inputs)
    A[5cm classified images]:::input
    B[Drone image Footprints]:::input
    C[Ecological Survey Footprints]:::input
    D[10m Sentinel-2 images]:::input
    E[Export_10m_DEM.js]:::script
    F[Export_Cloud_FeatureClass.js]:::script
    G[SRER Footprint]:::input

    %% Intermediate Assets (The First Manifold Stage)
    H[10m USGS Resampled DEM]:::asset
    I[Cloud-Free Feature Class]:::asset

    %% Feature Extraction Core
    J[Model_FeatureClass.js]:::script
    K[Training Feature Class]:::asset

    %% Analysis & Execution
    L[Model_Regressions.js]:::script

    %% Downstream Visualization & Comparison (Bottom Row Outputs)
    M[Model_Visualization.js]:::script
    N[RAP_Comparison.js]:::script
    O[RAP_Visualization.js]:::script

    %% --- CONNECTIONS & FLOW ---

    %% Pre-Processing Flows
    D --> D1--> H
    E --> E1--> H
    F --> F1--> I

    %% Central Extraction Manifold (5 Data Paths Merging Into Feature Class)
    A --> J
    B --> J
    C --> J
    H --> J
    I --> J

    %% Post-Extraction Path
    J --> K
    K --> L

    %% The Linear Regression Multi-Funnels
    L --> M
    L --> N
    L --> O

    %% Global Asset Bypasses
    H -.-> M
    H -.-> O
    I -.-> M
    I -.-> O
    G --> O
```

### Scripts
* **`Export_10m_DEM.js`**: Resamples the high-resolution 1m USGS 3DEP Digital Elevation Model to a 10m spatial resolution using a true spatial mean reduction to avoid aliasing errors, and exports it as a stripe-free asset.
* **`Export_Cloud_FeatureClass.js`**: Scans the 2018-2025 Sentinel-2 Harmonized record to locate the optimal 7-day imagery windows for each month based on a custom "True Obscuration" metric, exporting a feature class asset of the clearest dates.
* **`Model_FeatureClass.js`**: Generates 10m Sentinel-2 grid feature class asset that is filtered by spatial overlap with footprints. Adds BGR, LPI, and MFT. Adds Sentinel-2 band values. Adds NDVI, MCARI, BSI, and NBR2. Adds slope from 10m USGS DEM.
* **`Model_Regressions.js`**: Executes K-folds cross-validation and final models, visualizes prediction map, and exports predicted vs. true values.
* **`Model_Visualization.js`**: Loads the tuned Gradient Tree Boost models, resampled DEM, and optimal cloud windows to dynamically generate interactive, topography-aware predictive maps and local time-series charts.
* **`RAP_Comparison.js`**: Isolates drone-based ground truth data and compares it against the 10m Rangeland Analysis Platform (RAP) by computing Cumulative Distribution Functions (CDFs) to visualize fractional cover distributions.
* **`RAP_Visualization.js`**: Maps the spatial difference (residuals) between the custom SRER predictive model and the 10m RAP product, aggregating monthly predictions into an annual mean and providing an interactive pixel inspector for historical trends.
