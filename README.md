
```mermaid
graph LR
    classDef input fill:#f9f9f9,stroke:#333,stroke-width:2px;
    classDef script fill:#e1f5fe,stroke:#3182bd,stroke-width:2px;

    %% Data Sources
    DEM_Raw[(USGS 3DEP 1m)]:::input
    S2_Raw[(Sentinel-2 SR)]:::input
    UAS_Raw[(Classified sUAS Orthos)]:::input

    %% Processing Scripts
    DEM_Export[Export_10m_DEM.js]:::script
    Cloud_Export[Export_Cloud_FeatureClass.js]:::script
    FC_Model[Model_FeatureClass.js]:::script
    Reg_Model[Model_Regressions.js]:::script

    %% Visualization & Analysis
    Vis_Model[Model_Visualization.js]:::script
    Comp_RAP[RAP_Comparison.js]:::script
    Vis_RAP[RAP_Visualization.js]:::script

    %% Flow
    DEM_Raw --> DEM_Export
    S2_Raw --> Cloud_Export

    UAS_Raw --> FC_Model
    DEM_Export -->|10m DEM Asset| FC_Model
    Cloud_Export -->|Optimal Dates| FC_Model

    FC_Model -->|Training Data| Reg_Model

    Reg_Model -->|Hyperparameters| Vis_Model
    Reg_Model -->|Validation| Comp_RAP
    Reg_Model -->|Hyperparameters| Vis_RAP

    %% Asset Dependencies
    DEM_Export -.-> Vis_Model
    DEM_Export -.-> Vis_RAP
    Cloud_Export -.-> Vis_Model
    Cloud_Export -.-> Vis_RAP
```

### Scripts
* **`Export_10m_DEM.js`**: Resamples the high-resolution 1m USGS 3DEP Digital Elevation Model to a 10m spatial resolution using a true spatial mean reduction to avoid aliasing errors, and exports it as a stripe-free asset.
* **`Export_Cloud_FeatureClass.js`**: Scans the 2018-2025 Sentinel-2 Harmonized record to locate the optimal 7-day imagery windows for each month based on a custom "True Obscuration" metric, exporting a feature class asset of the clearest dates.
* **`Model_FeatureClass.js`**: Generates 10m Sentinel-2 grid feature class asset that is filtered by spatial overlap with footprints. Adds BGR, LPI, and MFT. Adds Sentinel-2 band values. Adds NDVI, MCARI, BSI, and NBR2. Adds slope from 10m USGS DEM.
* **`Model_Regressions.js`**: Executes K-folds cross-validation and final models, visualizes prediction map, and exports predicted vs. true values.
* **`Model_Visualization.js`**: Loads the tuned Gradient Tree Boost models, resampled DEM, and optimal cloud windows to dynamically generate interactive, topography-aware predictive maps and local time-series charts.
* **`RAP_Comparison.js`**: Isolates drone-based ground truth data and compares it against the 10m Rangeland Analysis Platform (RAP) by computing Cumulative Distribution Functions (CDFs) to visualize fractional cover distributions.
* **`RAP_Visualization.js`**: Maps the spatial difference (residuals) between the custom SRER predictive model and the 10m RAP product, aggregating monthly predictions into an annual mean and providing an interactive pixel inspector for historical trends.
