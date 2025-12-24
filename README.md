# GCM–RCM (Central Europe hydroclimate)

Scripts for downloading, preprocessing, and analyzing CMIP5/CMIP6 global climate model (GCM) simulations and EURO-CORDEX EUR-44 regional climate model (RCM) simulations for the Central Europe hydroclimate study.

## About
This repository contains the analysis workflow used in the manuscript:

> *A dialogue concerning hydroclimatic projections by global and regional climate models in Central Europe*

The focus is on hydroclimatic diagnostics (e.g., evapotranspiration, vapor pressure deficit, Budyko framework, effective conductance diagnostics) and on characterizing ensemble behavior across GCM and RCM simulations.

## Repository structure
- `downloading_and_preprocessing/`  
  Scripts for data retrieval and preprocessing (download, variable preparation, temporal aggregation, domain averaging, etc.).
- `data_analysis/`  
  Scripts for analysis, diagnostics, and figure generation.

## Data sources
Model output is retrieved from:
- Earth System Grid Federation (ESGF) for CMIP5/CMIP6
- CEDA for EURO-CORDEX
- Copernicus Climate Data Store (CDS) for reanalysis products (e.g., ERA5-Land)

**Note:** Input datasets are not distributed with this repository due to data volume and licensing/access constraints.

## Requirements
### Core software
- **R** (analysis and plotting)
- **CDO** (Climate Data Operators) for preprocessing and aggregation
- (Optional) **Python** for auxiliary utilities (if used in your local workflow)

### R packages (typical)
- `tidyverse`
- `ranger`
- `DEoptim`
- `MASS`
- (plus any additional packages used by scripts in `data_analysis/`)

## Typical workflow
1. **Download and preprocess**
   - Use scripts in `downloading_and_preprocessing/` to retrieve model data and derive required variables.
   - Aggregate daily data to monthly products if needed.
   - Create regional means over the Central European domain.

2. **Run diagnostics and generate figures**
   - Use scripts in `data_analysis/` to compute hydroclimatic metrics and diagnostics.
   - Generate manuscript and Extended Data figures.

## Configuration
Paths to local data directories and credentials (e.g., ESGF) are environment-specific.

Recommended approach:
- Store paths as environment variables (e.g., `DATA_DIR`, `OUTPUT_DIR`)
- Or maintain a local, untracked config file (e.g., `config_local.R`) and add it to `.gitignore`

## Outputs
Depending on the scripts executed, outputs may include:
- Cleaned/derived regional time series (monthly climatologies)
- Diagnostic metrics (e.g., Budyko-space quantities)
- Figures for the manuscript and Extended Data

## Reproducibility notes
- Analyses are conducted on **regional means** derived from spatially averaged model output.
- Time periods typically follow:
  - historical: 1981–2005
  - future: 2076–2100 (RCP8.5 for CMIP5; SSP5-8.5 for CMIP6)
- EUR-44 native horizontal resolution is ~0.44° (≈50 km).

## Citation
If you use code from this repository, please cite:
- the associated manuscript (once published / preprint link), and
- this GitHub repository (release tag preferred).

A suggested citation format will be added after publication.

## License
Add a license file if you intend others to reuse this code (e.g., MIT, BSD-3, GPL-3).  
If no license is provided, reuse is legally restricted by default.

## Contact
Open an issue on GitHub for questions, or contact the repository owner via the email listed in the paper.

