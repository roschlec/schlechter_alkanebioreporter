# The alkane 1-monooxigenase gene *alkB* of *Pseudomonas* sp. FF2 is upregulated during colonisation of *Arabidopsis thaliana* leaves

## Overview

This repository contains R scripts for analyzing bacterial single-cell fluorescence and colony-forming unit (CFU) data in both *in vitro* and *in planta* conditions. The analysis includes data preprocessing, classification models, and statistical evaluation of bacterial fitness and fluorescence. The associated microscopy and raw data (CFU and single-cell) can be found in the Zenodo repository.

## Repository Structure

```
📂 data/                  # Input datasets (not included in repo)
📂 results/               # Output files including model results and plots (must be created)
📂 src/                   # Additional helper scripts
  📄 01_invitro_cell_classifier.R   # Machine learning classifier for *in vitro* single cells
  📄 01_invitro_cfu.R                # CFU analysis for *in vitro* samples
  📄 01_invitro_microscopy.R         # *In vitro* microscopy and fluorescence analysis
  📄 02_inplanta_cell_classifier.R   # Machine learning classifier for *in planta* single cells
  📄 02_inplanta_cfu.R               # CFU analysis for *in planta* samples
  📄 02_inplanta_microscopy.R        # *In planta* microscopy and fluorescence analysis
  📄 macros/                # ImageJ macros for cell and background selection
```

## Installation

To run the scripts, ensure you have R installed with the following required packages:

```r
install.packages(c("tidyverse", "tidymodels", "here", "vip", "patchwork", 
                   "moments", "rstatix", "ggpubr", "ggprism", "scales", 
                   "RColorBrewer"))
```

## Workflow

### 1. *In Vitro* Analysis

- `01_invitro_cfu.R`: Analyzes bacterial growth based on CFU counts and optical density measurements.
- `01_invitro_cell_classifier.R`: Uses logistic regression and random forest models to classify bacterial cells from *in vitro* data.
- `01_invitro_microscopy.R`: Analyses *in vitro* fluorescence data.

### 2. *In Planta* Analysis

- `02_inplanta_cfu.R`: Evaluates bacterial load in the *A. thaliana* phyllosphere over time.
- `02_inplanta_cell_classifier.R`: Trains machine learning classifiers to distinguish bacterial cell from *in planta* microscopy images.
- `02_inplanta_microscopy.R`: Analyses *in planta* fluorescence data.

## Data Processing

- **Preprocessing**: Label cleanup, filtering, and normalization of fluorescence intensity.
- **Machine Learning**: Logistic regression and random forest models for bacterial classification.
- **Statistical Analysis**: Wilcoxon tests, ANOVA, and normality assessments for fluorescence distributions.
- **Visualization**: Cumulative probability distributions, ROC curves, feature importance plots, and CFU trends.

## Usage

Run the scripts sequentially to process the data and generate results:

```r
source("01_invitro_cfu.R")
source("01_invitro_cell_classifier.R")
source("01_invitro_microscopy.R")
source("02_inplanta_cfu.R")
source("02_inplanta_cell_classifier.R")
source("02_inplanta_microscopy.R")
```

Ensure that input datasets are stored in the `data/` directory with appropriate formatting.

## Results

- Machine learning models stored as `.rds` files in `results/`
- Plots saved in `results/` as `.png`
- Statistical summaries printed in the console
