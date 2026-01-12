# Affective EEG Microstate Syntax — Task Pipeline

This repository contains MATLAB code to reproduce the analyses for:  
**“Task-evoked EEG microstate syntax complements static metrics in bipolar disorder.”**

## Overview

### Analyses (Q1–Q3)
- **Q1:** Sequence-independence tests for transition rates (Poisson IRRs vs. independence offsets).
- **Q2:** Age-anchored rate-ratio contrasts vs. healthy controls with BH–FDR within edge families.
- **Q3:** Trial-level link between selected transition features (e.g., P300 Anchor-2 outflow) and response time on valenced trials; Holm-adjusted confirmatory tests; LPP as negative control.


- **Goal.** Quantify static microstate metrics (duration, coverage, occurrence) and **directed transition syntax** during an affective decision task; test group (BD-I, BD-II, siblings, HC), condition (negative/neutral/positive), and mood-state effects.; test **Q1** (independence), **Q2** (age-anchored contrasts vs HC), and **Q3** (trial-level behavior link).
- **Windows.** N200 (180–300 ms), P300 (300–500 ms), LPP (500–1000 ms).
- **Templates.** 7-class Custo2017 topographies; polarity ignored for clustering; back-fitting with temporal smoothing.
- **Statistics.**
  - **Metrics:** linear mixed-effects (LME) with transformations (logit for coverage, log(rate+c) for occurrence), omnibus LRTs, back-transformed EMMs/contrasts with BH–FDR. See `lme_posthoc_21FDR_transformsBT.m`.
  - **Transitions:** Q1 independence tests (Poisson IRRs vs sequence-independence null) and Q2 **age-anchored** rate-ratio contrasts vs. HC; FDR within principled edge families. Display routines re-compute BH–FDR per family at plot time.  
- **Outputs.** Publication figures (violin+EMM bars; circular graphs), and tables (Omnibus, Age, EMMs, Pairwise).

> Reproducibility details: code re-computes EMMs and FDR exactly as reported in the manuscript (see `make_report_tables.m`, `print_summary_measure.m`).

## Dependencies

- **MATLAB** R2021b+ (tested with R2024a)
  - Statistics and Machine Learning Toolbox (for `fitlme`, `anova`)
- **EEGLAB** (no-GUI usage supported)
- **MicrostateLab** plugin (Custo2017 templates)
- Helper functions (included): `daviolinplot_microstate_task.m`, `drawERPGridFigures.m`, etc.

## Installation

1. Install EEGLAB and MicrostateLab. Note the plugin root path.
2. Clone this repo.
3. Create a project-specific `config/config.m` from `config/config_template.m`:
   - set `project_root`, `data_raw`, `data_deriv`, `results_dir`, and the path to your EEGLAB/MicrostateLab installation.
4. Place **preprocessed & epoched** EEG files in `data/derivatives/…` as referenced by `Microstate_analysis.m`.

## Data layout

The scripts expect preprocessed `.set` files organized by **Group/Subject/Condition**. See `data/README.md` for the expected folder structure. You can optionally mirror BIDS-EEG (`sub-*/ses-*/eeg/`) but it is not required to run.

## Quick start

```matlab
% 1) Configure
edit config/config.m

% 2) Run the full pipeline
cd code/matlab
run_all

% 3) Key outputs
% - results/figures/*  : metrics plots and circular transition graphs
% - results/tables/*   : Omnibus/Age/EMM/Pairwise tables (CSV)
% - results/stats/*    : MATLAB .mat with OUT_* structs
