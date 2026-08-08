
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Project Work Architecture

**Biostatistica per la Ricerca Clinica e la Pubblicazione Scientifica
(BRCPS26 @UNIPD)**

pw_data_analysis_pipeline/

├─ \*.Rproj

├─ .gitignore

├─ main.R

├─ data_input/

├─ output/

└─ R/

    ├─ 00_utils.R                  # global utilities

    ├─ 01_data_prep.R              # load and data cleaning

    ├─ 02_propensity_score.R       # model estimation PS, weigths and balance 

    ├─ 03_outcome_models.R         # model outcome G-computation, Bootstrap

    ├─ 04_spatial_data_analysis.R  # fns maps

    ├─ 05_plots.R                  # fns plots

    ├─ 06_trend_mk_sen.R           # fns trend mann-kenndal, ses's slope NO2

    └─ 07_pop_analysis.R           # fns demography

Folders “data_input/” and “output/” are both “gitignored” but fully
available in local repo.
