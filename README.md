# Determinants of RSV resurgence after the COVID-19-related non-pharmaceutical interventions

Code for modelling RSV resurgence after COVID-19-related Non-Pharmaceutical Interventions (NPIs) and how the patterns of resurgence can be used to infer RSV epidemiology parameters.
**Figures in the main text and SI can be reproduced by the script [reproduce_results.R](https://github.com/mbkoltai/RSV-model/blob/master/reproduce_results.R).**

## Files and folders:

-   _param_scan_model_checks.R_: scripts to perform parameter scans and create plots (contains some extra scripts for which inputs might not be available. For reproducing results use [reproduce_results.R](https://github.com/mbkoltai/RSV-model/blob/master/reproduce_results.R).)
-  _RSV_model.R_: script to run individual simulations
-  _UK_rsv_data.R_: scripts to plot UK RSV data
-  _load_params.R_: loads constant parameters of the model
-  _fcns/_: functions for simulations and plotting
-  _data/_: data files used for model set-up or data plots
-  _repo_data/_: outputs from simulations required to reproduce figures in the article. Currently the full simulation outputs for the filtered paramete sets are uploaded to Google Drive, [available here](https://drive.google.com/file/d/12ohuGEPrVnOxazXnxEGZGwJIwj16frCc/view?usp=sharing).