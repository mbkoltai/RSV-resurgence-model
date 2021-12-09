# Determinants of RSV epidemiology following  suppression through pandemic contact restrictions

Code for modelling how the patterns of RSV resurgence after COVID-19-related Non-Pharmaceutical Interventions (NPIs) depend RSV epidemiology parameters. The code accompanies the article available at [to_be_added]().

**Figures in the main text and SI can be reproduced by the script [reproduce_results.R](https://github.com/mbkoltai/RSV-model/blob/master/reproduce_results.R).**

## Files and folders:

-  _fcns/_: functions for simulations and plotting
-  _data/_: data files used for model set-up or data plots
-  _repo_data/_: data files from simulations to reproduce figures in the article.
-  **[reproduce_results.R](https://github.com/mbkoltai/RSV-model/blob/master/reproduce_results.R): reproduce figures of the article**
-  _load_params.R_: loads constant parameters of the model
-   _param_scan_model_checks.R_: scripts to perform parameter scans and create plots; contains some extra scripts for which inputs might not be available. For reproducing results use [reproduce_results.R](https://github.com/mbkoltai/RSV-model/blob/master/reproduce_results.R). (not maintained)
-  _RSV_model.R_: script to run individual simulations (not maintained)
-  _UK_rsv_data.R_: scripts to plot UK RSV data

