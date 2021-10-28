#!/usr/bin/bash 
module load R/3.6.3
Rscript --vanilla fcns/parscan_runner_cmd_line.R 1178 1255 15 4 simul_output/parscan/parallel/partable_filtered.csv data/estim_attack_rates.csv SAVE > simul_output/parscan/parallel/nohup_16.out
