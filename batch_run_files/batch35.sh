#!/usr/bin/bash 
module load R/3.6.3
Rscript --vanilla fcns/parscan_runner_cmd_line.R 12395 12758 25 4 partable_full.csv data/estim_attack_rates.csv nosave 2018-09-01 > simul_output/parscan/parallel/nohup_35.out 
