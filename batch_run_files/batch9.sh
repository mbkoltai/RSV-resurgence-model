#!/usr/bin/bash 
module load R/3.6.3
Rscript --vanilla fcns/parscan_runner_cmd_line.R 6319 7108 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_9.out 
