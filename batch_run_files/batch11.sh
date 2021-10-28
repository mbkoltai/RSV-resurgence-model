#!/usr/bin/bash 
module load R/3.6.3
Rscript --vanilla fcns/parscan_runner_cmd_line.R 7899 8688 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_11.out 
