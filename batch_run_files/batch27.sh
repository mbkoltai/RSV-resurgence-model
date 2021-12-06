#!/usr/bin/bash 
module load R/3.6.3
Rscript --vanilla fcns/parscan_runner_cmd_line.R 5135 5331 25 4 0.5 repo_data/partable_full_linear_kage_kexp.csv NOSAVE 2016-09-01 broad_age > simul_output/parscan/nohup_27.out 
