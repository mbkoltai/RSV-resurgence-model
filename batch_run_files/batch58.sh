#!/usr/bin/bash 
module load R/3.6.3
Rscript --vanilla fcns/parscan_runner_cmd_line.R 17994 18309 25 4 0.95 repo_data/partable_full_lhs.csv SAVE 2016-09-01 broad_age > simul_output/parscan/nohup_58.out 
