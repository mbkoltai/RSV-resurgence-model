#!/bin/bash


nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1 38 15 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_1.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 39 75 15 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_2.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 76 112 15 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_3.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 113 149 15 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_4.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 150 186 15 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_5.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 187 223 15 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_6.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 224 260 15 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_7.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 261 297 15 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_8.out

PID=$!

wait $PID

Rscript fcns/collect_save_fullruns.R
