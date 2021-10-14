#!/bin/bash


nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1 9 10 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_1.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 10 17 10 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_2.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 18 25 10 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_3.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 26 32 10 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_4.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 33 40 10 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_5.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 41 48 10 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_6.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 49 56 10 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_7.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 57 64 10 3 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_8.out

PID=$!

wait $PID

Rscript fcns/collect_save_fullruns.R
