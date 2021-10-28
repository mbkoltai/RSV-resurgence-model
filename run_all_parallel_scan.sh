#!/bin/bash


nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1 1580 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_1.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1581 3160 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_2.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 3161 4739 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_3.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 4740 6318 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_4.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 6319 7898 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_5.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 7899 9477 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_6.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 9478 11057 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_7.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 11058 12636 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_8.out

PID=$!

wait $PID

Rscript fcns/collect_save_fullruns.R
