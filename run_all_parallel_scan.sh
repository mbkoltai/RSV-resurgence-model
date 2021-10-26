#!/bin/bash


nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1 1054 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_1.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1055 2107 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_2.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 2108 3160 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_3.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 3161 4213 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_4.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 4214 5266 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_5.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 5267 6318 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_6.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 6319 7371 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_7.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 7372 8424 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_8.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 8425 9477 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_9.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 9478 10530 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_10.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 10531 11583 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_11.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 11584 12636 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_12.out

PID=$!

wait $PID

Rscript fcns/collect_save_fullruns.R
