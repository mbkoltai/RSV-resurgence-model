#!/bin/bash


nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1 75 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_1.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 76 149 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_2.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 150 223 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_3.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 224 298 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_4.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 299 372 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_5.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 373 446 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_6.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 447 520 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_7.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 521 594 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv > simul_output/parscan/parallel/nohup_8.out

PID=$!

wait $PID

Rscript fcns/collect_save_fullruns.R
