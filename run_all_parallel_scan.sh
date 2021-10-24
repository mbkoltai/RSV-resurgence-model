#!/bin/bash


nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1 118 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_1.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 119 235 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_2.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 236 352 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_3.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 353 469 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_4.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 470 586 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_5.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 587 703 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_6.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 704 820 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_7.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 821 937 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_8.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 938 1054 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_9.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1055 1171 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_10.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1172 1288 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_11.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1289 1404 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_12.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1405 1521 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_13.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1522 1638 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_14.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1639 1755 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_15.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1756 1872 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_16.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1873 1989 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_17.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 1990 2106 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_18.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 2107 2223 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_19.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 2224 2340 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_20.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 2341 2457 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_21.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 2458 2574 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_22.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 2575 2691 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_23.out & 
nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R 2692 2808 15 4 simul_output/parscan/parallel/partable.csv data/estim_attack_rates.csv nosave > simul_output/parscan/parallel/nohup_24.out

PID=$!

wait $PID

Rscript fcns/collect_save_fullruns.R
