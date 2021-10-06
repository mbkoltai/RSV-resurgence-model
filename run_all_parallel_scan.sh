#!/bin/bash
nohup Rscript --vanilla parscan_starter_cmd_line.R 1 7 14 > simul_output/parscan/parallel/nohup_starter_1.out & 
nohup Rscript --vanilla parscan_starter_cmd_line.R 2 7 14 > simul_output/parscan/parallel/nohup_starter_2.out & 
nohup Rscript --vanilla parscan_starter_cmd_line.R 3 7 14 > simul_output/parscan/parallel/nohup_starter_3.out & 
nohup Rscript --vanilla parscan_starter_cmd_line.R 4 7 14 > simul_output/parscan/parallel/nohup_starter_4.out & 
nohup Rscript --vanilla parscan_starter_cmd_line.R 5 7 14 > simul_output/parscan/parallel/nohup_starter_5.out & 
nohup Rscript --vanilla parscan_starter_cmd_line.R 6 7 14 > simul_output/parscan/parallel/nohup_starter_6.out

PID=$!

wait $PID

Rscript collect_save_initconds.R

nohup Rscript --vanilla parscan_runner_cmd_line.R 1 31 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_1.out & 
nohup Rscript --vanilla parscan_runner_cmd_line.R 32 61 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_2.out & 
nohup Rscript --vanilla parscan_runner_cmd_line.R 62 91 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_3.out & 
nohup Rscript --vanilla parscan_runner_cmd_line.R 92 120 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_4.out & 
nohup Rscript --vanilla parscan_runner_cmd_line.R 121 150 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_5.out & 
nohup Rscript --vanilla parscan_runner_cmd_line.R 151 180 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_6.out & 
nohup Rscript --vanilla parscan_runner_cmd_line.R 181 210 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_7.out

PID=$!

wait $PID

Rscript collect_save_fullruns.R
