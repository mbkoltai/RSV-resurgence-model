#!/bin/bash
nohup Rscript --vanilla parscan_runner_cmd_line.R 1 31 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_1.out & 
nohup Rscript --vanilla parscan_runner_cmd_line.R 32 61 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_2.out & 
nohup Rscript --vanilla parscan_runner_cmd_line.R 62 91 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_3.out & 
nohup Rscript --vanilla parscan_runner_cmd_line.R 92 120 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_4.out & 
nohup Rscript --vanilla parscan_runner_cmd_line.R 121 150 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_5.out & 
nohup Rscript --vanilla parscan_runner_cmd_line.R 151 180 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_6.out & 
nohup Rscript --vanilla parscan_runner_cmd_line.R 181 210 7 3 simul_output/parscan/parallel/initconds_all.csv > simul_output/parscan/parallel/nohup_7.out
