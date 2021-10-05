#!/bin/bash
nohup Rscript --vanilla parscan_starter_cmd_line.R 1 7 14 > simul_output/parscan/parallel/nohup_starter_1.out & 
nohup Rscript --vanilla parscan_starter_cmd_line.R 2 7 14 > simul_output/parscan/parallel/nohup_starter_2.out & 
nohup Rscript --vanilla parscan_starter_cmd_line.R 3 7 14 > simul_output/parscan/parallel/nohup_starter_3.out & 
nohup Rscript --vanilla parscan_starter_cmd_line.R 4 7 14 > simul_output/parscan/parallel/nohup_starter_4.out & 
nohup Rscript --vanilla parscan_starter_cmd_line.R 5 7 14 > simul_output/parscan/parallel/nohup_starter_5.out & 
nohup Rscript --vanilla parscan_starter_cmd_line.R 6 7 14 > simul_output/parscan/parallel/nohup_starter_6.out & 
nohup Rscript --vanilla parscan_starter_cmd_line.R 7 7 14 > simul_output/parscan/parallel/nohup_starter_7.out
