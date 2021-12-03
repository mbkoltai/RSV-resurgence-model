#!/bin/bash
#$ -N ARRAY_TEST_JOB
#$ -cwd -V
#$ -q short.q
#$ -l mem_free=1G,h_vmem=4G
#$ -t 1-64
Rscript fcns/calc_interyear_diff.R  simul_output/parscan/parsets_filtered_1084_90pct_red/ 2018-09-01 2018-10-10 2020-03-15 42 9 ${SGE_TASK_ID} > simul_output/parscan/parsets_filtered_1084_90pct_red/nohup_interyear${SGE_TASK_ID}.out
