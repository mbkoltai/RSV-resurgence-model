#!/bin/bash
#$ -N ARRAY_TEST_JOB
#$ -cwd -V
#$ -q short.q
#$ -l mem_free=1G,h_vmem=4G
#$ -t 1-64

Rscript fcns/calc_interyear_diff.R FOLDER 2018-09-01 2018-10-01 ${SGE_TASK_ID} > FOLDER/nohup_interyear${SGE_TASK_ID}.out

