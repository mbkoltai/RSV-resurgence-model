# write file for checking interyear difference

# "Rscript fcns/write_interyear_calc_file.R simul_output/parscan/parsets_filtered_1084_90pct_red/  (1)
#                                           2018-09-01 (2)
#                                           2018-10-10 (3)
#                                           2020-03-15 (4)
#                                           42 9 (5-6)
#                                           64 4" (7-8)
foldername <- commandArgs(trailingOnly=TRUE)[1]
rscript_command <- paste0(c("Rscript fcns/calc_interyear_diff.R ",foldername,
        commandArgs(trailingOnly=TRUE)[2],commandArgs(trailingOnly=TRUE)[3],commandArgs(trailingOnly=TRUE)[4],
        commandArgs(trailingOnly=TRUE)[5],commandArgs(trailingOnly=TRUE)[6]),collapse=" ")
no_files <- as.numeric(commandArgs(trailingOnly=TRUE)[7])
memory_max <- as.numeric(commandArgs(trailingOnly=TRUE)[8])

master_file <- paste0(c("#!/bin/bash
#$ -N ARRAY_TEST_JOB
#$ -cwd -V
#$ -q short.q
#$ -l mem_free=1G,h_vmem=",memory_max,"G
#$ -t 1-",no_files,"\n",
  rscript_command," ${SGE_TASK_ID} > ",foldername,"nohup_interyear${SGE_TASK_ID}.out"),collapse="")
  
write.table(paste0(master_file,collapse = "\n"),file="start_batches_calc_interyear.sh",
            col.names=F,row.names=F,quote=F)