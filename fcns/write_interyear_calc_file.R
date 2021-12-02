# write file for checking interyear difference

# Rscript fcns/calc_interyear_diff.R simul_output/parscan/parallel/parsets_filtered_1084_35y/ 2018-09-01 2018-10-01
foldername <- commandArgs(trailingOnly=TRUE)[1]
rscript_command <- paste0("nohup Rscript fcns/calc_interyear_diff.R ",foldername," ",
                    commandArgs(trailingOnly=TRUE)[2]," ",commandArgs(trailingOnly=TRUE)[3])
no_files <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
memory_max <- as.numeric(commandArgs(trailingOnly=TRUE)[5])

master_file <- paste0(c("#!/bin/bash
#$ -N ARRAY_TEST_JOB
#$ -cwd -V
#$ -q short.q
#$ -l mem_free=1G,h_vmem=",memory_max,"G
#$ -t 1-",no_files,"\n",
  rscript_command," ${SGE_TASK_ID} > ",foldername,"nohup_${SGE_TASK_ID}.out &"),collapse="")
  
write.table(paste0(master_file,collapse = "\n"),file="start_batches_calc_interyear.sh",
            col.names=F,row.names=F,quote=F)