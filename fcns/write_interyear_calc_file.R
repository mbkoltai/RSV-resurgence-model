# write file for checking interyear difference

rscript_command <- paste0("Rscript fcns/calc_interyear_diff.R simul_output/parscan/parallel/", 
                    commandArgs(trailingOnly=TRUE)[1],commandArgs(trailingOnly=TRUE)[2], collapse = " ")
no_files <- as.numeric(commandArgs(trailingOnly=TRUE)[3])
memory_max <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
for (k in 1:no_files) {
  
write.table(paste0(rscript_command," ",k),
      file=paste0("batch_run_files/batch_calc_interyear",k,".sh",collapse = ""),col.names=F,row.names=F,quote=F)
}

# master sh file that'll launch all the sh files
master_file <- unlist(lapply(1:no_files, function(x_k)
  paste0("qsub -V -cwd -M lshmk17@lshtm.ac.uk -m ea -N batch_calc_interyear",x_k," -l mem_free=1G,h_vmem=",memory_max,
         "G -q short.q batch",x_k,".sh",collapse="")))
write.table(paste0(master_file,collapse = "\n"),file="batch_run_files/start_batches_calc_interyear.sh",
            col.names=F,row.names=F,quote=F)