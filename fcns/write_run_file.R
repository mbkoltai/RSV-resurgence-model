# nohup Rscript --vanilla --verbose my_script.R 1 5 1 5 &> nohup_1.out &
# "Rscript fcns/write_run_file.R 64 36450 25 4 0.5 repo_data/partable_full.csv NOSAVE sep_qsub_files 2018-09-01 8"

arg_vals=as.numeric(commandArgs(trailingOnly=TRUE)[1:4])
n_core=arg_vals[1]; n_row=arg_vals[2]; simul_dur=arg_vals[3]; 
post_npi=arg_vals[4]
contact_red <- commandArgs(trailingOnly=TRUE)[5]
partable_filename <- commandArgs(trailingOnly=TRUE)[6]
saveflag <- commandArgs(trailingOnly=TRUE)[7]; 
sep_flag <- commandArgs(trailingOnly=TRUE)[8]
start_date_dyn_save <- commandArgs(trailingOnly=TRUE)[9]; 
agegroup_res <- commandArgs(trailingOnly=TRUE)[10]
memory_max <- commandArgs(trailingOnly=TRUE)[11]
parscan_split <- lapply(1:n_core,function(x) matrix(sort(c(round(seq(1,n_row,length.out=n_core+1)),
                                round(seq(1,n_row,length.out=n_core+1))[2:n_core]+1)),nrow=2)[,x])

command_list <- unlist(lapply(1:n_core, # "#!/bin/bash\n",
                      function(x) 
                            paste0("nohup Rscript --vanilla fcns/parscan_runner_cmd_line.R",
                                  paste0(c(parscan_split[[x]],simul_dur,post_npi,
                                           contact_red,partable_filename,saveflag,
                                           start_date_dyn_save,agegroup_res),collapse=" "),
                                  " > simul_output/parscan/nohup_",x,".out",
                                  ifelse(x<n_core," & \n",""),collapse="") ))

if (grepl("sep",sep_flag)) {
  qsub_start_rows <- "#!/usr/bin/bash \nmodule load R/3.6.3\n"
  # full_strings <- paste0(qsub_start_rows,command_list); 
  full_strings <- gsub("& \n","",gsub("nohup ","", paste0(qsub_start_rows,command_list)))
  # write shell files for launching Rscript with different arguments
  for (k in 1:length(full_strings)) { 
    write.table(full_strings[k],file=paste0("batch_run_files/batch",k,".sh",collapse = ""),
                col.names=F,row.names=F,quote=F) }
  # master sh file that'll launch all the sh files
  master_file <- unlist(lapply(1:length(full_strings), function(x_k)
    paste0("qsub -V -cwd -M lshmk17@lshtm.ac.uk -m ea -N batch",x_k," -l mem_free=1G,h_vmem=",memory_max,
           "G -q short.q batch",x_k,".sh",collapse="")))
  write.table(paste0(master_file,collapse = "\n"),file="batch_run_files/start_batches.sh",col.names=F,row.names=F,quote=F)
} else {
write.table(paste0(c("#!/bin/bash\n",paste0(paste0(command_list,collapse="")),"PID=$!","wait $PID",
  "Rscript fcns/collect_save_fullruns.R"),collapse="\n\n"),file="run_all_parallel_scan.sh",col.names=F,row.names=F,quote=F)
}
