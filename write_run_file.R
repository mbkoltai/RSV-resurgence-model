# nohup Rscript --vanilla --verbose my_script.R 1 5 1 5 &> nohup_1.out &
#   nohup Rscript --vanilla --verbose my_script.R 1 5 6 10  &> nohup_2.out &
#   nohup Rscript --vanilla --verbose my_script.R 6 10 1 5 &> nohup_3.out &
#   nohup Rscript --vanilla --verbose my_script.R 6 10 6 10 &> nohup_4.out &
arg_vals=as.numeric(commandArgs(trailingOnly=TRUE)[1:4]); 
k_length=arg_vals[1]; n_row=arg_vals[2]; simul_dur=arg_vals[3]; post_npi=arg_vals[4]; 
initcond_filename=commandArgs(trailingOnly=TRUE)[5]
parscan_split <- lapply(1:k_length,function(x) matrix(sort(c(round(seq(1,n_row,length.out=k_length+1)),
                                round(seq(1,n_row,length.out=k_length+1))[2:k_length]+1)),nrow=2)[,x])
string_main_run <- paste0("#!/bin/bash\n",paste0(unlist(lapply(1:k_length,
  function(x) paste0("nohup Rscript --vanilla parscan_runner_cmd_line.R ",
                     paste0(c(parscan_split[[x]],simul_dur,post_npi,initcond_filename),collapse=" "),
  " > simul_output/parscan/parallel/nohup_",x,".out",ifelse(x<k_length," & \n",""),collapse="") )),collapse=""))
# write file
write.table(string_main_run,file="run_parallel_scan.sh",col.names=F,row.names=F,quote=F)

# starter simulations: for each parameter "group" there is a 1st simulation that the others use as init conds, we need to run these first
string_start_run <- paste0("#!/bin/bash\n",paste0(unlist(lapply(1:k_length, 
  function(x) paste0("nohup Rscript --vanilla parscan_starter_cmd_line.R ",paste0(c(x,k_length,2*simul_dur),collapse=" "),
 " > simul_output/parscan/parallel/nohup_starter_",x,".out",ifelse(x<k_length," & \n",""),collapse="") )),
                                                 collapse=""))
# write to file
write.table(string_start_run,file="start_parallel_scan.sh",col.names=F,row.names=F,quote=F)