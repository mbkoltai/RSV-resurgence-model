library(readr); library(dplyr)
parall_foldername <- commandArgs(trailingOnly=TRUE)[1] # "simul_output/parscan/parallel/"
write_csv(bind_rows(lapply(list.files(path=parall_foldername,pattern="summ_parsets*"),
              function(x) read_csv(paste0(parall_foldername,x)))),file=paste0(parall_foldername,"results_summ_all.csv"))
if (grepl("delete|del|remove",commandArgs(trailingOnly=TRUE)[2])){
  unlink(paste0(parall_foldername,list.files(path=parall_foldername,pattern="summ_parsets*")))
}

if (length(commandArgs(trailingOnly=TRUE))>2){
write_csv(bind_rows(lapply(list.files(path=parall_foldername,pattern="dyn_parsets*"),
          function(x) read_csv(paste0(parall_foldername,x)))),
          file=paste0(parall_foldername,"results_dyn_all.csv"))
  }
# unlink(paste0(parall_foldername,list.files(file=parall_foldername,pattern="parsets_start*")))