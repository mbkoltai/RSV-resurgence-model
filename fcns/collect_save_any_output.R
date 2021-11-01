library(readr); library(dplyr)
parall_foldername <- commandArgs(trailingOnly=TRUE)[1] # "simul_output/parscan/parallel/"
file_pattern<-commandArgs(trailingOnly=TRUE)[2]; file_name_save<-commandArgs(trailingOnly=TRUE)[3]
write_csv(bind_rows(lapply(list.files(path=parall_foldername,pattern=file_pattern),
         function(x) read_csv(paste0(parall_foldername,x)))),file=paste0(parall_foldername,file_name_save))
if (grepl("delete|del|remove",commandArgs(trailingOnly=TRUE)[4])){
  unlink(paste0(parall_foldername,list.files(path=parall_foldername,pattern=file_pattern)))
}