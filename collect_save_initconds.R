# collect and save outputs for initconds
library(readr); library(dplyr)
parall_foldername="simul_output/parscan/parallel/"
write_csv(bind_rows(lapply(list.files(path=parall_foldername,pattern="parsets_start*"),
   function(x) read_csv(paste0(parall_foldername,x)))),path=paste0(parall_foldername,"initconds_all.csv"))

unlink(paste0(parall_foldername,list.files(path=parall_foldername,pattern="parsets_start*")))