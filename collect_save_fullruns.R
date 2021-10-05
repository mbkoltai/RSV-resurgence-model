parall_foldername="simul_output/parscan/parallel/"
write_csv(bind_rows(lapply(list.files(path=parall_foldername,pattern="parsets*"),
              function(x) read_csv(paste0(parall_foldername,x)))),path=paste0(parall_foldername,"results_summ_all.csv"))

unlink(paste0(parall_foldername,list.files(path=parall_foldername,pattern="parsets_start*")))