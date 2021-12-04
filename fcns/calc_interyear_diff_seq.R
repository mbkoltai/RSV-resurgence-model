x1<-c("tidyverse","readr")
x2 <- x1 %in% row.names(installed.packages()); if (any(x2 == FALSE)) { install.packages(x1[! x2]) }
library(lubridate); library(dplyr)
# Load all packages    
lapply(x1,library,character.only=TRUE) # as.Date <- zoo::as.Date
foldername<-commandArgs(trailingOnly=TRUE)[1]; # print("FOLDERNAME: "); # print(foldername)
print("cmd line arguments: "); print(commandArgs(trailingOnly = TRUE))
print("# of cmd line args: "); print(length(commandArgs(trailingOnly=TRUE)))
file_list <- list.files(path=foldername,pattern="dyn_parsets*"); # print(file_list)
start_date_dyn_save <- commandArgs(trailingOnly=TRUE)[2]
start_date_calc<-as.Date(commandArgs(trailingOnly=TRUE)[3]) # print("start of season calc"); print(yday_start_end)
stop_date_calc<-as.Date(commandArgs(trailingOnly=TRUE)[4])
start_week <- as.numeric(commandArgs(trailingOnly=TRUE)[5]); stop_week <- as.numeric(commandArgs(trailingOnly=TRUE)[6])
print("starting loop")
# k_file <- as.numeric(commandArgs(trailingOnly=TRUE)[7]); print(file_list[k_file])
for (k_file in 1:length(file_list)){
dyn_df <- read_csv(paste0(foldername,file_list[k_file])); 
cntr=0
for (k_par in unique(dyn_df$par_id)) {
  x <- dyn_df %>% filter(par_id==k_par) %>% mutate(date=as.Date(start_date_dyn_save)+t-min(t)) %>% 
    # usual limits of calc: "2018-10-10" | "2020-03-15"
    filter(date>=start_date_calc & date<=stop_date_calc) %>% mutate(week=week(date)) %>%
    filter(week>=start_week | week<=stop_week) %>%
    mutate(day_of_year=yday(date),epi_year=ifelse(week>=start_week,paste0(year(date),"_",year(date)+1),
                      paste0(year(date)-1,"_",year(date))) ) %>% group_by(agegroup,infection,par_id,day_of_year) %>% 
    summarise(diff_interyr=abs(diff(value)),value=mean(value)) %>% group_by(agegroup,infection,par_id) %>% 
    summarise(cumul_mean_incid=round(sum(value)),sum_abs_diff=round(sum(diff_interyr)),
              sum_rel_diff=round(sum(diff_interyr)/sum(value),4) )
  print(paste0("done, file: ",k_file, ", param: ", k_par, "(",cntr,")"))
  write_csv(x,paste0(foldername,"summ_diff_interyr",k_file,".csv"),append=ifelse(k_par==unique(dyn_df$par_id)[1],F,T))
}
}

write_csv(bind_rows(lapply(list.files(path=foldername,pattern="summ_diff_interyr"),
              function(x) read_csv(paste0(parall_foldername,x)))),file=paste0(foldername,"summ_diff_interyr_all.csv"))
