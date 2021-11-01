x1<-c("tidyverse","readr")
x2 <- x1 %in% row.names(installed.packages()); if (any(x2 == FALSE)) { install.packages(x1[! x2]) }
library(lubridate); library(dplyr)
# Load all packages    
lapply(x1,library,character.only=TRUE) # as.Date <- zoo::as.Date
foldername<-commandArgs(trailingOnly=TRUE)[1]
file_list <- list.files(path=foldername,pattern="dyn_parsets*")
start_date_dyn_save <- commandArgs(trailingOnly=TRUE)[2]
yday_start_end<-yday(as.Date(commandArgs(trailingOnly=TRUE)[3])); print(yday_start_end)
for(k_par in 1:length(file_list)){
  dyn_df <- read_csv(paste0(foldername,file_list[k_par])); print(file_list[k_par])
  x <- dyn_df %>% mutate(date=as.Date(start_date_dyn_save)+t-min(t)) %>% 
    mutate(day_of_year=yday(date),epi_year=ifelse(day_of_year>=yday_start_end,paste0(year(date),"_",year(date)+1),
                    paste0(year(date)-1,"_",year(date))) ); print("created cols")
  x <- x %>% group_by(agegroup,infection,par_id,day_of_year) %>% 
    summarise(diff_interyr=abs(diff(value)),value=mean(value)) %>% group_by(agegroup,infection,par_id) %>% 
    summarise(cumul_mean_incid=sum(value),sum_abs_diff=sum(diff_interyr),sum_rel_diff=sum(diff_interyr)/sum(value))
  print("done, param: ",k_par)
  write_csv(x,paste0(foldername,"summ_diff_interyr.csv"),append=ifelse(k_par==1,F,T))
  # if (k_par %% 10 == 0) {     print(k_par) }
}