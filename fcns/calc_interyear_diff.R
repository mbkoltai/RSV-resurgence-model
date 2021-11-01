x1<-c("tidyverse","readr")
x2 <- x1 %in% row.names(installed.packages()); if (any(x2 == FALSE)) { install.packages(x1[! x2]) }
library(lubridate); library(dplyr)
# Load all packages    
lapply(x1,library,character.only=TRUE) # as.Date <- zoo::as.Date
foldername<-commandArgs(trailingOnly=TRUE)[1]
file_list <- list.files(path=foldername,pattern="dyn_parsets*")
start_date_dyn_save <- commandArgs(trailingOnly=TRUE)[2]
yday_start_end<-yday(as.Date(commandArgs(trailingOnly=TRUE)[3])); print(yday_start_end)

k_file <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
  dyn_df <- read_csv(paste0(foldername,file_list[k_file])); print(file_list[k_file])
  for (k_par in unique(dyn_df$par_id)) {
  x <- dyn_df %>% filter(par_id==k_par) %>% mutate(date=as.Date(start_date_dyn_save)+t-min(t)) %>% 
    filter(date>=as.Date(start_date_dyn_save) & date<=as.Date("2020-05-01")) %>%
    mutate(day_of_year=yday(date),epi_year=ifelse(day_of_year>=yday_start_end,paste0(year(date),"_",year(date)+1),
                    paste0(year(date)-1,"_",year(date))) ) %>% group_by(agegroup,infection,par_id,day_of_year) %>% 
    summarise(diff_interyr=abs(diff(value)),value=mean(value)) %>% group_by(agegroup,infection,par_id) %>% 
    summarise(cumul_mean_incid=round(sum(value)),
              sum_abs_diff=round(sum(diff_interyr)),sum_rel_diff=round(sum(diff_interyr)/sum(value),4) )
  print(paste0("done, file: ",k_file, ", param: ", k_par))
  write_csv(x,paste0(foldername,"summ_diff_interyr",k_file,".csv"),
            append=ifelse(k_par==unique(dyn_df$par_id)[1],F,T))
  }
