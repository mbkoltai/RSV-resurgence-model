### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# library(contactdata); library(fitdistrplus);  library(bbmle); library(Rcpp); library(GillespieSSA)
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","plotly","Rcpp","zoo","lubridate","tsibble","pracma","lubridate",
         "qs","ungeviz"),library,character.only=TRUE)
standard_theme=theme(plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90,vjust=1/2),
  axis.text.y=element_text(size=9),axis.title=element_text(size=14), text=element_text(family="Calibri"))
### UK RSV data  --------------------------------------------------------
### data with weekly resolution (no age resol) ----
season_weeks=c(9,41)
resp_detects_weekly_all_age=read_csv("data/Respiratory viral detections by any method UK.csv") %>% 
  mutate(year_week=factor(paste0(Year,"-",Week),unique(paste0(Year,"-",Week))), RSV_rolling_av=rollmean(RSV,k=7,align="center",fill=NA) ) %>% 
  select(-(contains("virus")|contains("flu"))) %>% mutate(epi_year=ifelse(Week>=season_weeks[2],Year+1,Year)-min(Year)+1) %>% 
  group_by(epi_year) %>% mutate(perc_yearly=RSV/sum(RSV)) %>% group_by(epi_year) %>% 
  mutate(season_share=sum(perc_yearly[Week>=season_weeks[2] | Week<=season_weeks[1]]),
         on_off_season=ifelse(findInterval(Week,season_weeks+c(1,0))==1,"off","on"))
# plot
ggplot(resp_detects_weekly_all_age %>% mutate(section=ceiling((as.numeric(year_week)*0.92)/100)) %>% 
         mutate(section=ifelse(section>3,3,section)),aes(x=year_week)) + geom_point(aes(y=RSV,group=1,color=factor(Year))) +
  geom_line(aes(y=RSV_rolling_av,group=1)) + labs(color="Year") + scale_y_continuous(expand=expansion(0,0.05)) + 
  facet_wrap(~section,nrow=3,scale="free_x") + geom_vline(data=subset(resp_detects_weekly_all_age, Week %in% c(40,11)),
             aes(xintercept=year_week),color="red",linetype="dashed",size=0.5) + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5,size=6))
# save
ggsave("simul_output/uk_rsv_data2020_allagegroups_weekly.png",width=32,height=20,units="cm")

# concentration of cases within 'season'
resp_detects_weekly_all_age_means_shares=left_join(
  resp_detects_weekly_all_age %>% group_by(epi_year,on_off_season) %>% summarise(mean_on_off=mean(RSV)) %>% 
    pivot_wider(names_from=on_off_season,values_from=mean_on_off,names_prefix="mean_"),  
  resp_detects_weekly_all_age %>% group_by(epi_year) %>% filter(epi_year>1&epi_year<8) %>% summarise(season_share=unique(season_share)) ) %>%
  mutate(on_off_ratio=mean_on/mean_off)
write_csv(resp_detects_weekly_all_age_means_shares,"data/resp_detects_weekly_all_age_means_shares.csv")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# with AGE RESOLUTION
resp_virus_data_uk=read_csv("data/Respiratory viral detections by any method UK Ages.csv")
resp_virus_data_uk_tidy = resp_virus_data_uk %>% pivot_longer(!c("Year","startweek","Age")) %>% 
  mutate(Age=factor(gsub(" Y","Y",Age),levels=unique(gsub(" Y","Y",Age))))
leaveout_year=c(2014); truthvals_rsv=resp_virus_data_uk_tidy$name %in% "RSV" & (!resp_virus_data_uk_tidy$Year %in% leaveout_year)
# means across years
averages_years=data.frame(resp_virus_data_uk_tidy[truthvals_rsv,] %>% group_by(startweek,Age) %>% 
                            summarise(value=mean(value)),Year='mean',name='RSV')
if (!any(resp_virus_data_uk_tidy$Year %in% "mean")){
  resp_virus_data_uk_tidy=rbind(resp_virus_data_uk_tidy,averages_years[,colnames(resp_virus_data_uk_tidy)])}
truthvals_rsv=(resp_virus_data_uk_tidy$name %in% "RSV") & resp_virus_data_uk_tidy$Year %in% c("mean",2020)
# (!resp_virus_data_uk_tidy$Year %in% c(leaveout_year,))
resp_virus_data_uk_tidy[,"type"]="indiv year"; resp_virus_data_uk_tidy[resp_virus_data_uk_tidy$Year %in% "mean","type"]="5-year average"
resp_virus_data_uk_tidy[,"width"]=1.01; resp_virus_data_uk_tidy[resp_virus_data_uk_tidy$Year %in% "mean","width"]=1.015
# convert weel number to date
resp_virus_data_uk_tidy[,"date"]=as.Date(paste(resp_virus_data_uk_tidy$Year,resp_virus_data_uk_tidy$startweek,1,sep="-"),"%Y-%U-%u")
resp_virus_data_uk_tidy[,"year_week"]=gsub("mean-","",paste0(resp_virus_data_uk_tidy$Year,"-",resp_virus_data_uk_tidy$startweek))
resp_virus_data_uk_tidy$year_week=factor(resp_virus_data_uk_tidy$year_week,unique(resp_virus_data_uk_tidy$year_week))
# PLOT all years
ggplot(subset(resp_virus_data_uk_tidy,name %in% "RSV" & grepl("indiv",type)),aes(x=date,y=value,group=Age)) + 
  geom_area(aes(fill=Age),position=position_stack(reverse=T),color="black",size=1/4) +
  # geom_line(aes(color=Age)) + geom_point(size=0.4) + facet_wrap(~Age,ncol=2,scales="free") + # 
  scale_x_date(breaks="2 month",expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.01,0)) + 
  ylab("number of reported RSV cases") + theme_bw() + standard_theme + xlab("")  + labs(fill="") +
  theme(axis.text.x=element_text(size=13),axis.text.y=element_text(size=14),legend.position="bottom",legend.text=element_text(size=15),
        axis.title.y=element_text(size=16))
#  labs(caption="source: gov.uk/health-and-social-care/health-protection-infectious-diseases")
ggsave("data/uk_rsv_data2014_2020.png",width=32,height=22,units="cm")
# ggsave("simul_output/uk_rsv_data2014_2020_age_facet.png",width=32,height=16,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Datamart with number of tests and positives
datamart_weekly_rsv_data <- read_csv("data/datamart_weekly_rsv_data_2012wk40-2021wk24.csv") %>% 
  mutate(`Week Commencing`=dmy(`Week Commencing`)) %>% pivot_longer(!c(`Week Commencing`,`Week No`)) %>%
  mutate(categ=case_when(grepl("RSV %",name) ~ "positivity", grepl("RSV pos",name) ~ "cases", grepl("RSV tot",name) ~ "tests"),
    age=case_when(grepl("0-<05",name) ~ "0-5", grepl("05-14",name) ~ "5-14", grepl("15-44",name) ~ "15-44",
            grepl("45-64",name) ~ "45-64", grepl("65+",name) ~ "65+"),age=ifelse(is.na(age),"total",age)) %>% 
  mutate(age=factor(age,levels=unique(age)))
# library(wesanderson)
# all years in a row
# ggplot(datamart_weekly_rsv_data %>% filter(name %in% "RSV tot_Total"),aes(x=`Week Commencing`,y=value)) +
#        geom_line() + geom_point(shape=21) +
#   geom_vline(data=datamart_weekly_rsv_data %>% filter(`Week No`==1 & name %in% "RSV tot_Total"),
#              aes(xintercept=`Week Commencing`),color="red",linetype="dashed",size=0.3) +
#   geom_vline(data=datamart_weekly_rsv_data %>% filter(`Week No`==36 & name %in% "RSV tot_Total"),
#              aes(xintercept=`Week Commencing`),color="gray45",size=0.3) +
#   scale_x_date(date_breaks="2 months",expand=expansion(0.01,0)) + theme_bw() + standard_theme + 
#   theme(axis.text.x=element_text(vjust=1/2)) + ggtitle("number of tests")
# # save
# ggsave(paste0("data/test_number_allyears.png"),width=38,height=22,units="cm")
###
# test vs case numbers
max_vals=datamart_weekly_rsv_data %>% group_by(name) %>% summarise(max_val=max(value,na.rm=T))
# plot
df_plot<-datamart_weekly_rsv_data %>% group_by(name) %>% 
  mutate(norm_val_basel_subtr=(value-min(value,na.rm=T))/(max(value,na.rm=T)-min(value,na.rm=T)),norm_val=value/max(value,na.rm=T))
###
# plot total
ggplot(df_plot %>% filter(age=="total" & categ != "positivity"), aes(x=`Week Commencing`,y=norm_val,color=categ) ) + geom_line() + 
  geom_point(data=df_plot %>% filter(grepl("tests",name)),shape=21) + labs(color="") +
  geom_vline(data=datamart_weekly_rsv_data %>% filter(`Week No`==50),aes(xintercept=`Week Commencing`),
             color="black",size=0.2,linetype="dashed") + 
  scale_x_date(date_breaks="2 months",expand=expansion(0.005,0)) + scale_y_continuous(breaks=(0:10)/10,expand=expansion(0.02,0)) + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=1/2),
      legend.text=element_text(size=16),legend.position="bottom") + ggtitle("tests vs positives") + ylab("normalised value") + xlab("")
# save
ggsave(paste0("data/tests_cases_allyears_normalised.png"),width=40,height=22,units="cm") # _baseline_subtr

# positivity and test number. by week, color-coded by year
ggplot(datamart_weekly_rsv_data %>% 
         mutate(year=ifelse(`Week No`==1 & month(`Week Commencing`)==12,year(`Week Commencing`)+1,year(`Week Commencing`))) %>%
         filter(categ %in% c("positivity","tests") & age=="total" & year>2012),
       aes(x=`Week No`,y=ifelse(grepl("tot_",name),value/1e2,value),color=categ)) + labs(color="") +
  geom_line() + geom_point(shape=21,size=0.5) + theme_bw() + standard_theme + ylab("test number (x100) / positivity %") + 
  facet_wrap(~year) + theme(legend.position="top",legend.text=element_text(size=14)) + 
  scale_y_continuous(expand=expansion(0.01,0)) + scale_x_continuous(breaks=(1:26)*2,expand=expansion(0.01,0))
# save
ggsave(paste0("data/tests_positivity_year_week.png"),width=32,height=22,units="cm")

### ### ### ### ### ### ### ### ### ### ###
# by age
datamart_plot=datamart_weekly_rsv_data %>% 
  mutate(year=ifelse(`Week No`==1 & month(`Week Commencing`)==12,year(`Week Commencing`)+1,year(`Week Commencing`))) %>%
  filter(categ %in% c("cases") & `Week Commencing`>as.Date("2013-06-01") & age!="total")
max_weeks=datamart_plot %>% filter(year<2020) %>% group_by(year,age) %>% 
  summarise(max_val=max(value,na.rm=T),week=unique(`Week No`[value==max(value,na.rm=T)]),
            `Week Commencing`=unique(`Week Commencing`[value==max(value,na.rm=T)])) %>% filter(!is.na(week))
# 
ggplot(datamart_plot,aes(x=`Week Commencing`,y=value)) + 
  geom_line() + geom_point(aes(y=ifelse(categ=="tests",NA,value)),shape=21,size=0.5) + #,color=categ
  facet_wrap(~age,scales="free") + theme(legend.position="top",legend.text=element_text(size=14)) + theme_bw() + standard_theme + 
  geom_vline(data=max_weeks,aes(xintercept=`Week Commencing`),color="red",size=1/4,linetype="dashed") +
  scale_x_date(date_breaks="6 month",expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.01,0)) + 
  xlab("") + ylab("# cases") + labs(color="")
# save
# ggsave(paste0("data/datamart_cases_tests_by_age.png"),width=32,height=22,units="cm")
ggsave(paste0("data/datamart_cases_by_age.png"),width=32,height=22,units="cm")

# plot all years on the same x axis
ggplot(datamart_plot %>% filter(year<2020) %>% group_by(year,age) %>% mutate(norm_val=value/max(value,na.rm=T)),
       aes(x=`Week No`,y=norm_val,color=factor(year))) + geom_line() + geom_point() +
  facet_wrap(~age,scales="free",ncol=2) + theme(legend.position="top",legend.text=element_text(size=14)) + theme_bw() + standard_theme + 
  geom_vline(data=max_weeks,aes(xintercept=week,color=factor(year)),size=1/3,linetype="dashed") +
  scale_x_continuous(breaks=1+(0:26)*2,expand=expansion(0.01,0)) + scale_y_continuous(breaks=(0:10)/10,expand=expansion(0.01,0)) + 
  xlab("") + ylab("# cases") + labs(color="")
# save
ggsave(paste0("data/datamart_cases_by_age_years_together_norm.png"),width=32,height=22,units="cm")

# plot normalised average
datamart_plot_multiyear_aver = datamart_weekly_rsv_data %>% 
  mutate(year=ifelse(`Week No`==1 & month(`Week Commencing`)==12,year(`Week Commencing`)+1,year(`Week Commencing`))) %>%
  filter(categ %in% c("cases") & `Week Commencing`>as.Date("2013-06-01") & age!="total") %>% 
  filter(year<2020) %>% group_by(year,age) %>% mutate(norm_val=value/max(value,na.rm=T)) %>% group_by(age,`Week No`) %>%
  summarise(mean_norm_val=mean(norm_val,na.rm=T),ci50_low_norm_val=quantile(norm_val,probs=c(0.25,0.75),na.rm=T)[1], 
            ci50_high_norm_val=quantile(norm_val,probs=c(0.25,0.75),na.rm=T)[2],
            mean=mean(value,na.rm=T),ci50_low=quantile(value,probs=c(0.25,0.75),na.rm=T)[1],
              ci50_high=quantile(value,probs=c(0.25,0.75),na.rm=T)[2])
# PLOT
ggplot(datamart_plot_multiyear_aver,aes(x=`Week No`,y=mean_norm_val)) + geom_line() + geom_point() + # _norm_val
  geom_ribbon(aes(ymin=ci50_low_norm_val,ymax=ci50_high_norm_val),fill="pink",alpha=0.3) + facet_wrap(~age,scales="free",ncol=2) + 
  theme(legend.position="top",legend.text=element_text(size=14)) + theme_bw() + standard_theme +
  geom_vline(data=max_weeks,aes(xintercept=week,color=factor(year)),size=1/3,linetype="dashed") +
  geom_hline(yintercept=0.1,size=1/5,linetype="dashed") + xlab("") + ylab("# cases") + labs(color="") +
  scale_x_continuous(breaks=1+(0:26)*2,expand=expansion(0.01,0)) + scale_y_continuous(breaks=(0:10)/10,expand=expansion(0.01,0))
# maximal week on average? 49 # 10% of max: week 5-9, week 41
# max_weeks %>% mutate(week=ifelse(week<20,week+52,week)) %>% group_by(age) %>% summarise(mean(week))
ggsave(paste0("data/datamart_cases_by_age_years_together_norm_mean.png"),width=32,height=22,units="cm")
# ggsave(paste0("data/datamart_cases_by_age_years_together_mean.png"),width=32,height=22,units="cm")

### ### ### ### ### ### ### ### ### ### ###
# compare datamart with Resp Infs data
# head(datamart_weekly_rsv_data %>% filter(`Week No`>45 & year(`Week Commencing`)==2014))
# head(resp_detects_weekly_all_age %>% filter(Year==2014))

datamart_resp_infs_combined <- right_join(resp_detects_weekly_all_age %>% rename(year=Year,week=Week,RSV_pos_resp_infs=RSV) %>% 
                        select(c(year,week,RSV_pos_resp_infs)), 
                        read_csv("data/datamart_weekly_rsv_data_2012wk40-2021wk24.csv") %>% 
                          mutate(`Week Commencing`=dmy(`Week Commencing`),year=year(`Week Commencing`)) %>%
                        select(`Week Commencing`,year,`Week No`,`RSV pos_Total of LabPK`,`RSV tot_Total`) %>% 
              rename(week=`Week No`,RSV_pos_datamart=`RSV pos_Total of LabPK`,RSV_total_datamart=`RSV tot_Total`), by=c("year","week"))
seas_lims_data=data.frame(lapply(c(11,44), function(x) 
  datamart_resp_infs_combined$`Week Commencing`[datamart_resp_infs_combined$week==x])); colnames(seas_lims_data)=c("on","off")
# plot
datamart_resp_infs_combined %>% filter(`Week Commencing`>as.Date("2014-06-01")) %>% pivot_longer(!c(week,year,`Week Commencing`,epi_year)) %>% 
  mutate(categ=ifelse(grepl("pos",name),"positives (normalised)","tests (datamart)")) %>% group_by(name) %>% 
  mutate(norm_value=ifelse(grepl("pos",name),value/max(value,na.rm=T),value)) %>% filter(grepl("pos",name)) %>%
  ggplot(aes(x=`Week Commencing`,y=norm_value,color=name)) + geom_line() + geom_point(shape=21) + # ifelse(value>=1,value,NA)
  facet_wrap(~categ,scales="free",nrow=2) + scale_x_date(date_breaks="2 month",expand=expansion(0.01,0)) + 
  scale_y_continuous(expand=expansion(0.01,0)) + theme_bw() + standard_theme +
  # geom_vline(data=datamart_weekly_rsv_data %>% filter(`Week No` %in% c(44,8) & categ=="cases" &age=="total"),
  #    aes(xintercept=`Week Commencing`),size=0.2,linetype="dashed") + 
  xlab("") + ylab("") + labs(color="") + theme(legend.position="top")
# SAVE
ggsave(paste0("data/datamart_resp_infs_comparison_norm_val_positives.png"),width=32,height=22,units="cm")
# ggsave(paste0("data/datamart_resp_infs_comparison_positives.png"),width=32,height=22,units="cm")
# ggsave(paste0("data/datamart_resp_infs_comparison_norm_val.png"),width=32,height=22,units="cm")

# phase plot: posit-test number
# ggplot(datamart_weekly_rsv_data %>% 
#          mutate(year=ifelse(`Week No`==1 & month(`Week Commencing`)==12,year(`Week Commencing`)+1,year(`Week Commencing`))) %>%
#          select(year,`Week No`,`RSV %_Total`,`RSV tot_Total`),
#        aes(x=`RSV %_Total`,y=`RSV tot_Total`,color=`Week No`)) + labs(color="week") +
#   geom_path() + geom_point(shape=21,size=0.5) + theme_bw() + standard_theme + xlab("positivity %") + ylab("test number") + 
#   facet_wrap(~year,scales="free") + theme(legend.position="top") + 
#   scale_color_gradientn(colours=wes_palette("Zissou1",60,type="continuous")) + 
#   scale_y_continuous(expand=expansion(0.01,0)) + scale_x_continuous(breaks=(1:26)*2,expand=expansion(0.01,0))
# # save
# ggsave(paste0("data/tests_positivity_year_week_phaseplot.png"),width=32,height=22,units="cm")
