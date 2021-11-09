rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# load constant parameters and functions
source("load_params.R") # ; library(wesanderson)
# options(dplyr.summarise.inform=FALSE)
# estimated attack rates
estim_attack_rates <- data.frame(agegroup_name=rsv_age_groups$agegroup_name, # paste0("age=",,"yr")
      median_est=c(rep(65,4),rep(40,4),10,8,5)) %>% mutate(min_est=median_est*0.25,max_est=median_est*2.5,
      median_all_inf=c(rep(70,4),rep(60,4),50,30,20),min_est_all_inf=median_all_inf/2,max_est_all_inf=median_all_inf*2)
# write estim_attack_rates # write_csv(estim_attack_rates,file="data/estim_attack_rates.csv")
# % cases within season (filtering parameter sets)
# seas_conc_lim=0.8
npi_dates<-as.Date(c("2020-03-26","2021-05-17"))
partable <- bind_rows(expand.grid(list(exp_dep=seq(1/4,2,1/8),age_dep=seq(1/8,1,1/16),seasforc_width_wks=c(3,5,7),
      R0=1+(0:5)/10,peak_week=c(48),seasforce_peak=c(3/4,1,5/4),omega=c(1/250,1/350,1/450)) ) )
l_delta_susc <- lapply(1:nrow(partable), function(n_p) {sapply(1:n_age,
                function(x) {(1*exp(-partable$exp_dep[n_p]*(1:3)))/(exp(partable$age_dep[n_p]*x))})} ) 
partable <- partable %>% mutate(par_id=row_number(),
      const_delta=R0/unlist(lapply(l_delta_susc, function(x) R0_calc_SIRS(C_m,x,rho,n_inf))),seas_conc_lim=0.85,
      npi_start=npi_dates[1],npi_stop=npi_dates[2],seas_start_wk=42,seas_stop_wk=8)
# filtering param sets, from 1st fit: 
# selected parsets are along the line `age=-exp/3+5/6` (and the point (age,exp)=(1/8,1.75))
partable <- partable %>% mutate(age_dep_fit=5/6-exp_dep/3) %>% filter(abs(age_dep-age_dep_fit)/age_dep<1/3) %>%
  select(!age_dep_fit)
write_csv(partable,"partable.csv"); write_csv(partable,"partable.csv")
# agegroup indices for maternal immunity
mat_imm_flag <- TRUE; mat_imm_inds<-list(fun_sub2ind(i_inf=1,j_age=1,"R",c("S","I","R"),n_age,3),
                                         fun_sub2ind(i_inf=c(1,2,3),j_age=9,"R",c("S","I","R"),n_age,3),
                                         fun_sub2ind(i_inf=c(1,2,3),j_age=9,"S",c("S","I","R"),n_age,3))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# hospitalisations data
source("fcns/calc_hosp_rates.R")

ggplot(uk_rsv_hospitalisation_estimate,aes(x=midpoint,y=rate_per_100e3_person_year,color=source)) + geom_point() + 
  geom_line(data=fit_vals_under18,aes(x=age,y=fitvals,color=source)) + ylab("hospitalisations/100e3 population") +
  scale_x_continuous(expand=expansion(0.01,0)) + theme_bw() + scale_y_log10(expand=expansion(0.01,0))
# save
# ggsave(paste0("data/uk_rsv_hospitalisation_fit_5sources_log10.png"),width=32,height=20,units="cm")
# plot rates/100e3 from population-based estimates
ggplot(fit_vals_under18_aver %>% pivot_longer(!c(age,upper_limit,size,perc_pop,average_fit)),
  aes(x=age,y=average_fit))+geom_line()+geom_point(fill=NA,shape=21)+ylab("hospitalisations/100K population/year")+theme_bw() 

# plots both estimates
left_join(hosp_probabilities,hosp_rates,by="agegroup_name") %>% 
  select(c(agegroup_name,hosp_num_from_per_inf_prob,hosp_number_from_pop_estim)) %>% pivot_longer(!agegroup_name) %>%
ggplot(aes(x=agegroup_name, y=value, fill=name)) + labs(fill="") + xlab("Age Group") + ylab("hospitalisations/year") + 
  geom_bar(stat="identity", width=0.5, position = "dodge")  + theme_bw() + standard_theme + theme(legend.position="top") +
  scale_y_continuous(expand=expansion(0,0))
ggsave(paste0("data/uk_rsv_hospitalisation_popul_probperinf_estims.png"),width=32,height=20,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# RUN SIMULATIONS
# write file that'll run scripts
simul_length_yr<-25; n_post_npi_yr<-4; n_core<-64; memory_max <- 8; start_date_dyn_save <- "2018-09-01"
partable_filename <- "simul_output/parscan/parallel/partable_filtered.csv"
# write_csv(partable,file=partable_filename); 
command_print_runs<-paste0(c("Rscript fcns/write_run_file.R",n_core,nrow(read_csv("partable_filtered.csv")),simul_length_yr,
    n_post_npi_yr,partable_filename,"data/estim_attack_rates.csv SAVE sep_qsub_files",start_date_dyn_save,memory_max),collapse=" ")
system(command_print_runs)
# run calculation
system("sh run_all_parallel_scan.sh")
# download results from cluster
# system("scp lshmk17@hpclogin:RSV-model/simul_output/parscan/parallel/results_summ_all.csv simul_output/parscan/parallel/")
# system("scp lshmk17@hpclogin:RSV-model/simul_output/parscan/parallel/results_dyn_all.csv simul_output/parscan/parallel/")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# READ IN RESULTS
peak_week_lims <- c(48,2)
foldername<-"simul_output/parscan/parallel/parsets_1255_filtered/" 
partable <- read_csv(paste0(foldername,"partable_filtered.csv")) # partable.csv
results_summ_all <- read_csv(paste0(foldername,"results_summ_all.csv")) %>%
  mutate(max_incid_week_check=ifelse(max_incid_week>=peak_week_lims[1]|max_incid_week<=peak_week_lims[2],TRUE,FALSE))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# plot attack rates by age group and years
check_crit=11/11; sel_yrs<-2019; n_sel_yr=length(sel_yrs)
# unlist(lapply(1:length(l_delta_susc), function(x) R0_calc_SIRS(C_m,partable$const_delta[x]*l_delta_susc[[x]],rho,n_inf)))
all_sum_inf_epiyear_age_filtered <- left_join(results_summ_all %>% filter(epi_year %in% sel_yrs),
              partable %>% rename(forcing_peak_week=peak_week),
              by=c("par_id","seasforce_peak","R0","exp_dep","age_dep","seasforc_width_wks")) %>% 
  group_by(seasforce_peak,exp_dep,age_dep,seasforc_width_wks,par_id) %>% 
  filter(sum(attack_rate_check)>=round(n_age*n_sel_yr*check_crit) & sum(seas_share_check)>=round(n_age*n_sel_yr*check_crit) )

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT attack rates, seasonal share, peak week
sel_var<-c("attack_rate_perc","seas_share","max_incid_week")
# plot_labels=c("attack rate % age group","seasonal share of infections","peak week")
estim_rates <- estim_attack_rates %>% select(agegroup_name,median_est,min_est,max_est) %>% 
  pivot_longer(!c(agegroup_name)) %>% rename(type=name) %>% mutate(name="attack_rate_perc")
estim_rates <- bind_rows(estim_rates, estim_rates %>% filter(type!="median_est") %>% 
        mutate(name="max_incid_week",value=ifelse(grepl("min",type),3,48)))
color_var<-"exp_dep" # R0 exp_dep age_dep
ggplot(all_sum_inf_epiyear_age_filtered %>% mutate(attack_rate_perc=ifelse(epi_year==2020,NA,attack_rate_perc),
  agegroup_name=factor(agegroup_name,levels=unique(agegroup_name))) %>% ungroup() %>% select(c(par_id,epi_year,
  agegroup_name,attack_rate_perc,seas_share,max_incid_week,exp_dep,age_dep,seasforc_width_wks,R0)) %>% 
    pivot_longer(!c(epi_year,agegroup_name,par_id,exp_dep,age_dep,seasforc_width_wks,R0)) ) +
  geom_hpline(aes(x=age_dep,y=value,color=get(color_var),group=par_id),width=0.1,size=1/2)+#position=position_dodge(width=1)
  facet_grid(name~agegroup_name,scales="free_y") + scale_y_continuous(expand=expansion(0.02,0))+
  scale_color_gradient2(midpoint=median(c(t(unique(all_sum_inf_epiyear_age_filtered[,color_var])))),low="blue",mid="white",high="red") +
  geom_hline(data=estim_rates,aes(yintercept=value),linetype="dashed",size=1/4)+ 
  xlab("age-dependence")+ylab("")+theme(legend.position="top")+theme_bw()+standard_theme+labs(color=color_var)
# save
ggsave(paste0(foldername,"parscan_attack_rates_filtered_",color_var,".png"),width=32,height=20,units="cm")

######
# which parsets selected?
# partable_filtered <- partable %>% filter(par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id))
# write_csv(partable_filtered,"partable_filtered.csv") # paste0(foldername,)
# plot
ggplot(partable %>% mutate(sel_par=TRUE),aes(x=exp_dep,y=age_dep)) + 
  geom_point(aes(color=factor(round(1/omega)),group=omega),size=3/4,position=position_dodge(width=0.3)) + 
  facet_grid(R0~seasforce_peak+seasforc_width_wks,labeller=labeller(seasforc_width_wks=label_both,
        seasforce_peak=label_both,R0=label_both)) +
  # geom_smooth(method="lm",color="black",size=1/2,se=F) + # geom_smooth(method="loess",se=F,size=1/2) + 
  # geom_point(data=partable,aes(x=exp_dep,y=age_dep),color="grey",size=1/2) + 
  scale_x_continuous(breaks=2*(1:8)/8) + theme(axis.text.x=element_text(size=9),
   axis.text.y=element_text(size=9),strip.text=element_text(size=7),legend.position="top") +
  labs(color="waning constant") + theme_bw() + xlab("exposure") + ylab("age") + standard_theme
# selected parameter sets
ggsave(paste0(foldername,"sel_parsets_scatterplot.png"),width=40,height=20,units="cm")

# PCA on parameter sets
# par_pca <- prcomp(partable_filtered %>% select(exp_dep,age_dep,seasforc_width_wks,R0,peak_week,seasforce_peak),
#               center=TRUE,scale.=TRUE)
# # library(devtools);install_github("vqv/ggbiplot"); library(ggbiplot)
# ggbiplot(par_pca,groups=factor(partable_filtered$seasforce_peak),ellipse=TRUE) # ,labels=partable_filtered$par_id

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# check DYNAMICS of SELECTED SIMUL
# results_dyn_all<-read_csv("simul_output/parscan/parallel/parsets_1255_filtered/dyn_parsets_main1021_1098.csv")
parsets_regular_dyn <- read_csv("simul_output/parscan/parallel/parsets_1255_filtered/partable_filtered_reg_dyn.csv")
results_dyn_all <- bind_rows(lapply(list.files(path=foldername,pattern="dyn_parsets*")[61],
  function(x) read_csv(paste0(foldername,x),col_types=cols()) %>% mutate(date=as.Date(start_date_dyn_save)+t-min(t)) %>%
   select(!name) %>% filter(date>=as.Date("2018-10-01") & date<=as.Date("2024-05-01") & par_id %in% parsets_regular_dyn$par_id) ) )

# start_date_dyn_save<-as.Date("2018-10-01")
sel_parsets<-unique(results_dyn_all$par_id) # [2] # unique(all_sum_inf_epiyear_age_filtered$par_id)
ggplot(results_dyn_all %>% filter((par_id %in% sel_parsets) & agegroup<=7) %>% 
   filter(date<as.Date("2023-04-15") & date>as.Date("2018-09-01") &  date<=as.Date("2022-04-01") )) + 
  geom_line(aes(x=date,y=value,color=factor(par_id))) + 
  facet_grid(infection~agegroup,scales="free_y",labeller=labeller(infection=label_both,agegroup=label_both)) +
  # scale_color_brewer(palette = "YlOrRd") + # scale_color_discrete() + 
  geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="grey",alpha=0.01) +
  geom_vline(xintercept=as.Date(paste0(2018:2022,"-12-13"))-56,linetype="dashed",size=1/4) + theme_bw() + 
  geom_vline(xintercept=as.Date(paste0(2018:2022,"-12-13"))+56,linetype="dashed",size=1/4) + standard_theme +
  theme(legend.position="none") + scale_x_date(date_breaks="3 month") + xlab("")+ylab("") + labs(color="# par ID")

### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT dynamics by age groups (one simulation)
n_sel=sel_parsets[1]; start_date<-start_date_dyn_save
sel_weeks <- results_dyn_all %>% filter(par_id==n_sel) %>% mutate(date=start_date+t-min(t),week=week(date),
    year=year(date)) %>% filter(week %in% c(9,41,49)) %>% group_by(year,agegroup,week) %>% filter(date==min(date) & infection==1)
fcn_plot_timecourse_by_agegr(results_dyn_all %>% filter(par_id==n_sel) %>% mutate(date=start_date+t-min(t)) %>%
    filter(t %% 7==0 & agegroup<=9 & date>as.Date("2019-07-01") & date<as.Date("2024-04-01")),
    agegroup_name=rsv_age_groups$agegroup_name,sel_agelim=9,varname="value",npidates=npi_dates,
    date_break_val="2 month",selweeks=sel_weeks,alphaval=0.01,vline_w=c(1/4,1/8))
# # sum of all cases
p<-fcn_plot_timecourse_sum(results_dyn_all %>% filter(par_id==n_sel) %>% mutate(date=start_date+t-min(t)) %>%
      filter(agegroup<=9 & date>as.Date("2019-09-01") & date<as.Date("2023-04-01")) %>% group_by(date,infection) %>%
      summarise(value=sum(value)) %>% mutate(infection=factor(infection)),npi_dates,n_peak_week=50)
p + scale_x_date(date_breaks="2 weeks",expand=expansion(0.01,0))

#####################################################################################
# check if 2018/19 and 2019/20 seasons are the same
yday_start_end<-yday(c(as.Date("2018-10-01"),as.Date("2019-04-01")))
sel_parsets<-unique(results_dyn_all$par_id)[2]
ggplot(results_dyn_all %>% mutate(day_of_year=yday(date),year=year(date)) %>% 
         filter(par_id %in% sel_parsets & # date %in% seq(as.Date("2018-10-01"),as.Date("2019-04-01"),1) & 
      agegroup>=8 & (day_of_year<yday_start_end[2]|day_of_year>yday_start_end[1] )) %>% # 
      mutate(epi_year=ifelse(day_of_year>yday_start_end[1],paste0(year,"_",year+1),paste0(year-1,"_",year)),
  day_of_year=ifelse(day_of_year>yday_start_end[1],day_of_year-yday_start_end[1],day_of_year+365-yday_start_end[1])))+
  geom_line(aes(x=day_of_year,y=value,color=factor(epi_year),linetype=epi_year)) +
  facet_grid(infection~agegroup,scales="free_y",labeller=labeller(infection=label_both,agegroup=label_both))+
  theme_bw()+standard_theme+xlab("")+ylab("")+labs(color="# par ID")+scale_x_continuous(expand=expansion(0.01,0)) 

# calculate diff between 2018/19 and 19/20 season - do this on cluster
# system("Rscript fcns/write_interyear_calc_file.R 2018-09-01 2018-10-01 64 8")
# system("qsub start_batches_calc_interyear.sh")  

# plot relative differences
summ_diff_interyr <- read_csv("simul_output/parscan/parallel/parsets_1255_filtered/summ_diff_interyr_reg_dyn.csv")
summ_diff_interyr <- left_join(summ_diff_interyr %>% mutate(par_id_sort=as.numeric(factor(par_id))),
        rsv_age_groups %>% mutate(agegroup=row_number()) %>% select(c(agegroup,stationary_popul)),by="agegroup") %>%
  mutate(attack_rate=cumul_mean_incid/stationary_popul)
# ggplot(summ_diff_interyr) + geom_point(aes(x=par_id_sort,y=sum_rel_diff),size=1/5) + 
#   facet_grid(infection~agegroup,scales = "free") + scale_y_log10() + geom_hline(yintercept=1/10,color="red") +
#   theme_bw() + standard_theme

signif_inf_types_by_age <- summ_diff_interyr %>% group_by(agegroup,infection) %>% 
  summarise(attack_rate=round(mean(attack_rate,na.rm=T),3)) %>% filter(attack_rate>0.01)

# histogram
ggplot(right_join(summ_diff_interyr,
    signif_inf_types_by_age %>% select(c(agegroup,infection)),by=c("agegroup","infection")),aes(sum_rel_diff)) + 
  stat_ecdf(aes(color=factor(infection),group=infection),geom="step") +
  facet_wrap(~agegroup,nrow=2,labeller=labeller(agegroup=label_both)) + ylab("CDF") +labs(color="# infection") +
  geom_vline(xintercept=1/10,linetype="dashed",size=1/2) + scale_x_log10() + theme_bw() + standard_theme + 
  theme(legend.position="top") + xlab("relative inter-year difference in cumulative incidence")
# save
ggsave(paste0(foldername,"interyear_difference_cumul_incid_reg_dyn.png"),width=35,height=30,units="cm")

# select parameter sets with less than 10% inter-year variation in cumul incid
parsets_regular_dyn <- right_join(summ_diff_interyr,signif_inf_types_by_age %>% select(c(agegroup,infection)),
  by=c("agegroup","infection")) %>% group_by(par_id) %>% summarise(score_reg_dyn=sum(sum_rel_diff<0.2)) %>%
  filter(score_reg_dyn==max(score_reg_dyn))

write_csv(partable %>% filter(par_id %in% parsets_regular_dyn$par_id),paste0(foldername,"partable_filtered_reg_dyn.csv"))

# compare with results_summ_all: it's ~ the same, OK
data.frame(age=unique(results_dyn_all$agegroup),
           results_dyn=(results_dyn_all %>% filter(par_id==33161 & (year(date) %in% c(2019,2020))&(week(date)>=42 | week(date)<=8 ) ) %>%
                          group_by(agegroup) %>% summarise(results_dyn=sum(value)))$results_dyn,
           results_summ=(results_summ_all %>% filter(par_id==33161 & epi_year==2019))$inf_in_seas) %>% pivot_longer(!age) %>%
ggplot(aes(x=age, y=value, fill=name)) + labs(fill="") + xlab("Age Group") + ylab("cases/year") +
  geom_bar(stat="identity",width=0.5,position="dodge") + theme_bw() + standard_theme + theme(legend.position="top") +
  scale_y_continuous(expand=expansion(0,0))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# calculate hospitalisations
# take 2019_2020 season -> calc # cases -> apply hosp probabilities -> reasonable numbers? if not, adjust hospitalisation rates
# check if 2018_2019 vs 2019_2020 season are fine
# for ONE PARAM SET
seas_cases_sum <- results_dyn_all %>% filter(par_id==33161 & (year(date) %in% 2018:2020) & date<as.Date("2020-07-01") ) %>%
  mutate(epi_year=ifelse(week(date)>=42,paste0(year(date),"_",year(date)+1),paste0(year(date)-1,"_",year(date)))) %>%
  filter(week(date)>=42 | week(date)<=14) %>% group_by(agegroup,epi_year) %>% summarise(sum_cases=sum(value)) %>%
  mutate(agegroup_name=rsv_age_groups$agegroup_name[agegroup])
# plot
ggplot(seas_cases_sum,aes(x=agegroup,y=sum_cases/1e6,color=epi_year)) + 
  geom_hpline(width=0.47,size=2,position=position_dodge(width=1)) + geom_vline(xintercept=(0:11)+1/2,linetype="dashed",size=1/2) +
  labs(fill="") + xlab("Age Group") + ylab("million total cases/season") + scale_x_continuous(expand=expansion(0,0),breaks=1:11) +
  scale_y_log10(breaks=round(10^seq(-2,1,by=1/4),3),expand=expansion(0.02,0))+theme_bw()+standard_theme+theme(legend.position="top")

# adjust hosp probs to get numbers close to LIT estimates
robabilities <- hosp_probabilities %>% mutate(prob_hosp_per_infection_adj=prob_hosp_per_infection*c(0.44,0.75,rep(1,7),4.5,4.5))

# compare hospitalisations from SIMULS to those predicted from (median attack rate)*(hosp prob estims from Hodgson)
simul_hosp <- left_join(results_summ_all %>% filter(par_id %in% parsets_regular_dyn$par_id & epi_year==2019),
        hosp_probabilities,by="agegroup_name") %>% 
  mutate(hosp_num_SIMUL_inf_tot=prob_hosp_per_infection_adj*inf_tot,hosp_num_SIMUL_inf_seas=prob_hosp_per_infection_adj*inf_in_seas) %>% 
  ungroup() %>% select(c(agegroup_name,hosp_num_from_per_inf_prob,hosp_num_SIMUL_inf_tot,hosp_num_SIMUL_inf_seas,par_id)) %>% 
  mutate(agegroup_name=factor(agegroup_name,levels=unique(agegroup_name)))  %>% pivot_longer(!c(agegroup_name,par_id)) %>%
  mutate(par_id=ifelse(grepl("per_inf",name) & par_id!=min(par_id),NA,par_id)) %>% filter(!is.na(par_id))
ggplot(simul_hosp %>% mutate(line_size=ifelse(grepl("per_inf",name),1/5,2)),aes(x=agegroup_name,y=ifelse(value>0,value/1e3,NA),color=name)) + 
  geom_hpline(size=1,width=0.32,position=position_dodge(width=1)) + geom_vline(xintercept=(0:11)+1/2,linetype="dashed",size=1/2) + 
  xlab("Age Group")+ylab("thousand hospitalisations in season/epi year") + #scale_x_continuous(expand=expansion(0,0),breaks=1:11)+
  scale_y_log10(breaks=round(10^seq(-2,2,by=1/4),2),expand=expansion(0.02,0)) + 
  labs(fill="") + theme_bw() + standard_theme + theme(legend.position="top")
# save

# how does this compare to attack rates?
attack_rates_simul_LIT <- left_join(results_summ_all %>% filter(par_id %in% parsets_regular_dyn$par_id & epi_year==2019) %>% 
      select(agegroup_name,par_id,inf_tot,inf_in_seas), data.frame(agegroup_name=rsv_age_groups$agegroup_name,
                 cumul_inf_LIT_ESTIM_median=rsv_age_groups$value*estim_attack_rates$median_est/100,
                 cumul_inf_LIT_ESTIM_min=rsv_age_groups$value*estim_attack_rates$min_est/100,
                 cumul_inf_LIT_ESTIM_max=rsv_age_groups$value*estim_attack_rates$max_est/100),by="agegroup_name") %>% 
  mutate(agegroup_name=factor(agegroup_name,levels=unique(agegroup_name))) %>% pivot_longer(!c(agegroup_name,par_id)) %>% 
  mutate(par_id=ifelse(grepl("LIT_ESTIM",name) & par_id!=min(par_id),NA,par_id),
         categ=ifelse(grepl("inf_tot|inf_in_seas",name),"SIMUL","LIT_estim"),
         name=ifelse(grepl("inf_tot|inf_in_seas",name),paste0(name,"_SIMUL"),"LIT_estim")) %>% filter(!is.na(par_id))
ggplot(attack_rates_simul_LIT,aes(x=agegroup_name,y=ifelse(value>0,value/1e3,NA),group=categ,color=name)) + 
  geom_hpline(size=1,width=0.47,position=position_dodge(width=1)) + geom_vline(xintercept=(0:11)+1/2,linetype="dashed",size=1/2) +
  xlab("Age Group") + ylab("thousand cases/season") + # scale_x_continuous(expand=expansion(0,0),breaks=1:11) +
  scale_y_log10(breaks=round(10^seq(-2,4,by=1/4)),expand=expansion(0.02,0)) + labs(color="") + theme_bw() + standard_theme + 
  theme(legend.position="top",axis.text.x=element_text(size=13),axis.text.y=element_text(size=13),legend.text=element_text(size=16))
# SAVE
ggsave(paste0(foldername,"attack_rates_comparison_with_lit.png"),width=25,height=20,units="cm")

# attack rate in newborns too high, can we lower it by reducing delta_susc?
# "Rscript --vanilla fcns/parscan_runner_cmd_line.R 1 1 25 4 partable_filtered_reg_dyn.csv 
# data/estim_attack_rates.csv SAVE 2018-09-01 > simul_output/manual.out"
# dyn_parsets_main2_2 <- read_csv("simul_output/parscan/parallel/dyn_parsets_main2_2.csv") %>%
#   mutate(date=as.Date(start_date_dyn_save)+t-min(t))
# summ_parsets_main2_2 <- dyn_parsets_main2_2 %>% filter((year(date) %in% 2018:2020) & date<as.Date("2020-07-01") ) %>%
#   mutate(epi_year=ifelse(week(date)>=42,paste0(year(date),"_",year(date)+1),paste0(year(date)-1,"_",year(date)))) %>%
#   filter(week(date)>=42 | week(date)<=14) %>% group_by(agegroup,epi_year) %>% summarise(sum_cases=sum(value)) %>%
#   mutate(agegroup_name=rsv_age_groups$agegroup_name[agegroup])
# left_join(summ_parsets_main2_2 %>% filter(epi_year %in% "2018_2019"), data.frame(agegroup_name=rsv_age_groups$agegroup_name,
#         cumul_inf_LIT_ESTIM_median=rsv_age_groups$value*estim_attack_rates$median_est/100,
#         cumul_inf_LIT_ESTIM_min=rsv_age_groups$value*estim_attack_rates$min_est/100,
#         cumul_inf_LIT_ESTIM_max=rsv_age_groups$value*estim_attack_rates$max_est/100),by="agegroup_name") %>% 
#   ungroup() %>% select(!agegroup) %>% pivot_longer(!c(agegroup_name,epi_year)) %>%
#   mutate(categ=ifelse(name %in% "sum_cases","SIMUL","LIT_ESTIM"),agegroup_name=factor(agegroup_name,levels=unique(agegroup_name))) %>%
# ggplot(aes(x=agegroup_name,y=ifelse(value>0,value/1e3,NA),group=categ,color=categ)) + 
#   geom_hpline(size=2,width=0.47,position=position_dodge(width=1)) + geom_vline(xintercept=(0:11)+1/2,linetype="dashed",size=1/2) +
#   xlab("Age Group") + ylab("thousand cases/season") + # scale_x_continuous(expand=expansion(0,0),breaks=1:11) +
#   scale_y_log10(breaks=round(10^seq(-2,4,by=1/4)),expand=expansion(0.02,0)) + labs(color="") + theme_bw() + standard_theme + 
#   theme(legend.position="top",axis.text.x=element_text(size=13),axis.text.y=element_text(size=13),legend.text=element_text(size=16))

# compare adj exp/age_dep rates
# exp_dep <- partable$exp_dep[k_par]; age_dep <- partable$age_dep[k_par]
# exp_dep <- 2; age_dep <- 1/4
# const_delta <- partable$const_delta[k_par]; delta_primary <- const_delta*exp(-exp_dep*(1:3)) # 1.24
# delta_susc <- sapply(1:n_age, function(x) {delta_primary/(exp(age_dep*x))})
# orig_pars <- data.frame(agegroup=factor(rsv_age_groups$agegroup_name,levels=rsv_age_groups$agegroup_name),
#                         t(delta_susc)) %>% mutate(categ="orig")
# ggplot(bind_rows(orig_pars) %>% pivot_longer(!c(agegroup,categ))) + 
#   geom_hpline(aes(x=agegroup,y=value,color=name,linetype=categ),width=1) + theme_bw() + standard_theme + 
#   geom_vline(xintercept = (0:11)+1/2) + scale_y_log10()

# ggplot(dyn_parsets_main2_2 %>% filter(agegroup<=7) %>% filter(date<as.Date("2023-04-15") &
#          date>as.Date("2018-09-01") &  date<=as.Date("2022-04-01") )) + geom_line(aes(x=date,y=value)) +
#   facet_grid(infection~agegroup,scales="free_y",labeller=labeller(infection=label_both,agegroup=label_both)) +
#   # scale_color_brewer(palette = "YlOrRd") + # scale_color_discrete() +
#   geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="grey",alpha=0.01) +
#   geom_vline(xintercept=as.Date(paste0(2018:2022,"-12-13"))-56,linetype="dashed",size=1/4) + theme_bw() +
#   geom_vline(xintercept=as.Date(paste0(2018:2022,"-12-13"))+56,linetype="dashed",size=1/4) + standard_theme +
#   theme(legend.position="none") + scale_x_date(date_breaks="3 month") + xlab("")+ylab("") + labs(color="# par ID")


# call in param scan results DYNAMICS
# foldername: "simul_output/parscan/parallel/parsets_1255_filtered/"
dyn_all_parsets <- bind_rows(lapply(list.files(foldername,pattern="dyn_parsets*"), 
      function(x) read_csv(file=paste0(foldername,x)) %>% filter(par_id %in% parsets_regular_dyn$par_id) %>% 
  mutate(date=t-min(t)+as.Date(start_date_dyn_save)) %>% select(!c(t,name)) %>% 
  filter(date>=as.Date("2018-10-01")&date<=as.Date("2024-04-01")) %>% group_by(agegroup,date,par_id) %>% summarise(value=sum(value)) %>%
  group_by(par_id,agegroup) %>% mutate(value=round(roll_sum(value,n=7,fill=NA,align="right",by=7))) %>% filter(!is.na(value))  )  )

# dyn_parsets_main1_21 <- dyn_parsets_main1_21_daily %>% group_by(par_id,agegroup) %>%
#   mutate(value=round(roll_sum(value,n=7,fill=NA,align="right"))) %>% filter(!is.na(value))
# # check if weekly sum correct, no gaps
# sel_parset <- unique(dyn_all_parsets$par_id)[1:6]
# ggplot(dyn_all_parsets %>% filter((par_id %in% sel_parset) & date<as.Date("2023-05-01")),
#        aes(x=date,y=value,group=par_id,color=factor(par_id))) + geom_line() + facet_wrap(~agegroup,scales="free_y") +
#   theme_bw() + standard_theme + scale_x_date(date_breaks="4 month") + theme(axis.text.x=element_text(size=10)  )
# ggsave(paste0(foldername,"dyn_weekly_sum_test.png"),width=25,height=15,units="cm")

# select parsets that produce regular annual dynamics
partable_regular_dyn <- partable %>% filter(par_id %in% parsets_regular_dyn$par_id)
ggplot(partable_regular_dyn) + 
  geom_point(aes(x=exp_dep,y=age_dep,color=factor(seasforce_peak),shape=factor(seasforce_peak)),size=2,fill=NA) + 
  theme_bw() + standard_theme

# for agedep-expdep, we should reduce this to a single parameter --> PCA
pca_age_exp <- prcomp(partable_regular_dyn %>% select(c(exp_dep,age_dep)))
pred_pca <- data.frame(predict(pca_age_exp, newdata=partable_regular_dyn %>% select(c(exp_dep,age_dep))),
       exp_dep=partable_regular_dyn$exp_dep,age_dep=partable_regular_dyn$age_dep,par_id=partable_regular_dyn$par_id)
ggplot(pred_pca %>% pivot_longer(!c(PC1,PC2,par_id)),aes(x=PC1,y=value,color=name)) +
  geom_point() + geom_smooth(method='lm') + scale_x_continuous(breaks=(-(2*3):(2*2))/4,limits=c(-1,1)) + theme_bw() + standard_theme
ggsave(paste0(foldername,"exp_age_PCA.png"),width=25,height=20,units="cm")

# cumulative infections per epi_year --> add hospitalisations
results_summ_all_hosp <- left_join(left_join(left_join(results_summ_all %>% filter(par_id %in% partable_regular_dyn$par_id), hosp_probabilities %>% 
  select(c(agegroup_name,prob_hosp_per_infection_adj)),by="agegroup_name"),pred_pca %>% select(par_id,PC1),by="par_id"),
  partable_regular_dyn %>% select(par_id,omega) %>% mutate(omega=1/omega) %>% rename(waning=omega), by="par_id") %>% 
  mutate(hosp_tot=prob_hosp_per_infection_adj*inf_tot,
    hosp_seas=prob_hosp_per_infection_adj*inf_in_seas) %>% select(!prob_hosp_per_infection_adj) %>%
  select(!c(median_est,min_est,max_est,median_all_inf,min_est_all_inf,max_est_all_inf,
            attack_rate_check,seas_share_check,max_incid_week_check,final)) %>% relocate(agegroup_name,.after=agegroup) %>%
 mutate(agegroup_broad=c("<1y","1-2y","2-5y","5+y")[findInterval(agegroup,c(2,4,7)+1)+1]) %>% 
  relocate(agegroup_broad,.after=agegroup_name) %>% relocate(c(inf_tot,inf_in_seas),.before=hosp_tot) %>% 
  relocate(par_id,.after=epi_year) %>% relocate(c(exp_dep,age_dep,PC1,seasforce_peak,R0,waning,seasforc_width_wks),.after=par_id)
# plot sums: pre-NPI, NPI+1 (2021/22), NPI+2 (2022/23)
parsets_broad_age_groups <- results_summ_all_hosp %>% group_by(par_id,epi_year,agegroup_broad) %>% 
  summarise(PC1=unique(PC1),seasforce_peak=unique(seasforce_peak),waning=round(unique(waning)),seasforc_width_wks=unique(seasforc_width_wks),
            R0=unique(R0),inf_tot=sum(inf_tot),inf_in_seas=sum(inf_in_seas),hosp_tot=sum(hosp_tot),hosp_seas=sum(hosp_seas),
      seas_share=mean(seas_share),attack_rate_perc=mean(attack_rate_perc),max_incid_week=mean(max_incid_week)) %>% 
  pivot_longer(!c(par_id,epi_year,agegroup_broad,PC1,seasforce_peak,waning,seasforc_width_wks,R0)) %>%
  group_by(par_id,agegroup_broad,name) %>% summarise(epi_year,PC1,seasforce_peak,waning,seasforc_width_wks,R0,value,
      value_norm=ifelse(epi_year==2019,1,ifelse(name %in% c("max_incid_week","seas_share","attack_rate_perc"),
          value-value[epi_year==2019],value/value[epi_year==2019]))) %>% relocate(name,.after=R0)

# summary plot (median, interquartile range)
summ_broad_age_groups <- parsets_broad_age_groups %>% 
  mutate(value_norm=ifelse(name %in% c("seas_share"), value_norm*100,value_norm)) %>% 
  group_by(agegroup_broad,epi_year,name) %>% summarise(mean=mean(value_norm),
      median=median(value_norm),ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
      ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% filter(epi_year>2019)

# plot summary statistics
sel_vars <- c("attack rate","in-season hospitalisations","annual hospitalisations",
             "in-season infections","annual infections","seasonal share of cases") # "season peak (calendar week)",
for (k_plot in 1:length(sel_vars)) {
dodge_val=0.9
ylab_tag <- ifelse(grepl("season peak|seasonal share|attack",sel_vars[k_plot])," (change from 2019 level)"," (normalised by 2019 level)")
p <- ggplot(summ_broad_age_groups %>% mutate(name=case_when(grepl("attack",name) ~ "attack rate",
          grepl("hosp_seas",name) ~ "in-season hospitalisations", grepl("hosp_tot",name) ~ "annual hospitalisations", 
          grepl("inf_in_seas",name) ~ "in-season infections",grepl("inf_tot",name) ~ "annual infections",
          grepl("max_incid_week",name) ~ "season peak (calendar week)",grepl("seas_share",name) ~ "seasonal share of cases")) %>%
    filter(epi_year>2020 & (name %in% sel_vars[k_plot])),aes(x=factor(epi_year),color=factor(agegroup_broad),group=agegroup_broad)) + 
  geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge(width=dodge_val),alpha=0.3,size=12,show.legend=FALSE) +
  geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),position=position_dodge(width=dodge_val),alpha=0.6,size=12) + 
  geom_hpline(aes(y=median),position=position_dodge(width=dodge_val),width=1/6,size=0.8,color="black") + # 
  geom_vline(xintercept=(0:4)+1/2,size=1/2) +
  scale_x_discrete(expand=expansion(0,0)) + xlab("") + ylab(paste0(sel_vars[k_plot],ylab_tag)) + 
  labs(color="age groups") + theme_bw() + standard_theme + theme(strip.text=element_text(size=15),legend.text=element_text(size=15),
      axis.text.x=element_text(size=13),axis.text.y=element_text(size=13),legend.title=element_text(size=15))
if (grepl("season peak|seasonal share|attack",sel_vars[k_plot])) {
  if (grepl("season peak",sel_vars[k_plot])) {break_vals <- (-5:15)*10} else {break_vals <- (-10:10)*10}
  p <- p + geom_hline(yintercept=0,linetype="dashed",size=1/2) + scale_y_continuous(breaks=break_vals) } else {
    p <- p + scale_y_log10() + geom_hline(yintercept=1,linetype="dashed",size=1/2)}
p
# save
sel_var_filename <- gsub("_calendar_week","",gsub("\\(|\\)","",gsub("-|\\s","_",sel_vars[k_plot])))
ggsave(paste0(foldername,"summ_stats_relative_2019_",paste0(sel_var_filename,collapse="_"),".png"),width=25,height=15,units="cm")
print(sel_vars[k_plot])
}

# plot indiv values, relative to 2019 season - as a function of (agedep - expdep)
sel_var <- c("hosp_seas","hosp_tot")[2] # "hosp_seas"
for (sel_par in c("PC1","seasforce_peak","waning","seasforc_width_wks","R0")){
n_par <- length(unlist(unique(parsets_broad_age_groups[,sel_par])))
p <- ggplot(parsets_broad_age_groups %>% filter(epi_year>2020 & name %in% sel_var),
       aes(x=factor(epi_year),y=value_norm,group=get(sel_par))) + facet_wrap(~agegroup_broad,scales="free_y") + 
  geom_hpline(aes(color=get(sel_par)),size=1/3,width=(1/n_par)*0.8,position=position_dodge(width=0.9)) +
  geom_vline(xintercept=(0:4)+1/2,size=1/4) + geom_hline(yintercept=1,linetype="dashed",size=1/2) + 
  ylab(paste0(ifelse(grepl("seas",sel_var),"in-season","annual")," hospitalisations compared to 2019 season")) + 
  scale_x_discrete(expand=expansion(0,0)) + theme_bw() + standard_theme + xlab("") + labs(color=sel_par) +
  theme(legend.position="top",strip.text=element_text(size=15),axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
        legend.text=element_text(size=15),legend.title=element_text(size=16))
if (sel_par %in% "PC1") {
  p <- p + scale_color_gradient2(low="blue",mid="grey",high="red",midpoint=0) + 
    labs(color="exposure-dependent (-1) <--> age-dependent (1)") };  p
# breaks=round(10^seq(1,6,by=1/2)), # + scale_y_log10(expand=expansion(0.02,0)) 
ggsave(paste0(foldername,sel_var,"_parsets_relative_2019_",sel_par,".png"),width=25,height=20,units="cm")
}

# segment age_exp_dep into x values, calculate summary statistics for each param value
summ_broad_age_groups_byvalue <- parsets_broad_age_groups %>% mutate(age_exp_par_bins=findInterval(PC1,seq(-1,1,by=1/5))) %>%
  group_by(age_exp_par_bins) %>% mutate(age_exp_par_bins=round(mean(PC1),1),waning=round(waning)) %>% select(!c(PC1,value)) %>% 
  relocate(age_exp_par_bins,.after=R0) %>%
  rename(varname=name) %>% pivot_longer(!c(par_id,agegroup_broad,epi_year,varname,value_norm)) %>% 
  rename(parname=name,parvalue=value) %>% relocate(c(varname,value_norm),.after=parvalue) %>%
  mutate(value_norm=ifelse(varname %in% c("seas_share"), value_norm*100,value_norm)) %>%
  group_by(agegroup_broad,epi_year,parname,parvalue,varname) %>% summarise(mean=mean(value_norm),median=median(value_norm),
      ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
      ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% filter(epi_year>2019)
#########  
# PLOT
sel_vars <- c("attack rate","annual hospitalisations", # ,"in-season hospitalisations"
              "annual infections","% cases in-season") # "season peak (calendar week)","in-season infections",
sel_pars <- c("age_exp_par_bins","R0","seasforc_width_wks","seasforce_peak","waning")
for (k_plot_var in 1:length(sel_vars)) {
  for (k_plot_par in 1:length(sel_pars)) {
    sel_par <- sel_pars[k_plot_par]; dodge_val=1
    ylab_tag <- ifelse(grepl("season peak|cases in-season|attack",sel_vars[k_plot_var]),
                     " (change from 2019 level)"," (normalised by 2019 level)")
  df_plot <- summ_broad_age_groups_byvalue %>% mutate(varname=case_when(grepl("attack",varname) ~ "attack rate",
        grepl("hosp_tot",varname) ~ "annual hospitalisations", # grepl("hosp_seas",varname) ~ "in-season hospitalisations",
        grepl("inf_tot",varname) ~ "annual infections", # grepl("inf_in_seas",varname) ~ "in-season infections",
        grepl("max_incid_week",varname) ~ "season peak (calendar week)",grepl("seas_share",varname) ~ "% cases in-season")) %>% 
    filter(epi_year>2020 & (varname %in% sel_vars[k_plot_var]) & (parname %in% sel_par)) %>%
    mutate(parname=case_when(grepl("age_exp_par_bins",parname) ~ "exposure (-1) <-> age (1)", 
                             grepl("seasforc_width_wks",parname) ~ "season width (weeks)", 
                             grepl("seasforce_peak",parname) ~ "seasonal forcing (above baseline)",
                             grepl("R0",parname) ~ "R0 (baseline)",grepl("waning",parname) ~ "waning (days)"))
  n_par_value <- length(unique(df_plot$parvalue))
  # colour palette
  if (!grepl("age_exp_par_bins",sel_par)){ colorpal=colorRampPalette(colors=c("orange","red"))(n_par_value)} else  {
    colorpal=colorRampPalette(colors=c("blue","grey","red"))(n_par_value) }
  p <- ggplot(df_plot,aes(x=factor(epi_year),color=factor(parvalue),group=parvalue)) + facet_wrap(~agegroup_broad,scales="free_y") + 
    # geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge(width=dodge_val),alpha=0.3,size=28/n_par_value) + #
    geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),position=position_dodge(width=dodge_val),alpha=0.6,size=24/n_par_value) + #
    geom_hpline(aes(y=median),position=position_dodge(width=dodge_val),width=(1/n_par_value)*0.75,size=0.8,color="black") + 
    geom_vline(xintercept=(0:4)+1/2,size=1/5) + labs(color=unique(df_plot$parname)) + # geom_hline(yintercept=) +
    scale_x_discrete(expand=expansion(0.02,0)) + xlab("") + ylab(paste0(sel_vars[k_plot_var],ylab_tag)) + 
    theme_bw() + standard_theme + theme(strip.text=element_text(size=15),axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=12),legend.text=element_text(size=11),legend.title=element_text(size=12),
        legend.position=ifelse(grepl("expos|forcing",unique(df_plot$parname)),"bottom","right"))+scale_color_manual(values=colorpal)
  if (grepl("season peak|cases in-season|attack",sel_vars[k_plot_var])) { #  # + scale_y_log10() 
    if (grepl("season peak",sel_vars[k_plot_var])) {break_vals <- (-5:15)*10} else {break_vals <- (-10:10)*10}
    p <- p + geom_hline(yintercept=0,linetype="dashed",size=1/2) + scale_y_continuous(breaks=break_vals) } else {
      p <- p + geom_hline(yintercept=1,linetype="dashed",size=1/2)}; p
  # save
  sel_var_filename <- gsub("%","share",gsub("_calendar_week","",gsub("\\(|\\)","",gsub("-|\\s","_",sel_vars[k_plot_var]))))
  subfldr_name <- "median_interquant_by_param_value/CI50_only/"
  ggsave(paste0(foldername,subfldr_name,"summ_stats_relative_2019_",paste0(sel_var_filename,collapse="_"),"_",sel_par,".png"),
         width=25,height=20,units="cm")
  print(paste0(c(sel_vars[k_plot_var],sel_pars[k_plot_par]),collapse=", "))
  }
}
##########################################################
# 2 more metrics needed: maximum hospitalisations + length of season (cases/hosp>n)
norm_seas_length_wk <- round(length(as.Date("2020-10-07"):as.Date("2021-03-03"))/7,1)
hosp_sum_prob_broad_agegr <- hosp_probabilities %>% mutate(
  agegroup_broad=c("<1y","1-2y","2-5y","5+y")[findInterval(factor(agegroup_name,levels=unique(agegroup_name)),c(2,4,7)+1)+1]) %>%
  group_by(agegroup_broad) %>% summarise(hosp_sum=sum(hosp_num_from_per_inf_prob)) %>% 
  mutate(hosp_per_week_season=hosp_sum/norm_seas_length_wk)
# simplify dynamics to broad age groups
dyn_all_parsets_broad_age <- left_join(dyn_all_parsets, hosp_probabilities %>% mutate(agegroup=as.numeric(factor(agegroup_name,levels=unique(agegroup_name)))) %>%
            select(agegroup,prob_hosp_per_infection_adj),by="agegroup") %>% mutate(incid_hosp=value*prob_hosp_per_infection_adj,
    agegroup_broad=c("<1y","1-2y","2-5y","5+y")[findInterval(agegroup,c(2,4,7)+1)+1]) %>% group_by(agegroup_broad,date,par_id) %>%
  summarise(value=sum(value),incid_hosp=sum(incid_hosp))
# summary by year
summ_dyn_max_incid_seas_length <- dyn_all_parsets_broad_age %>% mutate(epi_year=ifelse(week(date)>=42,year(date),year(date)-1)) %>% 
  rename(incid_case=value) %>% pivot_longer(!c(agegroup_broad,date,par_id,epi_year)) %>% group_by(epi_year,par_id,agegroup_broad,name) %>%
  mutate(seas_tot_2018=ifelse(epi_year==2018,sum(value),NA)) %>% group_by(par_id,agegroup_broad,name) %>%
  mutate(seas_tot_2018=min(seas_tot_2018,na.rm=T)) %>% ungroup() %>% mutate(above_baseline=value>seas_tot_2018/52) %>%
  group_by(epi_year,par_id,agegroup_broad,name) %>% summarise(max_value=max(value),sum_value=sum(value),
      seas_length_wk=sum(above_baseline)) %>% filter(epi_year>2017)

# left-join with partable to have parameters as inputs
parsets_max_incid_seas_length <- left_join(left_join(summ_dyn_max_incid_seas_length,pred_pca %>% select(par_id,PC1),by="par_id"),
          partable_regular_dyn %>% select(par_id,omega,seasforc_width_wks,R0,seasforce_peak) %>% mutate(omega=1/omega) %>% 
  rename(waning=omega), by="par_id") %>% relocate(c(name,max_value,sum_value,seas_length_wk),.after=seasforce_peak) %>%
  rename(varname=name)%>% pivot_longer(!c(epi_year,par_id,agegroup_broad,PC1,waning,seasforc_width_wks,R0,seasforce_peak,varname)) %>%
  rename(vartype=name) %>% group_by(par_id,agegroup_broad,PC1,waning,seasforc_width_wks,R0,seasforce_peak,varname,vartype) %>% 
  mutate(value_norm=ifelse(grepl("seas_length_wk",vartype),value-value[epi_year==2018],value/value[epi_year==2018]))

# plot summary statistics
# summary plot (median, interquartile range)
summ_max_incid_seas_length <- parsets_max_incid_seas_length %>% 
  group_by(agegroup_broad,epi_year,varname,vartype) %>% summarise(mean=mean(value_norm),
         median=median(value_norm),ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
         ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2])

# plot summary statistics
sel_vars <- c("incid_case","incid_hosp"); sel_vartypes <- c("max_value","seas_length_wk")
for (k_plot in 1:length(sel_vars)) {
  for (k_plot_type in 1:length(sel_vartypes)) {
  dodge_val=0.9
  df_plot <- summ_max_incid_seas_length %>% 
    filter(epi_year>2020 & (varname %in% sel_vars[k_plot]) & (vartype %in% sel_vartypes[k_plot_type])) %>%
    mutate(varname=case_when(grepl("incid_case",varname) ~ "cases",grepl("incid_hosp",varname) ~ "hospitalisations"),
           vartype=case_when(grepl("max_value",vartype) ~ "peak demand", grepl("seas_length_wk",vartype) ~ "above baseline (week)"))
  ylab_tag <- paste0(paste0(unique(df_plot$varname)," ",unique(df_plot$vartype)),
    ifelse(grepl("above",unique(df_plot$vartype))," (change from 2019 level)"," (normalised by 2019 level)"))
  p <- ggplot(df_plot,aes(x=factor(epi_year),color=factor(agegroup_broad),group=agegroup_broad)) +
    # geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge(width=dodge_val),alpha=0.3,size=16,show.legend=FALSE) +
    geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),position=position_dodge(width=dodge_val),alpha=0.6,size=16) + 
    geom_hpline(aes(y=median),position=position_dodge(width=dodge_val),width=1/6,size=0.8,color="black") + # 
    geom_vline(xintercept=(0:4)+1/2,size=1/2) + scale_x_discrete(expand=expansion(0,0)) + xlab("") + ylab(ylab_tag) + 
    labs(color="age groups") + theme_bw() + standard_theme + theme(strip.text=element_text(size=15),legend.text=element_text(size=15),
                     axis.text.x=element_text(size=13),axis.text.y=element_text(size=13),legend.title=element_text(size=15))
  if (grepl("seas_length_wk",sel_vartypes[k_plot_type])) {
    # if (grepl("season peak",sel_vars[k_plot])) {break_vals <- (-5:15)*10} else {break_vals <- (-10:10)*10}
    p <- p + geom_hline(yintercept=0,linetype="dashed",size=1/2) + scale_y_continuous(breaks=-15:5) } else {
      p <- p + geom_hline(yintercept=1,linetype="dashed",size=1/2)}; p # scale_y_log10() + 
  # save
  sel_var_filename <- gsub("_calendar_week","",gsub("\\(|\\)","",gsub("-|\\s","_",sel_vars[k_plot])))
  ggsave(paste0(foldername,"median_interquant_all_collapsed/CI50/summ_stats_relative_2019_linear_",
                paste0(c(sel_vars[k_plot],sel_vartypes[k_plot_type]),collapse="_"),".png"),width=25,height=20,units="cm")
  print(sel_vars[k_plot])
  }
}

################################################
# plot indiv parsets
for (sel_var in c("incid_case","incid_hosp")){
  for (sel_par in c("PC1","seasforce_peak","waning","seasforc_width_wks","R0")){
    for (sel_vartype in sel_vartypes){
  n_par <- length(unlist(unique(parsets_max_incid_seas_length[,sel_par])))
  df_plot <- parsets_max_incid_seas_length %>% 
    filter(epi_year>2020 & (varname %in% sel_var) & (vartype %in% sel_vartype)) %>%
    mutate(varname=case_when(grepl("incid_case",varname) ~ "cases",grepl("incid_hosp",varname) ~ "hospitalisations"),
           vartype=case_when(grepl("max_value",vartype) ~ "peak demand", grepl("seas_length_wk",vartype) ~ "above baseline (weeks)"))
  ylab_tag <- paste0(paste0(unique(df_plot$varname)," ",unique(df_plot$vartype)),
                     ifelse(grepl("above",unique(df_plot$vartype))," (change from 2019 level)"," (normalised by 2019 level)"))
  p <- ggplot(df_plot,aes(x=factor(epi_year),y=value_norm,group=get(sel_par))) + facet_wrap(~agegroup_broad,scales="free_y") + 
    geom_hpline(aes(color=get(sel_par)),size=1/3,width=(1/n_par)*0.8,position=position_dodge(width=0.9)) +
    geom_vline(xintercept=(0:4)+1/2,size=1/4) + geom_hline(yintercept=1,linetype="dashed",size=1/2) + 
    ylab(ylab_tag) + scale_x_discrete(expand=expansion(0,0)) + theme_bw() + standard_theme + xlab("") + labs(color=sel_par) +
    theme(legend.position="top",strip.text=element_text(size=15),axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
          legend.text=element_text(size=15),legend.title=element_text(size=16))
  if (sel_par %in% "PC1") {
    p <- p + scale_color_gradient2(low="blue",mid="grey",high="red",midpoint=0) + 
      labs(color="exposure-dependent (-1) <--> age-dependent (1)") };  p
  # breaks=round(10^seq(1,6,by=1/2)), # + scale_y_log10(expand=expansion(0.02,0)) 
  plot_filename<-paste0(foldername,"parsets_indiv_by_param/peak_duration/",sel_var,"_",sel_vartype,"_parsets_relative_2019_",sel_par,".png")
  print(plot_filename)
  ggsave(plot_filename,width=25,height=20,units="cm")
  }
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# # calculate total infections per epi year, 2019 vs 2021, GROUPED into 0-2y, 2-5y, 5-15y
# cumul_inf_pre_post_npi <- df_plot %>% filter(epi_year!=2020) %>% group_by(par_id,npi_str,epi_year,agegroup_large) %>% 
#   summarise(sum_inf=sum(value),max_inf=max(value),peak_week=unique(peak_week),resurg_week_20pct=week(min(date[value>=0.2*max(value)])),
#             dep_val_sel=unique(dep_val_sel),dep_val=unique(dep_val),R0=unique(R0),
#     seas_forc_peak=unique(seasforce_peak),dep_type=unique(dep_type)) %>% group_by(par_id,npi_str,agegroup_large) %>% 
#   mutate(sum_inf_norm=sum_inf/sum_inf[epi_year==2019],max_inf_norm=max_inf/max_inf[epi_year==2019],
#       peak_week_shift=peak_week[epi_year==2021]-ifelse(peak_week[epi_year==2019]>40,peak_week[epi_year==2019],peak_week[epi_year==2019]+52),
#       shift_seasonstart_week=ifelse(epi_year>=2021,ifelse(resurg_week_20pct[epi_year==2021]<10,
#       resurg_week_20pct[epi_year==2021]+52,resurg_week_20pct[epi_year==2021])-
#       mean(ifelse(resurg_week_20pct[epi_year<2020]<10,resurg_week_20pct[epi_year<2020]+52,resurg_week_20pct[epi_year<2020])),NA))
# 
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# ## plot sum of infections / shift in peak week / shift in start week by large agegroups
# # geom_point, NPI color coded, averaged over R0 and seas_forc_peak values (OR NOT)
# uniq_R0=(cumul_inf_pre_post_npi %>% group_by(dep_type) %>% summarise(n=length(unique(R0))))$n; dodge_val=1
# colorpal=c(colorRampPalette(colors=c("orange","red"))(uniq_R0[1]),colorRampPalette(colors=c("grey","black"))(uniq_R0[2]))
# sel_var="sum_inf_norm" # sum_inf_norm | peak_week_shift | shift_seasonstart_week
# ggplot(cumul_inf_pre_post_npi %>% mutate(dep_val=dep_val_sel) %>%
#          filter(epi_year==2021) %>% pivot_longer(c(sum_inf_norm,max_inf_norm,peak_week_shift,shift_seasonstart_week)) %>%
#      filter(grepl(sel_var,name)) %>% mutate(value=abs(value)) %>% group_by(agegroup_large,dep_type,dep_val,R0,name) %>% # ,seas_forc_peak
#          summarise(mean_val=mean(value),min_val=min(value),max_val=max(value)) %>%
#      mutate(name=ifelse(grepl("max",name),"PEAK","SUM"),dep_type=ifelse(grepl("age",dep_type),"~AGE","~IMMUNITY"),
#     type_R0=paste0(dep_type,", R0=",R0)) %>% ungroup(), #  %>% mutate(dep_val=as.numeric(factor(dep_val)))
#        aes(x=factor(dep_val),y=mean_val,color=type_R0,group=interaction(type_R0))) + # ,size=seas_forc_peak seas_forc_peak,
#   geom_point(position=position_dodge(width=dodge_val)) + scale_size(range=c(0.7,2)/3) +
#   geom_pointrange(aes(ymin=min_val,ymax=max_val),position=position_dodge(width=dodge_val)) +
#   facet_wrap(~agegroup_large,scales="free",nrow=3) + scale_color_manual(values=colorpal) + # 
#   geom_vline(xintercept=0.5+(0:8),linetype="dashed",size=1/3) + guides(color=guide_legend(ncol=5,byrow=TRUE)) + xlab("dependence on age/exp") + 
#   ylab(ifelse(grepl("sum",sel_var),"ratio of infections 2021 to 2019",paste0("forward shift in ",
#         ifelse(grepl("peak",sel_var),"peak","start")," week"))) + labs(size="seasonal forcing (>baseline)",color="") + theme_bw() + 
#   standard_theme + theme(legend.position="top") + scale_x_discrete(expand=expansion(0.1,0))
# # save
# # ggsave(paste0(foldername,"resurgence_",ifelse(grepl("sum",sel_var),"cumul_inf",sel_var),"_rel_to_2019_large_agegroups.png"),
# #        width=30,height=20,units="cm")
# # ggsave(paste0(foldername,"resurgence_",ifelse(grepl("sum",sel_var),"cumul_inf",sel_var),"_rel_to_2019_large_agegroups_seasfor_sep.png"),
# #         width=30,height=20,units="cm")
# 
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# # shift in mean age by large age groups
# output_files <- list.files(foldername,pattern="df_cases_infs_npi_contactlevel_")
# for (k_npi in 0:4){
# full_output <- read_csv(paste0(foldername,output_files[grepl(paste0("_",k_npi*10,"pct"),output_files)])) %>% 
#   mutate(npi_str=k_npi/10,epi_year=ifelse(date>ymd(paste(year(date),"-07-01")),year(date),year(date)-1),
#   in_out_season=ifelse(week(date)<=9 | week(date)>=41,"in","out")) %>% filter(par_id %in% parsets_filtered$par_id) %>% 
#   group_by(epi_year,par_id,dep_type,dep_val,R0,seasforce_peak,npi_str,agegroup) %>% 
#   summarise(sum_inf=sum(value),max_inf=max(value),peak_week=week(date[value==max(value)]),
#             resurg_week_20pct=week(min(date[value>=0.2*max(value)])) ) # %>% ungroup()
# full_output <- left_join(full_output,rsv_age_groups %>% select(agegroup_name,mean_age_weighted) %>% mutate(agegroup=row_number()),by="agegroup")
# if (k_npi==0) {npi_scan_sum_infs_by_age=full_output} else {npi_scan_sum_infs_by_age=rbind(npi_scan_sum_infs_by_age,full_output)} 
# }
# ### end of read-in loop
# 
# # mean age before/after NPIs
# full_output_mean_age = npi_scan_sum_infs_by_age %>% group_by(dep_type) %>% mutate(dep_val=as.numeric(factor(dep_val))) %>%
#   group_by(epi_year,par_id,dep_type,dep_val,R0,seasforce_peak,npi_str) %>% 
#   summarise(mean_age=sum((sum_inf/sum(sum_inf))*mean_age_weighted),
#             mean_age_under_5=sum((sum_inf[agegroup<=7]/sum(sum_inf[agegroup<=7]))*mean_age_weighted[agegroup<=7]),
#             mean_age_under_15=sum((sum_inf[agegroup<=8]/sum(sum_inf[agegroup<=8]))*mean_age_weighted[agegroup<=8])) %>%
#   group_by(par_id,npi_str) %>% mutate(
#    `shift in average age of infection <5y`=ifelse(epi_year>=2021,mean_age_under_5[epi_year==2021]-mean(mean_age_under_5[epi_year<2020]),NA),
#    `shift in average age of infection <15y`=ifelse(epi_year>=2021,mean_age_under_15[epi_year==2021]-mean(mean_age_under_15[epi_year<2020]),NA))
#    
# # subset data
# df_plot_mean_age_shift <- full_output_mean_age %>% pivot_longer(!c(epi_year,par_id,dep_type,dep_val,R0,seasforce_peak,npi_str)) %>%
#   group_by(name,epi_year,par_id,dep_val,dep_type,R0) %>% # ,seasforce_peak
#   summarise(mean_val=mean(value),min_npi_val=min(value),max_npi_val=max(value)) %>% filter(epi_year==2021 & grepl("shift",name)) %>% 
#   mutate(dep_type=ifelse(grepl("age",dep_type),"~age","~immunity")) %>% ungroup() 
# # plot params
# ratio_diff_flag <- "shift"; dodge_val=1
# title_str=ifelse(ratio_diff_flag=="shift","shift (in months) from 2019 to 2021","(mean age in 2021)/(mean age in 2019)")
# # uniq_R0=(df_plot_mean_age_shift %>% group_by(dep_type) %>% summarise(n=length(unique(R0))))$n; dodge_val=1
# # colorpal=c(colorRampPalette(colors=c("orange","red"))(uniq_R0[1]),colorRampPalette(colors=c("darkgrey","black"))(uniq_R0[2]))
# # PLOT mean age shift by large age groups
# ggplot(df_plot_mean_age_shift %>% filter(grepl("<",name)) %>% 
#     mutate(name=factor(name,levels=c("shift in average age of infection <5y","shift in average age of infection <15y"))) %>% 
#     mutate(type_R0=paste0(dep_type,", R0=",R0)), 
#     aes(x=factor(dep_val),y=mean_val*12,group=interaction(type_R0))) + # ,size=seasforce_peak seasforce_peak,
#   geom_point(aes(color=type_R0),position=position_dodge(width=dodge_val)) + scale_size(range=c(0.7,2)/2.5) +
#   geom_pointrange(aes(color=type_R0,ymin=min_npi_val*12,ymax=max_npi_val*12),position=position_dodge(width=dodge_val)) +
#   facet_wrap(~name,scales="free") + scale_color_manual(values=colorpal) + geom_vline(xintercept=0.5+0:8,linetype="dashed",size=1/3) +
#   xlab("strength of dependence on age/exposure") + ylab(title_str) + labs(color=paste0(""),size="seasonal forcing (>baseline)") + 
#   theme_bw() + standard_theme + scale_x_discrete(expand=expansion(0.12,0)) + guides(color=guide_legend(ncol=3)) + # ,byrow=TRUE
#   theme(legend.position="top",axis.title.x=element_text(size=15),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
#         strip.text=element_text(size=15),legend.text=element_text(size=13),legend.title=element_text(size=14))
# # save
# # ggsave(paste0(foldername,"resurgence_mean_age_",ratio_diff_flag,"_2021_2019.png"),width=43,height=25,units="cm") 
# # ggsave(paste0(foldername,"resurgence_mean_age_",ratio_diff_flag,"_2021_2019_seasfor_aver.png"),width=43,height=25,units="cm") 
# 
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# # resurgence features by YEARLY age groups
# agegrp_yr<-paste0(c("0-1","1-2","2-3","3-4","4-5","5-15",">15"),"y")
# sum_max_ratio_by_yearly_agegr <- npi_scan_sum_infs_by_age %>% group_by(dep_type) %>% mutate(dep_val=as.numeric(factor(dep_val))) %>%
#   mutate(age_yr=agegrp_yr[findInterval(agegroup,c(0,3,5,6,7,8,9))]) %>%
#   group_by(epi_year,par_id,dep_type,dep_val,R0,seasforce_peak,npi_str,age_yr) %>% 
#   summarise(sum_inf=sum(sum_inf),max_inf=sum(max_inf),peak_week=mean(peak_week),resurg_week_20pct=mean(resurg_week_20pct)) %>% 
#   group_by(par_id,npi_str,age_yr) %>% mutate(sum_ratio=ifelse(epi_year>=2021,sum_inf[epi_year==2021]/mean(sum_inf[epi_year<2020]),NA),
#          max_ratio=ifelse(epi_year>=2021,max_inf[epi_year==2021]/mean(max_inf[epi_year<2020]),NA),
#          shift_peak_week=ifelse(epi_year>=2021,
#           ifelse(peak_week[epi_year==2021]<10,peak_week[epi_year==2021]+52,peak_week[epi_year==2021])-
#             mean(ifelse(peak_week[epi_year<2020]<10,peak_week[epi_year<2020]+52,peak_week[epi_year<2020])),NA),
#          shift_seasonstart_week=ifelse(epi_year>=2021,
#                     ifelse(resurg_week_20pct[epi_year==2021]<10,resurg_week_20pct[epi_year==2021]+52,resurg_week_20pct[epi_year==2021])-
#                     mean(ifelse(resurg_week_20pct[epi_year<2020]<10,resurg_week_20pct[epi_year<2020]+52,resurg_week_20pct[epi_year<2020])),NA))
# #
# # plot ratio of CUMUL infections 2021/2019
# # subset data
# df_plot_yr_cum_inf <- sum_max_ratio_by_yearly_agegr %>% filter(epi_year==2021 & !grepl(">15y",age_yr)) %>% 
#   mutate(dep_type=ifelse(grepl("age",dep_type),"~age","~immunity")) %>% group_by(age_yr,dep_type,dep_val,R0) %>% # ,par_id 
#   summarise(mean_val=mean(sum_ratio),min_npi_val=min(sum_ratio),max_npi_val=max(sum_ratio) ) %>% ungroup() %>%
#   mutate(type_R0=factor(paste0(dep_type,", R0=",R0)) )
# # filter(dep_val %in% joint_dep_vals) %>% mutate(dep_val=as.numeric(factor(dep_val))) %>% 
# # plot params
# # uni_dep_vals=(df_plot_yr_cum_inf %>% group_by(dep_type) %>% summarise(n=length(unique(R0))))$n
# # colorpal=c(colorRampPalette(colors=c("orange","red"))(uni_dep_vals[1]),colorRampPalette(colors=c("grey","black"))(uni_dep_vals[2]))
# varname <- "sum_ratio"; dodge_val=1
# title_str <- ifelse(!grepl("ratio",varname),"shift in peak incidence week","ratio of cumulative infections 2021 to 2019")
# # plot
# ggplot(df_plot_yr_cum_inf,aes(x=factor(dep_val),y=mean_val,color=type_R0,group=interaction(seasforce_peak,type_R0),size=seasforce_peak)) + # 
#   geom_point(position=position_dodge(width=dodge_val)) + scale_size(range=c(0.7,2)/3) +
#   geom_pointrange(aes(ymin=min_npi_val,ymax=max_npi_val),position=position_dodge(width=dodge_val)) +
#   facet_wrap(~age_yr,scales="free",nrow=3) + scale_color_manual(values=colorpal) + 
#   geom_vline(xintercept=0.5+(0:8),linetype="dashed",size=1/3) + # guides(color=guide_legend(ncol=5,byrow=TRUE)) + 
#   xlab("dependence on age/immunity") + ylab(title_str) + labs(size="seasonal forcing (+ baseline)",color="") + 
#   theme_bw() + standard_theme  + scale_x_discrete(expand=expansion(0.1,0)) + theme(legend.position="top",
#       axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),
#       strip.text=element_text(size=15),legend.text=element_text(size=13),legend.title=element_text(size=14))
# # save
# # ggsave(paste0(foldername,"resurgence_",varname,"_2021_2019_agegr_yrly_seasfor_averaged.png"),width=35,height=25,units="cm")
# # ggsave(paste0(foldername,"resurgence_",varname,"_2021_2019_agegr_yrly_seasfor_sep.png"),width=45,height=30,units="cm") 
# 
# ### ### ### ### ### 
# # season timing
# # peak
# sel_var<-"shift_seasonstart_week" # !!sym(sel_var)
# df_plot_yr_seas_timing <- sum_max_ratio_by_yearly_agegr %>% filter(epi_year==2021 & !grepl(">15y",age_yr)) %>% 
#   mutate(dep_type=ifelse(grepl("age",dep_type),"~age","~immunity")) %>% ungroup() %>%
#   group_by(age_yr,dep_type,dep_val,R0,seasforce_peak) %>% # 
#   summarise(mean_val=abs(mean(!!sym(sel_var))),min_npi_val=abs(min(!!sym(sel_var))),max_npi_val=abs(max(!!sym(sel_var))) ) %>% ungroup() %>%
#   mutate(type_R0=factor(paste0(dep_type,", R0=",R0)) )
# # seas timing
# ggplot(df_plot_yr_seas_timing,aes(x=factor(dep_val),y=mean_val,color=type_R0,group=interaction(seasforce_peak,type_R0),size=seasforce_peak)) + #
#   geom_point(position=position_dodge(width=dodge_val)) + scale_size(range=c(0.7,2)/3.5) +
#   geom_pointrange(aes(ymin=min_npi_val,ymax=max_npi_val),position=position_dodge(width=dodge_val)) +
#   facet_wrap(~age_yr,scales="free",nrow=3) + scale_color_manual(values=colorpal) + 
#   geom_vline(xintercept=0.5+(0:8),linetype="dashed",size=1/3) + # guides(color=guide_legend(ncol=5,byrow=TRUE)) + 
#   xlab("dependence on age/exp") + ylab(paste0("forward shift in season ",ifelse(grepl("start",sel_var),"start","peak")," (weeks)")) + 
#   labs(size="seasonal forcing (+baseline)",color="") + theme_bw() + standard_theme  + scale_x_discrete(expand=expansion(0.1,0)) +
#   theme(legend.position="top",axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),legend.title=element_text(size=14),
#         axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),strip.text=element_text(size=15),legend.text=element_text(size=13))
# # SAVE
# # ggsave(paste0(foldername,"resurgence_shift_",ifelse(grepl("start",sel_var),"start","peak"),"week_2021_2019_agegr_yrly_seasfor_sep.png"),
# #       width=45,height=30,units="cm")
# #
# # ggsave(paste0(foldername,"resurgence_shift_",ifelse(grepl("start",sel_var),"start","peak"),"week_2021_2019_agegr_yrly_seasfor_aver.png"),
# #      width=40,height=30,units="cm")
# 
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# # plot seas forcing
# list_forcing_vector_npi=list()
# for (npi_red_str in c(0:4)/10){
# g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% fun_shutdown_seasforc(npi_dates,years_pre_post_npi=c(3,3),
#               season_width_wks=seasforc_width_wks,init_mt_day="06-01",ifelse(grepl("exp",parsets_filtered$dep_type[k_par_filt]),45,49),
#                                     forcing_above_baseline=1,npireduc_strength=0)
# npi_inds=as.numeric(npi_dates[1]-simul_start_end[1]):as.numeric(npi_dates[2]-simul_start_end[1])
# forcing_vector_npi[npi_inds]=1 + (forcing_vector_npi[npi_inds]-1)*npi_red_str
# list_forcing_vector_npi[[1+npi_red_str*10]]=forcing_vector_npi }
# # fcn_plot_seas_forc(simul_start_end,forcing_vector_npi,seas_lims_wks=c(7,42),npi_dates,date_resol="3 month")
# xx=as.data.frame(list_forcing_vector_npi); colnames(xx)=paste0("NPI strength=",100*(0:4)/10)
# df_seas_forc_npi=data.frame(date=seq(simul_start_end[1],simul_start_end[2],by=1),xx) %>% pivot_longer(!date) %>% 
#   mutate(year=year(date),week=week(date),name=paste0(gsub("NPI.strength.","",name),"%"))
# # plot
# ggplot(df_seas_forc_npi %>% filter(date<=as.Date("2022-07-01")),aes(x=date,y=value,color=name)) + geom_line(size=1.05) + 
#   geom_vline(data=df_seas_forc_npi %>% filter(week %in% c(7,42) & name=="1") %>% group_by(week,year) %>% filter(date==min(date)),
#              aes(xintercept=date),linetype="dashed",color="black",size=1/4,show.legend=F) + labs(color="contacts during NPIs (% normal)") + 
#   scale_x_date(date_breaks="2 month",expand=expansion(0.01,0)) + theme_bw() + standard_theme + xlab("") + ylab("strength of forcing") +
#   geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="grey",alpha=0.01,show.legend=F,color=NA) + 
#   theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),legend.text=element_text(size=14),legend.position="top",
#         legend.title=element_text(size=16),axis.title.y=element_text(size=15))
# # save
# ggsave(paste0(foldername,"seasonal_forcing_NPI_contactlevel.png"),width=32,height=22,units="cm")

