# This script is for reproducing the results and figures of the manuscript at [...]
# Mihaly Koltai, Nov/2021
####
# clear workspace
rm(list=ls())
# go to folder of file [need library(rstudioapi)!]
# currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# load constant parameters and functions for simulations, specify folder where inputs are stored
source("load_params.R"); foldername <- "repo_data/"
# estimated attack rates
estim_attack_rates <- read_csv(paste0(foldername,"estim_attack_rates.csv"))
# % cases within season (filtering parameter sets)
npi_dates<-as.Date(c("2020-03-26","2021-05-17"))
# set up the table of parameter vectors
partable <- bind_rows(expand.grid(list(exp_dep=seq(1/4,2,1/8),age_dep=seq(1/8,1,1/16),seasforc_width_wks=c(3,5,7),
                                       R0=1+(0:5)/10,peak_week=c(48),seasforce_peak=c(3/4,1,5/4),omega=c(1/250,1/350,1/450)) ) )
# calculating the susceptibility parameter (delta_i_j)
l_delta_susc <- lapply(1:nrow(partable), function(n_p) {sapply(1:n_age,
                                           function(x) {(1*exp(-partable$exp_dep[n_p]*(1:3)))/(exp(partable$age_dep[n_p]*x))})} ) 
partable <- partable %>% mutate(par_id=row_number(),
                                const_delta=R0/unlist(lapply(l_delta_susc, function(x) R0_calc_SIRS(C_m,x,rho,n_inf))),seas_conc_lim=0.85,
                                npi_start=npi_dates[1],npi_stop=npi_dates[2],seas_start_wk=42,seas_stop_wk=8); rm(l_delta_susc)
# filtering param sets: selected parsets are along the line `age=-exp/3+5/6` (and the point (age,exp)=(1/8,1.75))
partable <- partable %>% mutate(age_dep_fit=5/6-exp_dep/3) %>% filter(abs(age_dep-age_dep_fit)/age_dep<1/3) %>% select(!age_dep_fit)
# check the size (>x Mb) of objects you have in the workspace by `fcn_objs_mem_use(min_size=1)`
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# load hospitalisation data: `hosp_probabilities` contains hospitalisation/infection probabilities to convert cases to hospitalisations
source("fcns/calc_hosp_rates.R")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Figure 2: susceptibility as a function of age and exposure
age_exp_dep_uniqvals <- bind_rows(expand.grid(list(exp_dep=seq(1/4,2,1/8),age_dep=seq(1/8,1,1/16),age=1:11,exp=1:3))) %>%
  mutate(suscept_unscaled=exp(-(exp_dep*exp+age_dep*age)))
age_exp_dep_uniqvals <- age_exp_dep_uniqvals %>% 
                mutate(const_delta=1/unlist(lapply(lapply(1:nrow(age_exp_dep_uniqvals), function(n_p) {sapply(1:n_age,
                 function(x) {(1*exp(-age_exp_dep_uniqvals$exp_dep[n_p]*(1:3)))/(exp(age_exp_dep_uniqvals$age_dep[n_p]*x))})} ) , 
                 function(x) R0_calc_SIRS(C_m,x,rho,n_inf))),susc_scaled=suscept_unscaled*const_delta)
ggplot(age_exp_dep_uniqvals %>% filter(exp_dep %in% seq(1/4,2,3/4) & age_dep %in% seq(1/8,1,1/4)) %>%
         rename(`exposure-dependence`=exp_dep,`age-dependence`=age_dep) %>% mutate(age=factor(rsv_age_groups$agegroup_name[age],
      levels=unique(rsv_age_groups$agegroup_name))) )  + geom_line(aes(x=age,color=factor(exp),group=exp,y=susc_scaled),size=1.06) + 
  facet_grid(`exposure-dependence`~`age-dependence`,labeller=labeller(`exposure-dependence`=label_both,`age-dependence`=label_both)) + 
  scale_y_log10() + theme_bw() + standard_theme + labs(color="exposure") + xlab("age group") + ylab(expression(delta[exp]^(age)))
# ggsave
# ggsave(paste0(foldername,"age_exp_dep.png"),width=22,height=18,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# RUN SIMULATIONS with parallelisation (requires multiple cores)
# Go down to RESULTS if you want to read in the results from the parameter sampling already
#
# `start_date_dyn_save` will define the timepoint from which to save simulation outputs
# `simul_length_yr` is the length of simulations (in years), `n_post_npi_yr` is the number of years after the NPIs
# `n_core`: number of cores used, `memory_max`: allocated memory (GB); start_date_dyn_save: date from which to save results
simul_length_yr<-25; n_post_npi_yr<-4; n_core<-64; memory_max <- 8; start_date_dyn_save <- "2018-09-01" 
# this is the 1255 parameters selected based on the criteria of 1) attack rates 2) seasonal concentration of cases
partable_filename <- "repo_data/partable_filtered.csv"; n_row <- nrow(read_csv(partable_filename))
# we will split the parameter table into `n_core` batches and run them in parallel, the sh file will launch the jobs
# write the file that will launch jobs
command_print_runs<-paste0(c("Rscript fcns/write_run_file.R",n_core,n_row,simul_length_yr,n_post_npi_yr,partable_filename,
          "data/estim_attack_rates.csv SAVE sep_qsub_files",start_date_dyn_save,memory_max),collapse=" ")
system(command_print_runs)
# run calculation (this is for multiple cores) by:
# `start_batches.sh` in the folder `batch_run_files/`
####################################
# Calculations work by the file "fcns/parscan_runner_cmd_line.R" receiving command line arguments from the parameter table and calling the 
# function (sirs_seasonal_forc_mat_immun) that generates & solves ODEs. `sirs_seasonal_forc_mat_immun` is in "fcns/RSV_model_functions.R". 
# The force of infection term is partly constructed outside this function and passed as inputs, most importantly the input 
# `contmatr_rowvector` is the contact matrix stacked 3x times (3 infections) and normalised by population sizes to generate the 
# force of infection terms normalised by age group sizes. See the article's SI Methods
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# RESULTS
# results from the parameter sampling after filtering already stored in repo_data/
# peak_week_lims <- c(46,5)
results_summ_all <- read_csv(paste0(foldername,"results_summ_all.csv")) 
# %>% mutate(max_incid_week_check=ifelse(max_incid_week>=peak_week_lims[1]|max_incid_week<=peak_week_lims[2],TRUE,FALSE))
# filter for parameter sets where 10/11 (or 11/11) age groups satisfy criteria for attack rates and seasonal concentration
check_crit=11/11; sel_yrs<-2019; n_sel_yr=length(sel_yrs)
all_sum_inf_epiyear_age_filtered <- left_join(results_summ_all %>% filter(epi_year %in% sel_yrs),
      partable %>% rename(forcing_peak_week=peak_week),by=c("par_id","seasforce_peak","R0","exp_dep","age_dep","seasforc_width_wks")) %>% 
  group_by(seasforce_peak,exp_dep,age_dep,seasforc_width_wks,par_id) %>% 
  filter(sum(attack_rate_check)>=round(n_age*n_sel_yr*check_crit) & sum(seas_share_check)>=round(n_age*n_sel_yr*check_crit) 
    # & sum(max_incid_week_check)>=round(n_age*n_sel_yr*check_crit) # criteria for peak incidence week (not important as seas conc is ~same)
         ) %>% select(!c(median_est,min_est,max_est,median_all_inf,min_est_all_inf,max_est_all_inf,attack_rate_check,seas_share_check,
                         seas_conc_lim,npi_start,npi_stop,seas_start_wk,seas_stop_wk,final))
# this leads to 1100 (11/11) or 1153 (10/11) parameter sets: 
length(unique(all_sum_inf_epiyear_age_filtered$par_id))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT attack rates, seasonal share, peak week (this is not a figure in manuscript)
sel_var<-c("attack_rate_perc","seas_share","max_incid_week")
estim_rates <- estim_attack_rates %>% select(agegroup_name,median_est,min_est,max_est) %>% 
  pivot_longer(!c(agegroup_name)) %>% rename(type=name) %>% mutate(name="attack_rate_perc")
estim_rates <- bind_rows(estim_rates, estim_rates %>% filter(type!="median_est") %>% 
                           mutate(name="max_incid_week",value=ifelse(grepl("min",type),3,48)))
color_var<-"exp_dep" # R0 exp_dep age_dep
ggplot(all_sum_inf_epiyear_age_filtered %>% mutate(attack_rate_perc=ifelse(epi_year==2020,NA,attack_rate_perc),
          agegroup_name=factor(agegroup_name,levels=unique(agegroup_name))) %>% ungroup() %>% select(c(par_id,epi_year,
          agegroup_name,attack_rate_perc,seas_share,max_incid_week,exp_dep,age_dep,seasforc_width_wks,R0)) %>% 
         pivot_longer(!c(epi_year,agegroup_name,par_id,exp_dep,age_dep,seasforc_width_wks,R0)) ) +
  geom_hpline(aes(x=age_dep,y=value,color=get(color_var),group=par_id),width=0.1,size=1/2) +
  facet_grid(name~agegroup_name,scales="free_y") + scale_y_continuous(expand=expansion(0.02,0))+
  scale_color_gradient2(midpoint=median(c(t(unique(all_sum_inf_epiyear_age_filtered[,color_var])))),low="blue",mid="white",high="red") +
  geom_hline(data=estim_rates,aes(yintercept=value),linetype="dashed",size=1/4)+ 
  xlab("age-dependence")+ylab("")+theme(legend.position="top")+theme_bw()+standard_theme+labs(color=color_var)
# save
# ggsave(paste0(foldername,"parscan_attack_rates_filtered_",color_var,".png"),width=32,height=20,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# To remove model parameterisations that exhibit irregular patterns (varying from one year to another),
# we need to calculate relative difference (see SI Methods) between 2018/19 and 19/20 season
# write script for calculation:
# system("Rscript fcns/write_interyear_calc_file.R 2018-09-01 2018-10-01 64 8")
# run calculation:
# system("qsub start_batches_calc_interyear.sh")  

# Results already in the repo_data folder
summ_diff_interyr <- left_join(read_csv(paste0(foldername,"summ_diff_interyr_reg_dyn.csv")) %>% 
      mutate(par_id_sort=as.numeric(factor(par_id))),rsv_age_groups %>% mutate(agegroup=row_number()) %>% 
        select(c(agegroup,stationary_popul)),by="agegroup") %>% mutate(attack_rate=cumul_mean_incid/stationary_popul)

# SI Figure 3: CDF of interyear difference in incidence 
# (only taking those infection types with at least 1% attack rate. Eg. 1st infections for adults are negligible, so inter-year differences 
# in these small numbers we don't take into account)
ggplot(right_join(summ_diff_interyr,
                  summ_diff_interyr %>% group_by(agegroup,infection) %>% summarise(attack_rate=round(mean(attack_rate,na.rm=T),3)) %>% 
                    filter(attack_rate>0.01) %>% select(c(agegroup,infection)),by=c("agegroup","infection")),
       aes(sum_rel_diff)) + stat_ecdf(aes(color=factor(infection),group=infection),geom="step") +
  facet_wrap(~agegroup,nrow=2,labeller=labeller(agegroup=label_both))+ylab("CDF")+labs(color="# infection")+scale_x_log10() +
  geom_vline(xintercept=1/10,linetype="dashed",size=1/2) + theme_bw() + standard_theme + xlab("relative difference in incidence") +
  theme(legend.position="top",strip.text=element_text(size=14),legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),axis.title.x=element_text(size=16))
# save
ggsave(paste0(foldername,"interyear_difference_cumul_incid_reg_dyn.png"),width=30,height=25,units="cm")
##############################################################
# select the parameter sets with less than x% inter-year variation in cumul incid
parsets_regular_dyn <- right_join(summ_diff_interyr,
  summ_diff_interyr %>% group_by(agegroup,infection) %>% summarise(attack_rate=round(mean(attack_rate,na.rm=T),3)) %>% 
  filter(attack_rate>0.01) %>% select(c(agegroup,infection)),by=c("agegroup","infection")) %>% group_by(par_id) %>% 
  summarise(score_reg_dyn=sum(sum_rel_diff<0.2)) %>% filter(score_reg_dyn==max(score_reg_dyn))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# compare hospitalisations from SIMULATIONS to those predicted from (median attack rate in LITERATURE)*(hosp prob estims from Hodgson)
simul_hosp <- left_join(results_summ_all %>% filter(par_id %in% parsets_regular_dyn$par_id & epi_year==2019),
                        hosp_probabilities,by="agegroup_name") %>% 
  mutate(hosp_num_SIMUL_inf_tot=prob_hosp_per_infection_adj*inf_tot,hosp_num_SIMUL_inf_seas=prob_hosp_per_infection_adj*inf_in_seas) %>% 
  ungroup() %>% select(c(agegroup_name,hosp_num_from_per_inf_prob,hosp_num_SIMUL_inf_tot,hosp_num_SIMUL_inf_seas,par_id)) %>% 
  mutate(agegroup_name=factor(agegroup_name,levels=unique(agegroup_name)))  %>% pivot_longer(!c(agegroup_name,par_id)) %>%
  mutate(par_id=ifelse(grepl("per_inf",name) & par_id!=min(par_id),NA,par_id)) %>% filter(!is.na(par_id))
# plot number of hospitalisations from simulations compared to number expected from (literature estimate)*(hospit probab per infection)
ggplot(simul_hosp %>% mutate(line_size=ifelse(grepl("per_inf",name),1/5,2)),
       aes(x=agegroup_name,y=ifelse(value>0,value/1e3,NA),color=name)) + 
  geom_hpline(size=1,width=0.32,position=position_dodge(width=1)) + geom_vline(xintercept=(0:11)+1/2,linetype="dashed",size=1/2) + 
  xlab("Age Group")+ylab("thousand hospitalisations in season/epi year") + #scale_x_continuous(expand=expansion(0,0),breaks=1:11)+
  scale_y_log10(breaks=round(10^seq(-2,2,by=1/4),2),expand=expansion(0.02,0)) + labs(fill="",color="") + 
  theme_bw() + standard_theme + theme(legend.position="top")

# how do simulated attack rates compare to LIT estimates?
attack_rates_simul_LIT <- left_join(results_summ_all %>% filter(par_id %in% parsets_regular_dyn$par_id & epi_year==2019) %>% 
          select(agegroup_name,par_id,inf_tot,inf_in_seas), data.frame(agegroup_name=rsv_age_groups$agegroup_name,
          cumul_inf_LIT_ESTIM_median=rsv_age_groups$value*estim_attack_rates$median_est/100,
          cumul_inf_LIT_ESTIM_min=rsv_age_groups$value*estim_attack_rates$min_est/100,
          cumul_inf_LIT_ESTIM_max=rsv_age_groups$value*estim_attack_rates$max_est/100),by="agegroup_name") %>% 
  mutate(agegroup_name=factor(agegroup_name,levels=unique(agegroup_name))) %>% pivot_longer(!c(agegroup_name,par_id)) %>% 
  mutate(par_id=ifelse(grepl("LIT_ESTIM",name) & par_id!=min(par_id),NA,par_id),
         categ=ifelse(grepl("inf_tot|inf_in_seas",name),"SIMUL","LIT_estim"),
         name=ifelse(grepl("inf_tot|inf_in_seas",name),paste0(name,"_SIMUL"),"LIT_estim")) %>% filter(!is.na(par_id)) %>%
  mutate(name=case_when(name %in% "inf_in_seas_SIMUL" ~ "in-season cumulative incidence (SIMULATION)",
                        name %in% "LIT_estim" ~ "in-season cumulative incidence (LITERATURE ESTIMATE)"))
simul_median_cumul_inf <- attack_rates_simul_LIT %>% filter(name %in% "in-season cumulative incidence (SIMULATION)") %>% 
  group_by(agegroup_name) %>% summarise(value=median(value))
# SI FIGURE 2
ggplot(attack_rates_simul_LIT %>% filter(name %in% "in-season cumulative incidence (SIMULATION)"), 
       aes(x=agegroup_name,y=ifelse(value>0,value/1e3,NA),color=name)) + 
  geom_hpline(size=1/5,width=0.47,color="red") + # ,position=position_dodge(width=1)
  geom_hpline(data=attack_rates_simul_LIT %>% filter(grepl("LIT",name)) %>% group_by(agegroup_name) %>% 
                arrange(value) %>% mutate(min_med_max=row_number()),aes(linetype=ifelse(min_med_max==2,"solid","dashed")),
              size=1,width=1,color="black") + # position=position_dodge(width=1),
  geom_hpline(data=simul_median_cumul_inf,size=1.5,width=1,position=position_dodge(width=1),color="orange") +
  geom_vline(xintercept=(0:11)+1/2,linetype="dashed",size=1/2) + # scale_color_manual(values=c("black","red")) +
  xlab("Age Group") + ylab("thousand cases/season") + scale_y_log10(breaks=round(10^seq(-2,4,by=1/4)),expand=expansion(0.02,0)) + 
  labs(color="") + theme_bw() + standard_theme + theme(legend.position="none",axis.text.x=element_text(size=13),
                                                       axis.text.y=element_text(size=13),legend.text=element_text(size=16))
# SAVE
# ggsave(paste0(foldername,"attack_rates_comparison_with_lit_good_label.png"),width=25,height=20,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# reduce the two parameters exp_dep and age_dep to their value along their 1st principal component
partable_regular_dyn <- left_join(parsets_regular_dyn %>% select(!score_reg_dyn),partable,by="par_id")
pred_pca <- data.frame(predict(prcomp(partable_regular_dyn %>% select(c(exp_dep,age_dep))), newdata=partable_regular_dyn %>% 
      select(c(exp_dep,age_dep))),K_exp=partable_regular_dyn$exp_dep,K_age=partable_regular_dyn$age_dep,
      par_id=partable_regular_dyn$par_id)
# SI Figure 4 (linear relationship between K_age and K_exp for selected param sets)
ggplot(pred_pca %>% pivot_longer(!c(PC1,PC2,par_id)),aes(x=PC1,y=value)) + geom_point(aes(color=name)) + 
  geom_smooth(aes(group=name,color=name),fill=NA,method='lm') + scale_x_continuous(breaks=(-(2*3):(2*2))/4,limits=c(-1,1)) + theme_bw() + 
  standard_theme + ylab(expression(paste(kappa[exp],", ",kappa[age]))) + labs(color="",fill="") + theme(axis.title.y=element_text(size=14))
# ggsave(paste0(foldername,"exp_age_PCA.png"),width=25,height=20,units="cm")
# plot suscept ~ (age,exposure) for SELECTED parsets
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
suscept_sel_parsets <- left_join(pred_pca %>% rename(exp_dep=K_exp,age_dep=K_age),age_exp_dep_uniqvals,by=c("exp_dep","age_dep")) %>%
  rename(`exposure-dependence`=exp_dep,`age-dependence`=age_dep) %>% mutate(age=factor(rsv_age_groups$agegroup_name[age],
            levels=unique(rsv_age_groups$agegroup_name)),PC1_grouped=findInterval(PC1,seq(-1,1,by=1/5)) ) %>%
  group_by(PC1_grouped) %>% mutate(PC1_grouped=round(mean(PC1),1)) %>% group_by(PC1_grouped,age,exp) %>% 
  summarise(mean_val=mean(susc_scaled),med_val=median(susc_scaled),ci50_low=quantile(susc_scaled,c(0.25,0.75))[1],
            ci50_up=quantile(susc_scaled,c(0.25,0.75))[2],ci95_low=quantile(susc_scaled,c(0.025,0.975))[1],
            ci95_up=quantile(susc_scaled,c(0.025,0.975))[2]) %>% rename(exposure=exp)
colorpal <- colorRampPalette(colors=c("blue","grey","red"))(length(unique(suscept_sel_parsets$PC1_grouped)))
# facet by exposure level
ggplot(suscept_sel_parsets) + geom_line(aes(x=age,color=factor(PC1_grouped),group=PC1_grouped,y=med_val),size=1.06) + 
  facet_wrap(~exposure,labeller=labeller(exposure=label_both))+labs(color="exposure (-1) <-> age (1)") + scale_color_manual(values=colorpal)+
  scale_y_log10() + theme_bw() + standard_theme + theme(strip.text=element_text(size=16),legend.title=element_text(size=16),
      legend.text=element_text(size=15),legend.position="top",axis.text.x=element_text(size=14),axis.text.y=element_text(size=14)) + 
  xlab("age group") + ylab(expression(delta[exp]^(age)))
# save
ggsave(paste0(foldername,"suscept_by_exp_level.png"),width=25,height=20,units="cm")
# facet by 'dependence' parameter: FIGURE 2
label_parseall <- function(variable, value) {plyr::llply(value, function(x) parse(text = paste(variable, x, sep = "==")))}
ggplot(suscept_sel_parsets %>% rename(kappa=PC1_grouped) %>% ungroup() %>% mutate(kappa_num=as.numeric(factor(as.character(kappa)))) %>% 
    filter(kappa_num %in% c(1,3,5,7,9,10)),aes(x=age,color=factor(exposure),group=exposure,fill=factor(exposure))) +
  geom_line(aes(y=med_val),size=1.06) + geom_ribbon(aes(ymin=ci50_low,ymax=ci50_up),color=NA,alpha=0.2) +
  facet_wrap(~kappa,labeller=label_parseall)+labs(color="exposure",fill="exposure") + # scale_color_manual(values=colorpal) +
  scale_y_log10() + theme_bw() + standard_theme + theme(strip.text=element_text(size=18),legend.title=element_text(size=16),
      legend.text=element_text(size=15),legend.position="top",axis.text.x=element_text(size=14),axis.text.y=element_text(size=15),
      axis.title.x=element_text(size=16),axis.title.y=element_text(size=18)) + xlab("age group") + ylab(expression(delta[exp]^(age)))
# save
ggsave(paste0(foldername,"suscept_by_dep_level.png"),width=25,height=20,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Analyse results by epidemiologic parameters
# cumulative infections per epi_year --> add hospitalisations
results_summ_all_hosp <- left_join(left_join(left_join(results_summ_all %>% filter(par_id %in% partable_regular_dyn$par_id), 
                hosp_probabilities %>% select(c(agegroup_name,prob_hosp_per_infection_adj)),by="agegroup_name"),
                pred_pca %>% select(par_id,PC1),by="par_id"),
            partable_regular_dyn %>% select(par_id,omega) %>% mutate(omega=1/omega) %>% rename(waning=omega), by="par_id") %>% 
  mutate(hosp_tot=prob_hosp_per_infection_adj*inf_tot,
         hosp_seas=prob_hosp_per_infection_adj*inf_in_seas) %>% select(!prob_hosp_per_infection_adj) %>%
  select(!c(median_est,min_est,max_est,median_all_inf,min_est_all_inf,max_est_all_inf,
            attack_rate_check,seas_share_check,final)) %>% relocate(agegroup_name,.after=agegroup) %>%
  mutate(agegroup_broad=c("<1y","1-2y","2-5y","5+y")[findInterval(agegroup,c(2,4,7)+1)+1]) %>% 
  relocate(agegroup_broad,.after=agegroup_name) %>% relocate(c(inf_tot,inf_in_seas),.before=hosp_tot) %>% 
  relocate(par_id,.after=epi_year) %>% relocate(c(exp_dep,age_dep,PC1,seasforce_peak,R0,waning,seasforc_width_wks),.after=par_id)
# write_csv(results_summ_all_hosp,paste0(foldername,"results_summ_all_hosp.csv"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# outcomes by individual parsets # plot sums: pre-NPI, NPI+1 (2021/22), NPI+2 (2022/23)
parsets_broad_age_groups <- results_summ_all_hosp %>% group_by(par_id,epi_year,agegroup_broad) %>% 
  summarise(PC1=unique(PC1),seasforce_peak=unique(seasforce_peak),waning=round(unique(waning)),seasforc_width_wks=unique(seasforc_width_wks),
            R0=unique(R0),inf_tot=sum(inf_tot),inf_in_seas=sum(inf_in_seas),hosp_tot=sum(hosp_tot),hosp_seas=sum(hosp_seas),
            seas_share=mean(seas_share),attack_rate_perc=mean(attack_rate_perc),max_incid_week=mean(max_incid_week)) %>% 
  pivot_longer(!c(par_id,epi_year,agegroup_broad,PC1,seasforce_peak,waning,seasforc_width_wks,R0)) %>%
  group_by(par_id,agegroup_broad,name) %>% summarise(epi_year,PC1,seasforce_peak,waning,seasforc_width_wks,R0,value,
                              value_norm=ifelse(epi_year==2019,1,ifelse(name %in% c("max_incid_week","seas_share","attack_rate_perc"),
                              value-value[epi_year==2019],value/value[epi_year==2019]))) %>% relocate(name,.after=R0) %>%
  mutate(age_exp_par_bins=findInterval(PC1,seq(-1,1,by=1/5))) %>% group_by(age_exp_par_bins) %>% 
  mutate(age_exp_par_bins=round(mean(PC1),1)) %>% relocate(age_exp_par_bins,.after=PC1)
# write_csv(parsets_broad_age_groups,paste0(foldername,"parsets_broad_age_groups.csv"))
# parsets_broad_age_groups <- read_csv(paste0(foldername,"parsets_broad_age_groups.csv"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# summary plot (median, interquartile range)
summ_broad_age_groups <- parsets_broad_age_groups %>% 
  mutate(value_norm=ifelse(name %in% c("seas_share"), value_norm*100,value_norm)) %>% 
  group_by(agegroup_broad,epi_year,name) %>% summarise(mean=mean(value_norm),median=median(value_norm),
          ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
          ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% filter(epi_year>2019)
# write_csv(summ_broad_age_groups,paste0(foldername,"summ_broad_age_groups.csv"))
# summ_broad_age_groups <- read_csv(paste0(foldername,"summ_broad_age_groups.csv"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# segment age_exp_dep into x values, calculate summary statistics for each param value
summ_broad_age_groups_byvalue <- parsets_broad_age_groups %>% select(!c(PC1,value)) %>% mutate(waning=round(waning)) %>%
  relocate(age_exp_par_bins,.after=R0) %>% rename(varname=name) %>% pivot_longer(!c(par_id,agegroup_broad,epi_year,varname,value_norm)) %>% 
  rename(parname=name,parvalue=value) %>% relocate(c(varname,value_norm),.after=parvalue) %>%
  mutate(value_norm=ifelse(varname %in% c("seas_share"), value_norm*100,value_norm)) %>%
  group_by(agegroup_broad,epi_year,parname,parvalue,varname) %>% summarise(mean=mean(value_norm),median=median(value_norm),
          ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
          ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% filter(epi_year>2019)
# write_csv(summ_broad_age_groups_byvalue,paste0(foldername,"summ_broad_age_groups_byvalue.csv"))
# summ_broad_age_groups_byvalue <- read_csv(paste0(foldername,"summ_broad_age_groups_byvalue.csv"))
##############################################################
# PLOTS of summary statistics across ALL parameter values
# the plots below include: Figure 3,4,5 and SI Figure 5
sel_vars <- c("attack rate","in-season hospitalisations","cumulative hospitalisations",
              "in-season infections","cumulative infections","% cases in-season")
for (k_plot in 1:length(sel_vars)) {
  for (k_age_excl in 1:2) {
    dodge_val=0.9
    ylab_tag <- ifelse(grepl("season peak|% cases|attack",sel_vars[k_plot])," (change from 2019 level)"," (relative to pre-NPI)")
    df_plot <- summ_broad_age_groups %>% mutate(name=case_when(grepl("attack",name) ~ "attack rate",
                  grepl("hosp_seas",name) ~ "in-season hospitalisations", grepl("hosp_tot",name) ~ "cumulative hospitalisations", 
                  grepl("inf_in_seas",name) ~ "in-season infections",grepl("inf_tot",name) ~ "cumulative infections",
                  grepl("max_incid_week",name) ~ "season peak (calendar week)",grepl("seas_share",name) ~ "% cases in-season")) %>%
      filter(epi_year>2020 & (name %in% sel_vars[k_plot])); 
    if (k_age_excl>1) {df_plot <- df_plot %>% filter(!agegroup_broad %in% "5+y"); size_adj<-1.33} else {size_adj<-1}
    p <- ggplot(df_plot, aes(x=factor(epi_year),color=factor(agegroup_broad),group=agegroup_broad)) + 
      geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge(width=dodge_val),alpha=0.3,size=12*size_adj,show.legend=FALSE) +
      geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),position=position_dodge(width=dodge_val),alpha=0.6,size=12*size_adj) + 
      geom_hpline(aes(y=median),position=position_dodge(width=dodge_val),width=size_adj/6,size=0.8,color="black") + # 
      geom_vline(xintercept=(0:4)+1/2,size=1/2) + scale_x_discrete(expand=expansion(0,0)) +xlab("")+ylab(paste0(sel_vars[k_plot],ylab_tag)) +
      labs(color="age groups") + theme_bw() + standard_theme + theme(strip.text=element_text(size=15),legend.text=element_text(size=15),
      axis.text.x=element_text(size=13),axis.text.y=element_text(size=13),
        legend.title=element_text(size=15),panel.grid.major.x=element_blank())
    if (grepl("season peak|% cases|attack",sel_vars[k_plot])) {
      if (grepl("season peak",sel_vars[k_plot])) {break_vals <- (-5:15)*10} else {break_vals <- (-10:10)*10}
      p <- p + geom_hline(yintercept=0,linetype="dashed",size=1/2) + scale_y_continuous(breaks=break_vals) } else {
        p <- p + scale_y_log10(breaks=c(0.1,1/4,1/2,3/4,1,5/4,3/2,7/4,2,2.5,3,5)) + geom_hline(yintercept=1,linetype="dashed",size=1/2)}
    p
    # save # round(10^seq(-1,1,by=1/8),2)
    sel_var_filename <- gsub("%","share",gsub("_calendar_week","",gsub("\\(|\\)","",gsub("-|\\s","_",sel_vars[k_plot]))))
    subfldr_name<-"median_interquant_all_collapsed/"
    if (!dir.exists(paste0(foldername,subfldr_name))) {dir.create(paste0(foldername,subfldr_name))}
    if (!dir.exists(paste0(foldername,subfldr_name,"until_5y/"))) {dir.create(paste0(foldername,subfldr_name,"until_5y/"))}
    ggsave(paste0(foldername,subfldr_name,ifelse(k_age_excl>1,"until_5y/",""),"summ_stats_relative_2019_",
                  paste0(sel_var_filename,collapse="_"),".png"),width=25,height=20,units="cm")
    print(sel_vars[k_plot])
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT summary stats disaggregated by parameter VALUES
# the figures produced by this loop include: Figure 6-8 (main text) and SI Fig 6-11
sel_vars <- c("attack rate","cumulative hospitalisations", "cumulative infections","% cases in-season")
# "season peak (calendar week)","in-season infections",
sel_pars <- c("age_exp_par_bins","R0","seasforc_width_wks","seasforce_peak","waning")
for (k_plot_var in 1:length(sel_vars)) {
  for (k_plot_par in 1:length(sel_pars)) {
    sel_par <- sel_pars[k_plot_par]; dodge_val=1
    ylab_tag <- ifelse(grepl("season peak|cases in-season|attack",sel_vars[k_plot_var]),
                       " (change from 2019 level)"," (relative to pre-NPI)")
    df_plot <- summ_broad_age_groups_byvalue %>% mutate(varname=case_when(grepl("attack",varname) ~ "attack rate",
                     grepl("hosp_tot",varname) ~ "cumulative hospitalisations", # grepl("hosp_seas",varname) ~ "in-season hospitalisations",
                     grepl("inf_tot",varname) ~ "cumulative infections", # grepl("inf_in_seas",varname) ~ "in-season infections",
                     grepl("max_incid_week",varname) ~ "season peak (calendar week)",grepl("seas_share",varname) ~ "% cases in-season")) %>% 
      filter(epi_year>2020 & (varname %in% sel_vars[k_plot_var]) & (parname %in% sel_par) & (!agegroup_broad %in% "5+y") ) %>%
      mutate(parname=case_when(grepl("age_exp_par_bins",parname) ~ "exposure (-1) <-> age (1)", 
                               grepl("seasforc_width_wks",parname) ~ "season width (weeks)", 
                               grepl("seasforce_peak",parname) ~ "seasonal forcing (above baseline)",
                               grepl("R0",parname) ~ "R0 (baseline)",grepl("waning",parname) ~ "waning (days)"))
    n_par_value <- length(unique(df_plot$parvalue))
    # colour palette
    if (!grepl("age_exp_par_bins",sel_par)){ colorpal=colorRampPalette(colors=c("orange","red"))(n_par_value)} else  {
      colorpal=colorRampPalette(colors=c("blue","grey","red"))(n_par_value) }
    p <- ggplot(df_plot,aes(x=factor(epi_year),color=factor(parvalue),group=parvalue)) + facet_wrap(~agegroup_broad,scales="free_y") + 
      geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge(width=dodge_val),alpha=0.3,size=28/n_par_value) + #
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
    if (!dir.exists("median_interquant_by_param_value")) {dir.create(paste0(foldername,"median_interquant_by_param_value"))}
    subfldr_name <- "median_interquant_by_param_value/attackrate_sum_seas_share/below5y/"
    if (!dir.exists(paste0(foldername,subfldr_name))) {dir.create(paste0(foldername,subfldr_name))}
    ggsave(paste0(foldername,subfldr_name,"summ_stats_relative_2019_",paste0(sel_var_filename,collapse="_"),"_",sel_par,".png"),
           width=28,height=16,units="cm")
    print(paste0(c(sel_vars[k_plot_var],sel_pars[k_plot_par]),collapse=", "))
  }
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Calculate the average age of infection (by indiv parsets)
parsets_mean_age_inf <- left_join(results_summ_all_hosp,rsv_age_groups %>% select(agegroup_name,mean_age_weighted),by="agegroup_name") %>% 
  group_by(par_id,epi_year) %>% summarise(PC1=unique(PC1),seasforce_peak=unique(seasforce_peak),waning=round(unique(waning)),
                                          seasforc_width_wks=unique(seasforc_width_wks),R0=unique(R0), 
        mean_age_inf_tot_under_5=sum(12*inf_tot[agegroup<=7]*mean_age_weighted[agegroup<=7]/sum(inf_tot[agegroup<=7])),
        mean_age_inf_seas_under_5=sum(12*inf_in_seas[agegroup<=7]*mean_age_weighted[agegroup<=7]/sum(inf_in_seas[agegroup<=7])),
        mean_age_hosp_tot_under_5=sum(12*hosp_tot[agegroup<=7]*mean_age_weighted[agegroup<=7]/sum(hosp_tot[agegroup<=7])),
        mean_age_hosp_seas_under_5=sum(12*hosp_seas[agegroup<=7]*mean_age_weighted[agegroup<=7]/sum(hosp_seas[agegroup<=7]))) %>% 
  pivot_longer(!c(par_id,epi_year,PC1,seasforce_peak,waning,seasforc_width_wks,R0)) %>%
  group_by(par_id,name) %>% summarise(epi_year,PC1,seasforce_peak,waning,seasforc_width_wks,R0,value,
                                      value_norm=ifelse(epi_year==2019,1,value-value[epi_year==2019])) %>% relocate(name,.after=R0) %>%
  mutate(age_exp_par_bins=findInterval(PC1,seq(-1,1,by=1/5))) %>% group_by(age_exp_par_bins) %>% 
  mutate(age_exp_par_bins=round(mean(PC1),1)) %>% relocate(age_exp_par_bins,.after=PC1)
###########################################################
# summary plot (median, interquartile range)
summ_mean_age_infs <- parsets_mean_age_inf %>% group_by(epi_year,name) %>% summarise(mean=mean(value_norm),
                  median=median(value_norm),ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
                  ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% filter(epi_year>2019)
###########################################################
# segment age_exp_dep into x values, calculate summary statistics for each param value
summ_mean_age_inf_byvalue <- parsets_mean_age_inf %>% select(!c(PC1,value)) %>% mutate(waning=round(waning)) %>%
  relocate(age_exp_par_bins,.after=R0) %>% rename(varname=name) %>% pivot_longer(!c(par_id,epi_year,varname,value_norm)) %>% 
  rename(parname=name,parvalue=value) %>% relocate(c(varname,value_norm),.after=parvalue) %>%
  mutate(value_norm=ifelse(varname %in% c("seas_share"), value_norm*100,value_norm)) %>%
  group_by(epi_year,parname,parvalue,varname) %>% summarise(mean=mean(value_norm),median=median(value_norm),
                  ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
                  ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% filter(epi_year>2019)
###########################################################
# PLOT summary stats (mean age inf)
# Figure 4
sel_vars <- c("mean_age_hosp_seas_under_5","mean_age_hosp_tot_under_5")
for (k_plot in 1:length(sel_vars)) {
  dodge_val=0.9; df_plot <- summ_mean_age_infs %>% filter(epi_year>2020 & (name %in% sel_vars[k_plot])) %>% 
    mutate(name=case_when(grepl("mean_age_hosp_seas_under_5",name) ~ "mean age of in-season hospitalisations", 
                          grepl("mean_age_hosp_tot_under_5",name) ~ "mean age of cumulative hospitalisations"))
  ggplot(df_plot,aes(x=factor(epi_year))) + 
    geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge(width=dodge_val),alpha=0.3,size=12,show.legend=FALSE) +
    geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),position=position_dodge(width=dodge_val),alpha=0.6,size=12) + 
    geom_hpline(aes(y=median),position=position_dodge(width=dodge_val),width=1/6,size=1.25,color="black") + # 
    geom_vline(xintercept=(0:4)+1/2,size=1/2) + geom_hline(yintercept=0,linetype="dashed",size=1/2) +
    scale_x_discrete(expand=expansion(0,0)) + scale_y_continuous(breaks=(-10:15)) +
    xlab("") + ylab(paste0(unique(df_plot$name)," (change from 2019 in months)")) + 
    theme_bw() + standard_theme + theme(strip.text=element_text(size=15),legend.text=element_text(size=15),
                  axis.text.x=element_text(size=13),axis.text.y=element_text(size=13),legend.title=element_text(size=15))
  # create folders
  subfldr_name <- "median_interquant_by_param_value/mean_age_all_param/"
  if (!dir.exists(paste0(foldername,subfldr_name))) {dir.create(paste0(foldername,subfldr_name))}
  # save
  ggsave(paste0(foldername,subfldr_name,gsub("-|\\s","_",sel_vars[k_plot]),".png"),width=28,height=16,units="cm")
  print(sel_vars[k_plot])
}

###########################################################
# Plot of shift in average age of inf toggled by param values
# Figure 8, SI Fig 10, 11
sel_vars <- c("mean_age_hosp_seas_under_5","mean_age_hosp_tot_under_5")
sel_pars <- c("age_exp_par_bins","R0","seasforc_width_wks","seasforce_peak","waning")
for (k_plot_var in 1:length(sel_vars)) {
  for (k_plot_par in 1:length(sel_pars)) {
    sel_par <- sel_pars[k_plot_par]; dodge_val=1
    df_plot <- summ_mean_age_inf_byvalue %>% filter(epi_year>2020 & (varname %in% sel_vars[k_plot_var]) & (parname %in% sel_par)) %>%
      mutate(varname=case_when(grepl("mean_age_hosp_seas_under_5",varname) ~ "mean age of in-season hospitalisations",
                               grepl("mean_age_hosp_tot_under_5",varname) ~ "mean age of hospitalisations")) %>% 
      mutate(parname=case_when(grepl("age_exp_par_bins",parname) ~ "exposure (-1) <-> age (1)", 
                               grepl("seasforc_width_wks",parname) ~ "season width (weeks)",
                               grepl("seasforce_peak",parname) ~ "seasonal forcing (above baseline)",
                               grepl("R0",parname) ~ "R0 (baseline)",grepl("waning",parname) ~ "waning (days)"))
    n_par_value <- length(unique(df_plot$parvalue))
    # colour palette
    if (!grepl("age_exp_par_bins",sel_par)){ colorpal=colorRampPalette(colors=c("orange","red"))(n_par_value)} else  {
      colorpal=colorRampPalette(colors=c("blue","grey","red"))(n_par_value) }
    ggplot(df_plot,aes(x=factor(epi_year),color=factor(parvalue),group=parvalue)) + 
      geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),position=position_dodge(width=dodge_val),alpha=0.6,size=52/n_par_value) + #
      geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge(width=dodge_val),alpha=0.3,size=28/n_par_value) + #
      geom_hpline(aes(y=median),position=position_dodge(width=dodge_val),width=(1/n_par_value)*0.75,size=1.25,color="black") + 
      geom_vline(xintercept=(0:4)+1/2,size=1/5) + labs(color=unique(df_plot$parname)) + 
      scale_x_discrete(expand=expansion(0.02,0)) + geom_hline(yintercept=0,linetype="dashed",size=1/2) +
      xlab("") + ylab(paste0(unique(df_plot$varname)," (change from 2019 level)")) + scale_y_continuous(breaks=-10:10) +
      theme_bw() + standard_theme + theme(strip.text=element_text(size=15),axis.text.x=element_text(size=14),
              axis.text.y=element_text(size=14),legend.text=element_text(size=14),legend.title=element_text(size=16),
              legend.position=ifelse(grepl("expos|forcing",unique(df_plot$parname)),"bottom","right")) + scale_color_manual(values=colorpal)
    subfldr_name <- "median_interquant_by_param_value/mean_age_by_paramval/h35_w18/"
    if (!dir.exists(paste0(foldername,subfldr_name))) {dir.create(paste0(foldername,subfldr_name))}
    # save
    ggsave(paste0(foldername,subfldr_name,sel_vars[k_plot_var],"_",sel_par,".png"),width=35,height=18,units="cm")
    print(paste0(c(sel_vars[k_plot_var],sel_pars[k_plot_par]),collapse=", "))
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Calculate peak cases/hospitalisations and length of period with cases > some threshold value

norm_seas_length_wk <- round(length(as.Date("2020-10-07"):as.Date("2021-03-03"))/7,1)
hosp_sum_prob_broad_agegr <- hosp_probabilities %>% mutate(
  agegroup_broad=c("<1y","1-2y","2-5y","5+y")[findInterval(factor(agegroup_name,levels=unique(agegroup_name)),c(2,4,7)+1)+1]) %>%
  group_by(agegroup_broad) %>% summarise(hosp_sum=sum(hosp_num_from_per_inf_prob)) %>% 
  mutate(hosp_per_week_season=hosp_sum/norm_seas_length_wk)
# download `dyn_parsets_main.zip` from https://drive.google.com/file/d/12ohuGEPrVnOxazXnxEGZGwJIwj16frCc/view?usp=sharing
# and unzip into the folder
# LOAD dynamics of simulations (this takes a few minutes)
dyn_all_parsets_broad_age <- left_join( bind_rows(lapply(list.files(foldername,pattern="dyn_parsets.*csv"), 
                   function(x) read_csv(file=paste0(foldername,x)) %>% filter(par_id %in% parsets_regular_dyn$par_id) %>% 
      mutate(date=t-min(t)+as.Date(start_date_dyn_save)) %>% select(!c(t,name)) %>% 
      filter(date>=as.Date("2018-10-01")&date<=as.Date("2024-04-01")) %>% group_by(agegroup,date,par_id) %>% summarise(value=sum(value)) %>%
      group_by(par_id,agegroup) %>% mutate(value=round(roll_sum(value,n=7,fill=NA,align="right",by=7))) %>% filter(!is.na(value))  )  ), 
      hosp_probabilities %>% mutate(agegroup=as.numeric(factor(agegroup_name,levels=unique(agegroup_name)))) %>%
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
# left-join with partable to have input parameters
parsets_max_incid_seas_length <- left_join(left_join(summ_dyn_max_incid_seas_length,pred_pca %>% select(par_id,PC1),by="par_id"),
            partable_regular_dyn %>% select(par_id,omega,seasforc_width_wks,R0,seasforce_peak) %>% mutate(omega=1/omega) %>% 
            rename(waning=omega), by="par_id") %>% relocate(c(name,max_value,sum_value,seas_length_wk),.after=seasforce_peak) %>%
  rename(varname=name)%>% pivot_longer(!c(epi_year,par_id,agegroup_broad,PC1,waning,seasforc_width_wks,R0,seasforce_peak,varname)) %>%
  rename(vartype=name) %>% group_by(par_id,agegroup_broad,PC1,waning,seasforc_width_wks,R0,seasforce_peak,varname,vartype) %>% 
  mutate(value_norm=ifelse(grepl("seas_length_wk",vartype),value-value[epi_year==2018],value/value[epi_year==2018])) %>%
  mutate(age_exp_par_bins=findInterval(PC1,seq(-1,1,by=1/5))) %>% relocate(age_exp_par_bins,.after=PC1) %>%
  group_by(age_exp_par_bins) %>% mutate(age_exp_par_bins=round(mean(PC1),1))
##############################################################
# summary statistics collapse across all param values
summ_max_incid_seas_length <- parsets_max_incid_seas_length %>% 
  group_by(agegroup_broad,epi_year,varname,vartype) %>% summarise(mean=mean(value_norm),
            median=median(value_norm),ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
            ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2])
##############################################################
# calculate summary statistics for each param value
summ_max_incid_seas_length_byvalue <- parsets_max_incid_seas_length %>% mutate(waning=round(waning)) %>% select(!c(PC1,value)) %>% 
  relocate(age_exp_par_bins,.after=R0) %>% pivot_longer(!c(par_id,agegroup_broad,epi_year,varname,vartype,value_norm)) %>% 
  rename(parname=name,parvalue=value) %>% relocate(c(varname,value_norm),.after=parvalue) %>%
  group_by(agegroup_broad,epi_year,parname,parvalue,varname,vartype) %>% summarise(mean=mean(value_norm),median=median(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
            ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% filter(epi_year>2019)
##############################################################
# PLOT peak value of cases/hospitalisations + season length, disaggregated by parameter values
# Figure 7-8, SI 7-10
sel_vars <- c("incid_case","incid_hosp"); sel_pars <- c("age_exp_par_bins","R0","seasforc_width_wks","seasforce_peak","waning")
sel_vartypes <- c("max_value","seas_length_wk")
for (k_plot_var in 1:length(sel_vars)) {
  for (k_plot_par in 1:length(sel_pars)) {
    for (k_plot_vartype in 1:length(sel_vartypes)) {
      sel_vartype<-sel_vartypes[k_plot_vartype]; sel_par <- sel_pars[k_plot_par]; dodge_val=1
      df_plot <- summ_max_incid_seas_length_byvalue %>% 
        filter(epi_year>2020 & (varname %in% sel_vars[k_plot_var]) & (parname %in% sel_par) & (vartype %in% sel_vartype) & 
                 (!agegroup_broad %in% "5+y")) %>%
        mutate(varname=case_when(grepl("incid_case",varname) ~ "cases",grepl("incid_hosp",varname) ~ "hospitalisations"),
               vartype=case_when(grepl("max_value",vartype) ~ "peak", grepl("seas_length_wk",vartype) ~ "above baseline (weeks)")) %>% 
        mutate(parname=case_when(grepl("age_exp_par_bins",parname) ~ "exposure (-1) <-> age (1)", 
                                 grepl("seasforc_width_wks",parname) ~ "season width (weeks)", 
                                 grepl("seasforce_peak",parname) ~ "seasonal forcing (above baseline)",
                                 grepl("R0",parname) ~ "R0 (baseline)",grepl("waning",parname) ~ "waning (days)"))
      ylab_tag <- paste0(paste0(unique(df_plot$vartype)," ",unique(df_plot$varname)),
                         ifelse(grepl("above",unique(df_plot$vartype))," (change from 2019 level)"," (relative to pre-NPI)"))
      n_par_value <- length(unique(df_plot$parvalue))
      # colour palette
      if (!grepl("age_exp_par_bins",sel_par)){ colorpal=colorRampPalette(colors=c("orange","red"))(n_par_value)} else  {
        colorpal=colorRampPalette(colors=c("blue","grey","red"))(n_par_value) }
      p <- ggplot(df_plot,aes(x=factor(epi_year),color=factor(parvalue),group=parvalue)) + facet_wrap(~agegroup_broad,scales="free_y") + 
        # geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge(width=dodge_val),alpha=0.3,size=28/n_par_value) + #
        geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),position=position_dodge(width=dodge_val),alpha=0.6,size=24/n_par_value) + #
        geom_hpline(aes(y=median),position=position_dodge(width=dodge_val),width=(1/n_par_value)*0.75,size=0.8,color="black") + 
        geom_vline(xintercept=(0:4)+1/2,size=1/5) + labs(color=unique(df_plot$parname)) + # geom_hline(yintercept=) +
        scale_x_discrete(expand=expansion(0.02,0)) + xlab("") + ylab(ylab_tag) + theme_bw() + standard_theme + 
        theme(strip.text=element_text(size=15),axis.text.x=element_text(size=13),axis.text.y=element_text(size=12),
              legend.text=element_text(size=11),legend.title=element_text(size=12),
              legend.position=ifelse(grepl("expos|forcing",unique(df_plot$parname)),"bottom","right")) + scale_color_manual(values=colorpal)
      if (grepl("season peak|cases in-season|attack",sel_vars[k_plot_var])) { #  # + scale_y_log10() 
        if (grepl("season peak",sel_vars[k_plot_var])) {break_vals <- (-5:15)*10} else {break_vals <- (-10:10)*10}
        p <- p + geom_hline(yintercept=0,linetype="dashed",size=1/2) + scale_y_continuous(breaks=break_vals) } else {
          p <- p + geom_hline(yintercept=1,linetype="dashed",size=1/2)}; p
      # save
      subfldr_name<-"median_interquant_by_param_value/peak_duration/below5y/h28_w16/"
      if (!dir.exists(paste0(foldername,subfldr_name))) {dir.create(paste0(foldername,subfldr_name))}
      plot_filename <- paste0(foldername,subfldr_name,"summ_stats_",
                              paste0(c(sel_vars[k_plot_var],sel_vartype),collapse="_"),"_",sel_par,".png")
      ggsave(plot_filename,width=28,height=16,units="cm")
      print(paste0(c(sel_vars[k_plot_var],sel_pars[k_plot_par],sel_vartype),collapse=", "))
    }
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
cumul_peak_meanage_hosp_byvalue <- bind_rows(
  # cumulative hospitalisations
  parsets_broad_age_groups %>% filter(name %in% "hosp_tot" & epi_year>2020) %>% ungroup() %>% rename(varname=name) %>% 
    select(!c(value,PC1,par_id)) %>% pivot_longer(c(age_exp_par_bins,seasforce_peak,waning,seasforc_width_wks,R0)) %>% 
    rename(parname=name,parvalue=value) %>% relocate(value_norm,.after=parvalue) %>% group_by(agegroup_broad,epi_year,parname,varname) %>%
    summarise(mean_val=mean(value_norm),med_val=median(value_norm),
              min_val=mean(value_norm[parvalue==ifelse(parname %in% c("waning","age_exp_par_bins"),max(parvalue),min(parvalue))]),
              max_val=mean(value_norm[parvalue==ifelse(parname %in% c("waning","age_exp_par_bins"),min(parvalue),max(parvalue))])),
  # ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
  # ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2]
  # peak hospitalisations
  parsets_max_incid_seas_length %>% filter(varname %in% "incid_hosp" & vartype %in% "max_value" & epi_year>2020) %>% ungroup() %>% 
    select(!c(value,PC1,par_id)) %>% pivot_longer(c(age_exp_par_bins,seasforce_peak,waning,seasforc_width_wks,R0)) %>% 
    rename(parname=name,parvalue=value) %>% relocate(value_norm,.after=parvalue) %>% group_by(agegroup_broad,epi_year,parname,varname) %>%
    summarise(mean_val=mean(value_norm),med_val=median(value_norm),min_val=mean(value_norm[parvalue==min(parvalue)]),
              max_val=mean(value_norm[parvalue==max(parvalue)])),
  # mean age
parsets_mean_age_inf %>% filter(name %in% "mean_age_hosp_tot_under_5" & epi_year>2020) %>% ungroup() %>%
 select(!c(value,PC1,par_id)) %>% rename(varname=name) %>% pivot_longer(c(age_exp_par_bins,seasforce_peak,waning,seasforc_width_wks,R0)) %>%
    rename(parname=name,parvalue=value) %>% relocate(value_norm,.after=parvalue) %>% group_by(epi_year,parname,varname) %>%
    summarise(mean_val=mean(value_norm),med_val=median(value_norm),min_val=mean(value_norm[parvalue==min(parvalue)]),
              max_val=mean(value_norm[parvalue==max(parvalue)])) ) %>% filter(!agegroup_broad %in% "5+y")

# PLOT
# install.packages("ggdist"); library(ggdist)
# change in CUMUL and PEAK HOSP
p_cumul_peak_summ <- ggplot(cumul_peak_meanage_hosp_byvalue %>% mutate(epi_year=as.character(epi_year)) %>%
         filter(epi_year<2023 & (!parname %in% "seasforc_width_wks") & !is.na(agegroup_broad)) %>%
       mutate(parname=case_when(grepl("age_exp_par_bins",parname) ~ "exposure vs age",grepl("seasforce_peak",parname) ~ "seasonal forcing",
        grepl("R0",parname) ~ "baseline R0",grepl("waning",parname) ~ "waning rate"),
    varname=case_when(grepl("hosp_tot",varname) ~ "cumulative",grepl("incid_hosp",varname) ~ "peak"),
    epi_year=case_when(epi_year %in% "2021" ~ "2021-2022",epi_year %in% "2022" ~ "2022-2023")) %>% filter(epi_year %in% "2021-2022"), 
    aes(x=mean_val,y=agegroup_broad)) + 
  geom_interval(aes(group=parname,xmin=min_val,xmax=max_val,color=factor(parname)),position=position_dodge(width=dodge_val),size=10) +
  geom_vpline(aes(group=parname),position=position_dodge(width=dodge_val),color="black",size=1.6,height=0.2) + # 
  facet_grid(~varname) + # agegroup_broad
  geom_hline(yintercept=(0:4)+1/2,size=1/2) + geom_vline(xintercept=1,size=1/2,linetype="dashed") + 
  scale_x_log10(breaks=c(0.3,0.5,0.75,1,1.5,2,3)) + scale_y_discrete(expand=expansion(0,0)) + theme_bw() + standard_theme + 
  theme(strip.text=element_text(size=18),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),
        legend.text=element_text(size=17),legend.position="top",axis.title.x=element_text(size=16)) + 
  xlab("relative hospitalisation risk compared to pre-pandemic years") + ylab("") + labs(color=""); p_cumul_peak_summ
# SAVE
ggsave(paste0(foldername,"median_interquant_by_param_value/summary_range/cumul_peak_hosp_summary_plot_4params.png"),
       width=28,height=16,units="cm")

# shift in average age of hosp
p_aver_age_shift_summ <- ggplot(cumul_peak_meanage_hosp_byvalue %>% mutate(epi_year=as.character(epi_year)) %>%
                              filter(epi_year==2021 & (!parname %in% "seasforc_width_wks") & is.na(agegroup_broad)) %>%
  mutate(parname=case_when(grepl("age_exp_par_bins",parname) ~ "exposure vs age",
         grepl("seasforce_peak",parname) ~ "seasonal forcing",grepl("R0",parname) ~ "baseline R0",grepl("waning",parname) ~ "waning rate"),
  varname=case_when(grepl("mean_age",varname) ~ "shift in age of hospitalisations (<5y)"),
  epi_year=case_when(epi_year %in% "2021" ~ "2021-2022",epi_year %in% "2022" ~ "2022-2023")), 
  aes(x=mean_val,y=1)) + facet_grid(~varname) + # ,y=factor(epi_year)
  geom_interval(aes(group=parname,xmin=min_val,xmax=max_val,color=factor(parname)),position=position_dodge(width=dodge_val),size=9) +
  geom_vpline(aes(group=parname),position=position_dodge(width=dodge_val),color="black",size=1.5,height=0.2) + # 
  geom_vline(xintercept=0,size=1/3,linetype="dashed") + # facet_grid(agegroup_broad~varname,scales="free_x") + 
  geom_hline(yintercept=(0:2)+1/2,size=1/2) + scale_x_continuous(breaks=-3:5) + scale_y_discrete(expand=expansion(0,0)) +
  theme_bw() + standard_theme + theme(axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),axis.title.x=element_text(size=16),
      legend.position="null",strip.text=element_text(size=18)) + xlab("shift in average age (months)") + ylab("") + 
  labs(color="")
p_aver_age_shift_summ
# legend.text=element_text(size=13),legend.title=element_text(size=14),
ggsave(paste0(foldername,"median_interquant_by_param_value/summary_range/aver_age_hosp_summary_plot_4params.png"),
  width=28,height=6,units="cm")

# concat figures
# FIGURE 3 in main text
# plot_grid(p_cumul_peak_summ, p_aver_age_shift_summ,labels="AUTO",rel_widths=c(1.5,1))
if (!any(grepl("cowplot",row.names(installed.packages())))) { install.packages("cowplot") } else {library(cowplot)}
theme_set(theme_cowplot(font_size=18))
plot_grid(p_cumul_peak_summ, p_aver_age_shift_summ,labels="AUTO",ncol=1,rel_heights=c(3,1)) # rel_widths=c(2,1)
# save
ggsave(paste0(foldername,"median_interquant_by_param_value/summary_range/comb_summary_plot_4params_1col.png"),width=40,height=32,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Evaluating the dynamics as a function of parameters
summ_dyn_all_parsets_broad_age <- left_join(dyn_all_parsets_broad_age, 
                                            left_join(partable_regular_dyn %>% select(par_id,seasforc_width_wks,R0,seasforce_peak,omega),
                                                      pred_pca %>% select(par_id,PC1),by="par_id"), by="par_id") %>% 
  mutate(age_exp_par_bins=findInterval(PC1,seq(-1,1,by=1/5))) %>% group_by(age_exp_par_bins) %>% 
  mutate(age_exp_par_bins=round(mean(PC1),1)) %>% rename(incid_case=value) %>% pivot_longer(c(incid_case,incid_hosp)) %>%
  group_by(agegroup_broad,date,age_exp_par_bins,name) %>% summarise(mean=mean(value),median=median(value),
                ci50_low=quantile(value,c(0.25,0.75))[1],ci50_up=quantile(value,c(0.25,0.75))[2],
                ci95_low=quantile(value,c(0.025,0.975))[1],ci95_up=quantile(value,c(0.025,0.975))[2]) %>% rename(varname=name) %>% 
  pivot_longer(c(mean,median,ci50_low,ci50_up,ci95_low,ci95_up)) %>% mutate(epi_year=ifelse(week(date)>=42,year(date),year(date)-1)) %>% 
  rename(metric=name) %>% group_by(agegroup_broad,age_exp_par_bins,varname,metric) %>% # n_year=year(date)
  mutate(peak_2019=max(value[epi_year %in% 2019])) %>% ungroup() %>% mutate(value_norm=value/peak_2019) %>% select(!peak_2019)

##############################################################
# Plot dynamics faceted by the age vs immunity dependence of susceptibility: Figure 9
sel_vars <- c("incid_case","incid_hosp"); n_par<-length(unique(summ_dyn_all_parsets_broad_age$age_exp_par_bins))
colorpal <- colorRampPalette(colors=c("blue","grey","red"))(n_par); lim_dates<-as.Date(c("2022-04-01","2023-04-01"))
for (k_age in 1:4) {
  for (k_date in 1:2) {
    for (value_type in c("value","value_norm")) {
      val_type_excl <- c("value","value_norm")[!(c("value","value_norm") %in% value_type)]
      sel_agegr<-unique(summ_dyn_all_parsets_broad_age$agegroup_broad)[k_age]
      df_plot <- summ_dyn_all_parsets_broad_age %>% filter(varname %in% sel_vars[2] & date>=as.Date("2021-05-01") &
                  date<=lim_dates[k_date] & agegroup_broad %in% sel_agegr & (metric %in% c("median","ci50_low","ci50_up"))) %>% 
        select(!c(varname,!!val_type_excl)) %>% pivot_wider(names_from = metric,values_from=!!value_type)
      p<-ggplot(df_plot,aes(x=date,group=age_exp_par_bins,color=factor(age_exp_par_bins))) + 
        geom_line(aes(y=median),size=1.1) + geom_ribbon(aes(ymin=ci50_low,ymax=ci50_up,fill=factor(age_exp_par_bins)),color=NA,alpha=0.2) + 
        scale_color_manual(values=colorpal) + scale_fill_manual(values=colorpal) + # facet_wrap(~agegroup_broad,scales="free_y",nrow=3) + 
        scale_x_date(date_breaks="1 month",expand=expansion(0.01,0)) + theme_bw() + standard_theme + 
        theme(legend.position="top",axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
              axis.title.y=element_text(size=16),legend.text=element_text(size=15),legend.title=element_text(size=16)) + 
        labs(color="exposure (-1) <-> age (1)",fill="exposure (-1) <-> age (1)")+xlab("")+
        ylab(paste0("weekly hospitalisations (",sel_agegr,ifelse(grepl("norm",value_type),", normalised by 2019 peak",""),")"))
      if (grepl("norm",value_type)) {p<-p+geom_hline(yintercept=1,linetype="dashed",size=1/3)}; p
      # save
      if (!dir.exists(paste0(foldername,"dynamics/"))) {dir.create(paste0(foldername,"dynamics/"))}
      abs_val_folder <- ifelse(!grepl("norm",value_type),"abs_val/","")
      if (!dir.exists(paste0(foldername,"dynamics/",abs_val_folder))) {dir.create(paste0(foldername,"dynamics/",abs_val_folder))}
      plot_fn<-paste0(foldername,"dynamics/",abs_val_folder,"weekly_hosp_by_ageexpdep_",ifelse(grepl("norm",value_type),"norm_2019_peak_",""),
                      gsub("-","_",gsub("<","below",sel_agegr)),"_",substr(lim_dates,1,4)[k_date],".png")
      ggsave(plot_fn,width=25,height=20,units="cm"); print(gsub(foldername,"",plot_fn))
    }
  }
}
##############################################################
# Faceted by age/exp dependency: SI Fig 12-13
value_type<-"value_norm"
for (k_date in 1:2) {
  sel_agegr<-unique(summ_dyn_all_parsets_broad_age$agegroup_broad)[k_age]
  df_plot <- summ_dyn_all_parsets_broad_age %>% filter(varname %in% sel_vars[2] & date>=as.Date("2021-05-01") &
    !(agegroup_broad %in% "5+y")&date<=lim_dates[k_date]&(metric %in% c("median","ci50_low","ci50_up"))) %>% select(!c(varname,value)) %>%
     pivot_wider(names_from=metric,values_from=!!value_type) %>% filter(as.numeric(factor(age_exp_par_bins)) %in% c(1,2,4,7,9,10)) %>%
    rename(`exposure (-1) <-> age (1)`=age_exp_par_bins)
  p<-ggplot(df_plot,aes(x=date,group=agegroup_broad,color=factor(agegroup_broad))) + 
    geom_line(aes(y=median),size=1.1) + geom_ribbon(aes(ymin=ci50_low,ymax=ci50_up,fill=factor(agegroup_broad)),color=NA,alpha=0.2) + 
    facet_wrap(~`exposure (-1) <-> age (1)`,nrow=2,labeller=labeller(`exposure (-1) <-> age (1)`=label_both)) + 
    geom_vline(xintercept=as.Date(c("2021-10-15","2022-03-04","2022-10-15","2023-03-04")),size=1/3,linetype="dashed") +
    scale_x_date(date_breaks="2 month",expand=expansion(0.01,0)) + theme_bw() + standard_theme + 
    theme(legend.position="top",strip.text=element_text(size=15),axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
          axis.title.y=element_text(size=16),legend.text=element_text(size=16)) + 
    labs(color="",fill="")+ xlab("") +ylab("weekly hospitalisations (normalised by 2019 peak)")
  if (grepl("norm",value_type)) {p<-p+geom_hline(yintercept=1,linetype="dashed",size=1/3)}; p
  # save
  if (!dir.exists(paste0(foldername,"dynamics/facet_ageexpdep/"))) {dir.create(paste0(foldername,"dynamics/facet_ageexpdep/"))}
  plot_fn<-paste0(foldername,"dynamics/facet_ageexpdep/","weekly_hosp_by_ageexpdep_",ifelse(grepl("norm",value_type),"norm_2019_peak_",""),
                  "until",substr(lim_dates,1,4)[k_date],".png")
  ggsave(plot_fn,width=30,height=20,units="cm"); print(gsub(foldername,"",plot_fn))
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Plotting individual trajectories from the parameter sampling

ggplot(dyn_all_parsets_broad_age %>% filter(par_id %in% unique(par_id)[1:11]),
       aes(x=date,y=incid_hosp,group=par_id,color=factor(par_id))) + geom_line() + facet_wrap(~agegroup_broad,scales="free_y") + 
  scale_x_date(date_breaks="3 month",expand=expansion(0.01,0)) + theme_bw() + standard_theme