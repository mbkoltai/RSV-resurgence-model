# This script is for reproducing the results and figures of the manuscript at [...]
# Mihaly Koltai, May/2022
####

# clear workspace
rm(list=ls())
# To set the path we need the "here" package
if (!any(row.names(installed.packages()) %in% "here")) {install.packages("here")}; library(here)
# load constant parameters and functions for simulations, specify folder where inputs are stored
source(here("load_params.R"))
foldername <- "repo_data/"
# load estimated attack rates
attack_rate_tol_factor = 2
estim_attack_rates <- read_csv(here("repo_data/kenya_attack_rates.csv")) %>% 
  mutate(min_est=1e2*sympt_attack_rate/attack_rate_tol_factor,
         max_est=1e2*attack_rate_tol_factor*sympt_attack_rate)
# % cases within season (filtering parameter sets)
npi_dates<-as.Date(c("2020-03-26","2021-05-17"))
# set up the table of parameter vectors by Latin Hypercube Sampling (LHS)
partable <- data.frame(randomLHS(n=32e3,k=7))
colnames(partable) <- c("exp_dep","age_dep","seasforc_width_wks","R0","seasforce_peak","omega","peak_week")
# convert to relevant ranges, columns: 
# 1) expdep 2) agedep 3) peak_width 4) R0 5) peak_height 6) waning 7) peak week
partable[,1] <- qunif(partable[,1],min=3/10,max=1.25) # 1) expdep
partable[,2] <- qunif(partable[,2],min=1/15,max=1/3) # 2) agedep
partable[,3] <- qgamma(partable[,3],shape=10,rate=2) # 3) peak_width
partable[,4] <- qgamma(partable[,4],shape=14,rate=8) # 4) R0
partable[,5] <- qunif(partable[,5],min=0.2,max=1.5) # 5) peak height
partable[,6] <- qnorm(partable[,6],mean=350,sd=50); partable[,6] <- 1/partable[,6] # 6) waning
partable[,7] <- qunif(partable[,7],min=43,max=48) # 7) peak week
if (any(partable<0)) { partable=abs(partable); message("negative values") }

# initial sampling showed that accepted parsets satisfy the condition: exp_dep < 1.65 - 4.5*age_dep
partable <- partable %>% filter(exp_dep+4.5*age_dep<1.85)

# grid search:
# list(exp_dep=seq(1/4,2,1/8),
#                age_dep=seq(1/8,1,1/16),
#                seasforc_width_wks=c(3,5,7),
#                R0=1+(0:5)/10,
#                seasforce_peak=c(3/4,1,5/4),
#                omega=c(1/250,1/350,1/450)) %>% expand.grid %>% bind_rows

# calculating the susceptibility parameter (delta_i_j)
l_delta_susc <- 1:nrow(partable) %>%
  lapply( function(n_p) { 
    sapply(1:n_age, function(x) { (1 * exp(-partable$exp_dep[n_p] * (1:3))) / (exp(partable$age_dep[n_p]*x)) }) } )
# create partable with scaling parameters (this parameter scales 'expdep' and 'agedep' to have the desired R0 value)
partable <- partable %>% 
  mutate(par_id=row_number(), 
         const_delta=R0/unlist(lapply(l_delta_susc, function(x) R0_calc_SIRS(C_m,x,rho,n_inf)))) %>% 
         relocate(par_id,.before=exp_dep)
# clear list
rm(l_delta_susc)
#
# filtering param sets: selected parsets are along the line `age=-exp/3+5/6` (and the point (age,exp)=(1/8,1.75))
narrow_par_search <- FALSE
if (narrow_par_search) {
partable_full_linear_kage_kexp <- partable %>% mutate(age_dep_fit=5/6-exp_dep/3) %>% 
                                    filter(abs(age_dep-age_dep_fit)/age_dep<1/3) %>% 
                                    select(!age_dep_fit)
}
# SAVE filtered parameter table
# write_csv(partable,"repo_data/partable_full_lhs.csv")
# check the size of objects (>x Mb) in the workspace by: fcn_objs_mem_use(min_size=1)
# start date of simulations
start_date_dyn_save <- "2016-09-01"
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# load hospitalisation data: `hosp_probabilities` contains hospitalisation/infection probabilities to 
# convert cases to hospitalisations
source(here("fcns", "calc_hosp_rates.R"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Susceptibility as a function of age and exposure
age_exp_dep_uniqvals <- list(exp_dep=seq(1/4,2,1/8),age_dep=seq(1/8,1,1/16),age=1:11,exp=1:3) %>% 
                          expand.grid %>% bind_rows %>% 
                          mutate(suscept_unscaled=exp(-(exp_dep*exp+age_dep*age)))
age_exp_dep_uniqvals <- age_exp_dep_uniqvals %>% 
                mutate(const_delta=1/unlist(lapply(lapply(1:nrow(age_exp_dep_uniqvals), 
                 function(n_p) {sapply(1:n_age,function(x) { 
                   (1*exp(-age_exp_dep_uniqvals$exp_dep[n_p]*(1:3)))/(exp(age_exp_dep_uniqvals$age_dep[n_p]*x))})}),
                 function(x) R0_calc_SIRS(C_m,x,rho,n_inf))),
                 susc_scaled=suscept_unscaled*const_delta)
# PLOT: 
ggplot(age_exp_dep_uniqvals %>% 
          filter(exp_dep %in% seq(1/4,2,3/4) & age_dep %in% seq(1/8,1,1/4)) %>%
          rename(`exposure-dependence`=exp_dep,`age-dependence`=age_dep) %>% 
          mutate(age=factor(rsv_age_groups$agegroup_name[age],levels=unique(rsv_age_groups$agegroup_name))) )  + 
      geom_line(aes(x=age,color=factor(exp),group=exp,y=susc_scaled),size=1.06) + 
      facet_grid(`exposure-dependence`~`age-dependence`,
             labeller=labeller(`exposure-dependence`=label_both,`age-dependence`=label_both)) + 
      scale_y_log10() + theme_bw() + standard_theme + labs(color="exposure") + xlab("age group") + 
      ylab(expression(delta[exp]^(age)))
# ggsave
# ggsave(here::here(foldername,"age_exp_dep.png"),width=22,height=18,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# RUN SIMULATIONS with parallelisation (requires multiple cores)
# The below code does not have to be re-run as results are available: 
# go down to RESULTS and read in the results from the parameter sampling
#
# # `start_date_dyn_save` will define the timepoint from which to save simulation outputs
# # `simul_length_yr` is the length of simulations (in years), `n_post_npi_yr` is the number of years after NPIs
# # `n_core`: number of cores used, `memory_max`: allocated memory (GB)
# # start_date_dyn_save: date to save results from
# simul_length_yr <- 30; n_post_npi_yr<-4; n_core<-64; memory_max <- 8; start_date_dyn_save <- "2016-09-01"
# contact_red<-0
# agegroup_res<-"broad_age" # this if summary stats should be saved in 0-1, 1-2, 2-5, 65+ groups only
# # these are parameters selected by criteria of 1) attack rates 2) seasonal concentration of cases
# partable_filename <- "repo_data/partable_full_lhs.csv"; n_row <- nrow(read_csv(partable_filename))
# # we split the parameter table into `n_core` batches and run them in parallel, the sh file will launch the jobs
# # write the files launching jobs
# command_print_runs <- paste0(c("Rscript fcns/write_run_file.R",n_core,n_row,
#                                 simul_length_yr,n_post_npi_yr,contact_red,
#                                 partable_filename,"SAVE sep_qsub_files",start_date_dyn_save,
#                                 agegroup_res,memory_max),collapse=" ")
# # run this command by:
# system(command_print_runs)
# # run calculation (on a cluster, requires multiple cores) by:
# sh master_start.sh
# # collect & merge results (for summary statistics) by:
# nohup Rscript fcns/collect_save_any_output.R simul_output/parscan/ summ_parsets* results_summ_all.csv keep &
####################################
# Simulations are launched by the file "fcns/parscan_runner_cmd_line.R" receiving command line arguments 
# from the parameter table and calling the function `sirs_seasonal_forc_mat_immun` that generates & solves ODEs. 
# The fcn for `sirs_seasonal_forc_mat_immun` is in "fcns/essential_fcns.R".
# The structure of the force of infection is built outside this function and passed as a set of inputs:
# the input `contmatr_rowvector` is the contact matrix stacked 3x times (bc of 3 levels of infection) 
# and normalised by population sizes to generate the force of infection terms normalised by age group sizes. 
# See SI Methods
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# RESULTS
# results from the parameter sampling after filtering already stored in repo_data/
# starting date of simulations
# start_date_dyn_save <- "2016-09-01"
results_summ_all <- read_csv("simul_output/2e4_parsets/results_summ_all.csv")
# read_csv(here::here(foldername,"results_summ_all.csv")) 
# param sets filtered out bc of attack rates or seasonal concentr
AR_crit=9; seas_conc_crit=9
# score: how many age groups satisfy
fullscan_score_AR_seasconc <- left_join(results_summ_all %>% 
                    filter(epi_year %in% c(2017,2018)) %>% group_by(par_id,agegroup) %>%
                    summarise(attack_rate_perc=mean(attack_rate_perc),seas_share=mean(seas_share)),
                    estim_attack_rates %>% select(!c(attack_rate,sympt_attack_rate,n_test,RSV_posit,RSV_sympt_posit))  %>%
                      mutate(agegroup=as.numeric(factor(agegroup_name,levels=agegroup_name))),by="agegroup") %>%
                    mutate(attack_rate_check=(attack_rate_perc>=min_est & attack_rate_perc<=max_est),
                           seas_share_check=(seas_share>=0.85)) %>% group_by(par_id) %>%
                  summarise(n_attack_rate_check=sum(attack_rate_check), n_seas_share_check=sum(seas_share_check)) 
# summary stats of filtering
filtering_summ_stats <- fullscan_score_AR_seasconc %>% 
                    mutate(attack_rate_fail=n_attack_rate_check<AR_crit,
                        seas_conc_fail=n_seas_share_check<seas_conc_crit,both_fail=attack_rate_fail&seas_conc_fail) %>% 
                    summarise(n_both_fail=sum(both_fail),n_attack_rate_fail=sum(attack_rate_fail)-n_both_fail,
                      n_seas_conc_fail=sum(seas_conc_fail)-n_both_fail,n_accepted=sum(!attack_rate_fail&!seas_conc_fail))
# filter for parameter sets where x/11 age groups satisfy criteria for attack rates and seasonal concentration
check_crit=9/11; sel_yrs<-2019; n_sel_yr=length(sel_yrs)
all_sum_inf_epiyear_age_filtered <- results_summ_all %>% 
            filter(par_id %in% (fullscan_score_AR_seasconc %>% 
                            filter(n_attack_rate_check>=AR_crit & n_seas_share_check>=seas_conc_crit))$par_id)
# this leads to xxx parameter sets:
length(unique(all_sum_inf_epiyear_age_filtered$par_id))
# select parameter sets matching the first two criteria
partable_filtered_AR_seasconc <- partable %>% 
                                   filter(par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id))
# write_csv(partable_filtered_AR_seasconc,here(foldername,"partable_filtered_AR_seasconc.csv"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SI Figure 14: plot total infections in epi-years 2017,17,18 compared to 2019

# Plot cumul incid (relative to 2019) of accepted parsets
results_summ_all %>% 
  filter(epi_year<2020 & par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id)) %>%
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],levels=rsv_age_groups$agegroup_name)) %>%
  group_by(par_id,agegroup_name) %>% mutate(inf_in_seas=inf_in_seas/(inf_in_seas[epi_year==2019]*1+(1/26))) %>%
  filter(epi_year<2019) %>%
ggplot() + geom_jitter(aes(x=factor(epi_year),y=inf_in_seas,group=par_id),alpha=1/8,size=1/2) +
  # geom_segment(data=epiyear_means,aes(x=epi_year-0.4,xend=epi_year+0.4,y=inf_in_seas,yend=inf_in_seas),size=1,color="red") +
  facet_wrap(~agegroup_name,scales = "free") + scale_y_log10() +
  xlab("epi-year") + ylab("number of infections compared to 2019 (1=2019 level)") + 
  theme(legend.position="top") + theme_bw() + standard_theme
# save
# ggsave(here(foldername,"2020_RSV_check_log.png"),width=32,height=20,units="cm")

# Biennial patterns? Number of param sets where cumul incidence is within 15% of 2019 level (or not)
irreg_toler=0.15; scaling_2019=1+1/26
results_summ_all %>% 
  filter(epi_year<2020 & par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id)) %>%
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],levels=rsv_age_groups$agegroup_name)) %>%
  group_by(par_id,agegroup_name) %>% mutate(inf_in_seas=inf_in_seas/(scaling_2019*inf_in_seas[epi_year==2019])) %>%
  filter(epi_year<2019) %>% group_by(par_id,epi_year) %>% 
  summarise(regular_by_agegr=sum(inf_in_seas>1-irreg_toler & inf_in_seas<1+irreg_toler)) %>%
  group_by(epi_year) %>% summarise(regular=sum(regular_by_agegr>9),irregular=n()-regular)

# parsets with regular (annual) vs irregular patterns
reg_irreg_parsets <- results_summ_all %>% 
  filter(epi_year<2020 & par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id)) %>%
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],levels=rsv_age_groups$agegroup_name)) %>%
  group_by(par_id,agegroup_name) %>% mutate(inf_in_seas=inf_in_seas/(scaling_2019*inf_in_seas[epi_year==2019])) %>%
  filter(epi_year<2019) %>% group_by(par_id,epi_year) %>% 
  summarise(regular_by_agegr=sum(inf_in_seas>1-irreg_toler & inf_in_seas<1+irreg_toler)) %>% group_by(par_id) %>%
  summarise(regular=all(regular_by_agegr==11))
# number of regular/irregular 
reg_irreg_parsets %>% summarise(regular=sum(regular),irreg=n()-regular)

# CDF of percentage deviation from 2019
results_summ_all %>% 
  filter(epi_year<2020 & par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id)) %>%
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],levels=rsv_age_groups$agegroup_name)) %>%
  group_by(par_id,agegroup_name) %>% mutate(inf_in_seas=inf_in_seas/(scaling_2019*inf_in_seas[epi_year==2019])) %>%
  filter(epi_year<2019) %>% select(par_id,agegroup,inf_in_seas) %>% mutate(inf_in_seas=abs(1-inf_in_seas)) %>%
ggplot(aes(inf_in_seas)) + stat_ecdf(geom="step") + # ,group=agegroup,color=agegroup
  facet_wrap(~agegroup,nrow=2,labeller=labeller(agegroup=label_both)) + 
  geom_vline(xintercept=0.15,linetype="dashed",size=1/2,color="red") +
  scale_x_log10(limits=c(0.001,1)) + # ,expand=expansion(0.02,0)
  xlab("relative difference in incidence") + ylab("CDF") + labs(color="agegroups")+
  theme_bw() + standard_theme + theme(legend.position="top",strip.text=element_text(size=14),
                                      legend.title=element_text(size=15),legend.text=element_text(size=14),
                                      axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
                                      axis.title.x=element_text(size=16))
# save
ggsave(here("simul_output/interyear_difference_cumul_incid_reg_dyn.png"),width=35,height=25,units="cm")

# they haven't settled or stable biennial? read in one of the dynamic batches to check
dyn_parsets_main1_317 <- read_csv("simul_output/2e4_parsets/dyn_parsets_main1_317.csv") %>%
  filter(par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id)) %>% 
  mutate(date=as.Date("2016-09-01")+t-min(t)) %>% filter(date<as.Date("2020-04-01"))
# plot
ggplot(dyn_parsets_main1_317 %>% 
         filter(date<as.Date("2020-04-01") & par_id %in% reg_irreg_parsets$par_id & agegroup==1)) +
  geom_line(aes(x=date,y=value),alpha=1/2) + facet_wrap(~par_id,scales="free") + xlab("") + ylab("infections") +
  theme_bw() + standard_theme

# % difference in dynamic curves, calculated as: sum(abs(diff(2016,2018)))/sum(mean(2016,2018)) and same for 2017/2019
left_join(dyn_parsets_main1_317 %>% filter(date<as.Date("2020-04-01")),reg_irreg_parsets) %>%
  mutate(week=isoweek(date)) %>% filter(week>=35|week<=9) %>% 
  mutate(epi_year=ifelse((week>=35 & yday(date)>245)|yday(date)>245,year(date),year(date)-1)) %>% 
  filter(epi_year>=2016) %>%
  group_by(agegroup,epi_year,par_id,regular) %>% arrange(epi_year,date) %>% 
  mutate(epi_day=row_number()) %>% filter(epi_day<=180) %>% group_by(par_id,epi_day,agegroup,regular) %>%
  summarise(`diff_2016_18`=abs(value[epi_year==2016]-value[epi_year==2018]),
      `diff_2017_19`=abs(value[epi_year==2017]-value[epi_year==2019]),
      mean_2016_18=mean(value[epi_year %in% c(2016,2018)]),
      mean_2017_19=mean(value[epi_year %in% c(2017,2019)])) %>%
  group_by(agegroup,par_id,regular) %>% summarise(rel_diff_16_18=sum(diff_2016_18)/sum(mean_2016_18),
                                          rel_diff_17_19=sum(diff_2017_19)/sum(mean_2017_19)) %>%
  group_by(par_id,regular) %>% 
  summarise(rel_diff_16_18=mean(rel_diff_16_18)*100,rel_diff_17_19=mean(rel_diff_17_19)*100) %>%
  group_by(regular) %>% summarise(perc_diff_16_18=mean(rel_diff_16_18),perc_diff_17_19=mean(rel_diff_17_19))
# all under 4%, on average under 2% difference
# Plot curves comparing 2016 to 2018 and 2017 to 2019
left_join(dyn_parsets_main1_317 %>% filter(date<as.Date("2020-04-01")),reg_irreg_parsets) %>%
  mutate(week=isoweek(date)) %>% filter(week>=35|week<=9) %>%
  mutate(epi_year=ifelse((week>=35 & yday(date)>245)|yday(date)>245,year(date),year(date)-1)) %>% filter(epi_year>=2016) %>%
  group_by(agegroup,epi_year,par_id,regular) %>% arrange(epi_year,date) %>% 
  mutate(epi_day=row_number()) %>% filter(epi_day<=180 & agegroup==2 & !regular) %>%
  mutate(year_type=epi_year %% 2) %>% group_by(year_type) %>% 
  mutate(year_rank=ifelse(epi_year==min(epi_year),0,1)) %>%
ggplot(aes(x=epi_day,y=value,color=factor(epi_year),group=epi_year,linetype=factor(year_rank))) + 
  labs(linetype="",color="epi-year") +
  geom_line(size=1.1) + facet_wrap(~par_id,scales="free_y") + xlab("day of season") + theme_bw() + standard_theme
# save
# ggsave(here("simul_output/interyear_difference_irregular_parsets.png"),width=30,height=25,units="cm")
ggsave(here("simul_output/interyear_difference_regular_parsets.png"),width=30,height=25,units="cm")

# Plot again the comparison to 2019 for regular parameter sets only
results_summ_all %>% 
  filter(epi_year<2020 & par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id)) %>%
  filter(par_id %in% reg_irreg_parsets$par_id[reg_irreg_parsets$regular]) %>%
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],levels=rsv_age_groups$agegroup_name)) %>%
  group_by(par_id,agegroup_name) %>% mutate(inf_in_seas=inf_in_seas/(inf_in_seas[epi_year==2019]*scaling_2019)) %>%
  filter(epi_year<2019) %>%
  ggplot() + geom_jitter(aes(x=factor(epi_year),y=inf_in_seas,group=par_id),alpha=1/8) +
  facet_wrap(~agegroup_name) + xlab("epi-year") + ylab("number of infections (1=2019 level)") + 
  theme(legend.position="top") + theme_bw() + standard_theme
# save
ggsave(here("simul_output/cumul_incid_compared2019_by_agegrs.png"),width=30,height=25,units="cm")

# keep only parameters with regular annual patterns
partable_regular_dyn <- partable_filtered_AR_seasconc %>% 
  filter(par_id %in% reg_irreg_parsets$par_id[reg_irreg_parsets$regular])
# keep outputs with correct attack rate
results_summ_all_reg_dyn <- results_summ_all %>% filter(par_id %in% partable_regular_dyn$par_id)

# seasonal concentration density plot
results_summ_all_reg_dyn %>% filter(epi_year<2019) %>%
ggplot() + # stat_ecdf(aes(seas_share,group=epi_year,color=factor(epi_year)),geom="step") +
  geom_density(aes(x=seas_share,group=epi_year,color=factor(epi_year))) + labs(color="epi-year") +
  xlab("seasonal (w40≤x≤w13) share of cumulative incidence") + theme_bw() + standard_theme
# save
ggsave(here("simul_output/seas_conc_density_plot.png"),width=30,height=25,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# After filtering by attack rate and seasonal concentration we need full dynamics for *selected* parsets
# --> likelihood calculations

# accepted parsets
all_dynamics_accepted <- bind_rows(lapply(
  list.files("simul_output/2e4_parsets/",pattern="dyn_parsets",full.names=T), 
       function(x) read_csv(x) %>% filter(par_id %in% unique(results_summ_all_reg_dyn$par_id)) ))
write_csv(all_dynamics_accepted,"simul_output/2e4_parsets/all_dynamics_accepted.csv")
# remove timepoints after SARI-Watch data 
all_dynamics_accepted <- all_dynamics_accepted %>% mutate(date=as.Date(start_date_dyn_save)+t-min(t)) %>% 
  filter(date<=as.Date("2020-05-11"))

# we need to also take rejected parsets to make comparison (likelihood vs filtering)
all_dynamics_rejected <- bind_rows(lapply(
  list.files("simul_output/2e4_parsets/",pattern="dyn_parsets",full.names=T),
  function(x) read_csv(x) %>% filter(!par_id %in% unique(results_summ_all_reg_dyn$par_id)) %>%
    filter(par_id %in% sample(par_id,sample(c(16,17),1))) ))
write_csv(all_dynamics_rejected,"simul_output/2e4_parsets/all_dynamics_rejected_1e3.csv")
all_dynamics_rejected <- all_dynamics_rejected %>% mutate(date=as.Date(start_date_dyn_save)+t-min(t)) %>% 
  filter(date<=as.Date("2020-05-11"))

# load RSV hospitalisation counts
SARIwatch_RSVhosp_under5_2018_2020_weekly_counts <- 
  read_csv("data/SARIwatch_RSVhosp_under5_2018_2020_weekly_counts.csv",col_types="ffddd") %>% 
  mutate(wk_n=gsub("\\.","-W",wk_n),wk_n=factor(wk_n,levels=unique(wk_n)),
         date=ISOweek2date(paste0(gsub("\\.","-W",wk_n),"-1")))
# compare SARI-Watch with literature estimates 
# annual RATE of hospitalisations per 100K population: 500/1e5 (18-19), 494/1e5
SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% group_by(year) %>% 
  summarise(annual_cumul_rate=sum(rate_under5yrs))

# under-reporting for under-5 hospitalisations
# annual RATE of hospitalisations from Reeves 2017: <1y: 35.1/1000, 1-4y: 5.31/1000 = 1083.8/1e5
# annual RATE from Taylor 2016: <0.5y: 4184/1e5, 6-23mts: 1272/1e5, 2-4y: 114/1e5 = 822.7/1e5
reported_hosp_rate_per_100k <- c( reeves_2017=(35.1*sum(rsv_age_groups$value[1:2]) + 
                                    5.31*sum(rsv_age_groups$value[3:7]))*100/sum(rsv_age_groups$value[1:7]),
                                taylor_2016=(4184*sum(rsv_age_groups$value[1]) + 1272*sum(rsv_age_groups$value[2:4]) + 
                                    114*sum(rsv_age_groups$value[5:7]))/sum(rsv_age_groups$value[1:7]) )
under_report_factor_under5 = mean( (SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% group_by(year) %>% 
          summarise(annual_cumul_rate=sum(rate_under5yrs)))$annual_cumul_rate)/mean(reported_hosp_rate_per_100k)

# hospitalisation counts for 65+
# for 65+y
SARIwatch_RSVhosp_over65_2018_2020_weekly_counts <- 
  read_csv("data/SARIwatch_RSVhosp_over65y_2018_2020_weekly_counts.csv") %>%
  mutate(wk_n=gsub("-","-W",wk_n),wk_n=factor(wk_n,levels=unique(wk_n)),date=ISOweek2date(paste0(wk_n,"-1")))

# plot
# SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% pivot_longer(!c(year,wk_n,pop_AGE65PLUS,date)) %>% 
#   ggplot(aes(x=wk_n,y=value,group=1)) + geom_line() + geom_point() + facet_grid(name~year,scales="free") + 
#   xlab("") + ylab("count or rate per 100K") + theme_bw() + standard_theme


# rate per 100K from SARI-Watch
SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% group_by(year) %>% summarise(annual_rate=sum(rate_65yplus))
# Estimates from literature (Fleming 2015, https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-015-1218-z/tables/3)
# 65-74: 86/100e3 (62-101); 75+: 234/100e3 (180-291)
over65_hosp_rate_100k_lit_estim = (
  sum(ons_2020_midyear_estimates_uk$value[ons_2020_midyear_estimates_uk$age %in% 65:74])*86 + 
    sum(ons_2020_midyear_estimates_uk$value[76:nrow(ons_2020_midyear_estimates_uk)])*234)/
  sum(ons_2020_midyear_estimates_uk$value[66:nrow(ons_2020_midyear_estimates_uk)])
# under-reporting rate
under_report_factor_over65y <- mean((SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% 
                  group_by(year) %>% summarise(annual_rate=sum(rate_65yplus)))$annual_rate) / 
                  over65_hosp_rate_100k_lit_estim

# calculate weekly hospitalisations of accepted
agegroup_mapping=findInterval(1:11,c(3,5,8,11))+1; names(agegroup_mapping)=paste0("agegroup",1:11)
# calculate hospitalisations
simul_hosp_rate_weekly <- left_join(
  bind_rows(all_dynamics_accepted %>% mutate(accepted=TRUE),all_dynamics_rejected %>% mutate(accepted=FALSE)) %>% 
    filter( (agegroup<=3 | agegroup==5) & date>as.Date("2018-09-15")) %>%
    group_by(t,date,agegroup,par_id) %>% summarise(value=sum(value),accepted=unique(accepted)),
  # data on hospitalisation probs
  hosp_probabilities %>% mutate(broad_age=findInterval(agegroup,c(2,4,7,10)+1)+1) %>% group_by(broad_age) %>%
    # this is not entirely correct, hosp has to be calculated for smaller age groups and then summed
    summarise(prob_hosp_per_infection=sum(prob_hosp_per_infection*size/sum(size)),agegroup=broad_age) %>%
    select(agegroup,prob_hosp_per_infection),by=c("agegroup")) %>% # end of left_join
  mutate(hosp=value*prob_hosp_per_infection) %>% select(!c(prob_hosp_per_infection)) %>%
  mutate(broad_age=ifelse(agegroup<=3,"<5y",">65y")) %>%
  group_by(date,broad_age,par_id) %>% summarise(value=sum(value),hosp=sum(hosp),accepted=unique(accepted)) %>% 
  mutate(year_week=paste0(isoyear(date),"-", isoweek(date)),
                      year_week=factor(year_week,levels=unique(year_week))) %>%
  group_by(year_week,par_id,broad_age) %>% 
  summarise(simul_cases=sum(value),simul_hosp_sum_full=sum(hosp),date=min(date),accepted=unique(accepted) ) %>%
  mutate(simul_hosp_rate_100k=ifelse(grepl("65y",broad_age),
                                     1e5*simul_hosp_sum_full/rsv_age_groups$stationary_popul[11],
                                     1e5*simul_hosp_sum_full/sum(rsv_age_groups$stationary_popul[1:7]) ) )
# concatenate with data, calculate poisson likelihood
simul_hosp_rate_weekly <- left_join(
  left_join(SARIwatch_RSVhosp_under5_2018_2020_weekly_counts,
            SARIwatch_RSVhosp_over65_2018_2020_weekly_counts,by=c("wk_n","date","year")), 
  simul_hosp_rate_weekly, by="date") %>%
  mutate(simul_hosp_scaled=ifelse(grepl("65y",broad_age),
                                  simul_hosp_sum_full*under_report_factor_over65y*(
                                    pop_AGE65PLUS/rsv_age_groups$stationary_popul[11]),
                                  simul_hosp_sum_full*under_report_factor_under5*(
                                    pop_AGEUNDER5/sum(rsv_age_groups$stationary_popul[1:7])) ),
         log_lklh_poiss=dpois(x=ifelse(grepl("65y",broad_age),cases65plustotal,casesunder5total),
                              lambda=simul_hosp_scaled,log=T)) %>% relocate(date,.after=wk_n) %>%
  group_by(par_id,broad_age,accepted) %>% mutate(sum_neg_llh=-sum(log_lklh_poiss,na.rm=T))

# `simul_hosp_scaled` means that simulated hospitalisations were scaled by the under-reporting factor

# calculate likelihoods by parameter set (accepted/rejected)
all_likelihoods_hosp <- simul_hosp_rate_weekly %>% group_by(par_id,accepted,broad_age) %>% 
  summarise(sum_neg_llh=unique(sum_neg_llh)) %>%
  group_by(par_id,accepted) %>% mutate(sum_neg_llh_all=sum(sum_neg_llh)) %>% 
  pivot_longer(!c(par_id,accepted,broad_age)) %>%
  mutate(broad_age=ifelse(grepl("all",name),"all",broad_age )) %>% select(!name) %>% distinct()

# plot likelihoods: <5y, 65+y, both combined
ggplot(all_likelihoods_hosp,aes(x=accepted,y=value,color=accepted)) + 
  geom_boxplot(fill=NA,size=1/2,width=1/2,outlier.colour=NA) + 
  geom_point(position=position_jitterdodge(seed=1,dodge.width=0.9),alpha=1/2,size=1.5,shape=21) +
  facet_wrap(~broad_age,scales = "free_y") +
  scale_color_manual(values=c("black","red")) + labs(color="accepted parameterisations") + scale_y_log10() +
  xlab("") + ylab("negative log-likelihood") + theme_bw() + standard_theme # +
  # theme(legend.position="none")

# calculate likelihoods for attack rates
results_summ_rejected <- read_csv("simul_output/2e4_parsets/results_summ_all.csv") %>% 
  filter(par_id %in% unique(all_likelihoods_hosp$par_id[!all_likelihoods_hosp$accepted]) )

likelihoods_attackrates <- bind_rows(results_summ_all_reg_dyn %>% mutate(accepted=TRUE),
                                     results_summ_rejected %>% mutate(accepted=FALSE)) %>%
  filter(epi_year<2019) %>% select(c(epi_year,par_id,agegroup,inf_in_seas_AR,attack_rate_perc,accepted)) %>%
  mutate(AR_log_binom_LLH=dbinom(x=estim_attack_rates$RSV_sympt_posit[agegroup],
         size=estim_attack_rates$n_test[agegroup],prob=attack_rate_perc/100,log=T)) %>%
  group_by(par_id,accepted) %>% summarise(AR_neglog_binom_LLH=-sum(AR_log_binom_LLH,na.rm=T))

# plot likelihoods calculated from attack rates
ggplot(likelihoods_attackrates,aes(x=accepted,y=AR_neglog_binom_LLH,color=accepted)) + 
  geom_boxplot(fill=NA,size=1/2,width=1/2,outlier.colour=NA) + 
  geom_point(position=position_jitterdodge(seed=1,dodge.width=0.9),alpha=1/2,size=1.5,shape=21) +
  # facet_wrap(~broad_age,scales="free_y") +
  scale_color_manual(values=c("black","red")) + labs(color="accepted parameterisations") + scale_y_log10() +
  xlab("") + ylab("negative log-likelihood") + theme_bw() + standard_theme # +

# attack rates >500%?!!?
dyn_parsets_main1_317 <- read_csv("simul_output/2e4_parsets/dyn_parsets_main1_317.csv") %>%
    mutate(date=as.Date("2016-09-01")+t-min(t)) %>% filter(date<as.Date("2020-04-01"))

dyn_parsets_main1_317 %>% filter(par_id %in% 21 & date<as.Date("2019-08-01")) %>%
  mutate(year_week=paste0(year(date),"-",week(date))) %>% group_by(year_week,agegroup) %>% 
  summarise(value=sum(value),date=median(date)) %>%
ggplot() + geom_bar(aes(x=date,y=value/1e3),stat="identity") + ylab("thousand cases") + 
  facet_wrap(~agegroup,scales = "free_y") + theme_bw() + standard_theme

# calculate from dynamics
dyn_parsets_main1_317 %>% filter(par_id %in% 21) %>% filter(isoweek(date)>=40 | isoweek(date)<=13) %>%
  mutate(epi_year=ifelse(week(date)>=40,year(date),year(date)-1)) %>% group_by(epi_year,agegroup) %>%
  summarise(sum(value))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# THESE STEPS CAN BE SKIPPED AS RESULTS ARE ALREADY IN repo_data/ ----> 
# GO TO: `CONTINUE: RESULTS: interyear variation`
#
#
# simul_length_yr<-25; n_post_npi_yr<-4; n_core<-20; memory_max <- 8; start_date_dyn_save <- "2016-09-01"; contact_red=0.95
## these 1046 parameters selected by criteria of 1) attack rates 2) seasonal concentration of cases
# partable_filename <- "partable_filtered_AR_seasconc.csv"; n_row <- nrow(read_csv(partable_filename))
## we split parameter table into `n_core` batches and run them in parallel, the sh file will launch jobs
## write the file that will launch jobs
# command_print_runs<-paste0(c("Rscript fcns/write_run_file.R",n_core,n_row,
# simul_length_yr,n_post_npi_yr,contact_red,
# partable_filename,"SAVE sep_qsub_files",start_date_dyn_save,
# agegroup_res,memory_max),collapse=" ")

## write file to run all simulations:
# write.table(paste0("# master start file \n#!/usr/bin/bash \n",command_print_runs,
#   "\nscp batch_run_files/start_batches.sh . \nscp batch_run_files/batch*.sh . \nmodule load R/3.6.3
# sh start_batches.sh \nrm batch*.sh \nrm start_batches.sh",collapse = "\n"),
#   file="master_start.sh",col.names=F,row.names=F,quote=F)
## run calculation (this is for multiple cores) by 
# `sh master_start.sh`
## collect summary stat results:
# nohup Rscript fcns/collect_save_any_output.R FOLDERNAME summ_parsets* results_summ_all.csv keep &
##########################################
## To remove model parameterisations that exhibit irregular patterns (varying from one year to another),
## need to calculate relative difference (see SI Methods) between 2018/19 and 19/20 season, to do this, run:
#
# FOLDERNAME <- "simul_output/parscan/parsets_filtered_1084_50pct_red/"; n_file<-64; mem_max<-4; 
# start_date_calc<-"2018-10-10"; stop_date_calc<-"2020-03-15"; start_week<-42; stop_week<-9
# command_interydiff_calc<-paste0(c("Rscript fcns/calc_interyear_diff_seq.R",FOLDERNAME,start_date_dyn_save,
#                                   start_date_calc,stop_date_calc,start_week,stop_week),collapse=" ")
# system(command_interydiff_calc)
## collect outputs of cumul difference between incidence rates:
# cmd_collect_interyear_diff<-paste0(c("Rscript fcns/collect_save_any_output.R",FOLDERNAME,
#                                       "summ_diff_interyr* summ_diff_interyr_all.csv"),collapse=" ")
# system(cmd_collect_interyear_diff)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# CONTINUE: RESULTS: interyear variation
# (we only take infection levels (1,2,3) with at least 1% attack rate. 
# Eg. 1st infections for adults are negligible, 
# so inter-year differences in these small numbers we don't take into account)
# summ_diff_interyr <- left_join(read_csv(here(foldername,"summ_diff_interyr_all.csv")) %>% 
#       mutate(par_id_sort=as.numeric(factor(par_id))),rsv_age_groups %>% 
#       mutate(agegroup=row_number()) %>% 
#       select(c(agegroup,stationary_popul)),by="agegroup") %>% 
#       mutate(attack_rate=cumul_mean_incid/stationary_popul)
# 
# ## SI Figure 3: CDF of interyear difference in incidence 
# ggplot(right_join(summ_diff_interyr,
#                       summ_diff_interyr %>% 
#                       group_by(agegroup,infection) %>% 
#                       summarise(attack_rate=round(mean(attack_rate,na.rm=T),3)) %>% 
#                       filter(attack_rate>0.01) %>% 
#                       select(c(agegroup,infection)),by=c("agegroup","infection")),
#        aes(sum_rel_diff)) + stat_ecdf(aes(color=factor(infection),group=infection),geom="step") + 
#   facet_wrap(~agegroup,nrow=2,labeller=labeller(agegroup=label_both)) + 
#   geom_vline(xintercept=2/10,linetype="dashed",size=1/2) +
#   scale_x_log10() +
#   xlab("relative difference in incidence") + ylab("CDF") + labs(color="# infection")+
#   theme_bw() + standard_theme + theme(legend.position="top",strip.text=element_text(size=14),
#         legend.title=element_text(size=15),legend.text=element_text(size=14),
#         axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
#         axis.title.x=element_text(size=16))
# # save
# # ggsave(here(foldername,"interyear_difference_cumul_incid_reg_dyn.png"),width=30,height=25,units="cm")
# ##############################################################
# # select the parameter sets with less than x% inter-year variation in last 2 pre-NPI years
# cutoff_intery_var<-20/100
# parsets_regular_dyn <- right_join(summ_diff_interyr,
#                         summ_diff_interyr %>% 
#                         group_by(agegroup,infection) %>% 
#                         summarise(attack_rate=round(mean(attack_rate,na.rm=T),3)) %>% 
#                         filter(attack_rate>0.01) %>% 
#                         select(c(agegroup,infection)),by=c("agegroup","infection")) %>% 
#                         group_by(par_id) %>% 
#                         summarise(score_reg_dyn=sum(sum_rel_diff<cutoff_intery_var)) %>% 
#                         filter(score_reg_dyn==max(score_reg_dyn))
# # parameter sets with regular dynamics (796)
# partable_regular_dyn <- partable %>% 
#   filter(par_id %in% parsets_regular_dyn$par_id)
# nrow(partable_regular_dyn)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# compare hospitalisations from SIMULATIONS to those predicted from 
# (median attack rate in LITERATURE)*(hosp prob estims from Hodgson)
# append age group name column to `results_summ_all` (summary stats results)
results_summ_all <- results_summ_all %>% 
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],levels=rsv_age_groups$agegroup_name))
# convert infections into hospitalisations
simul_hosp <- left_join(results_summ_all %>% 
                  filter(par_id %in% parsets_regular_dyn$par_id & 
                               epi_year==2019),hosp_probabilities,by="agegroup_name") %>% 
                  mutate(hosp_num_SIMUL_inf_tot=prob_hosp_per_infection*inf_tot,
                         hosp_num_SIMUL_inf_seas=prob_hosp_per_infection*inf_in_seas) %>% ungroup() %>% 
                  select(c(agegroup_name,hosp_num_from_per_inf_prob,hosp_num_SIMUL_inf_tot,par_id)) %>% 
                  mutate(agegroup_name=factor(agegroup_name,levels=unique(agegroup_name))) %>% 
                  pivot_longer(!c(agegroup_name,par_id)) %>%
                  mutate(par_id=ifelse(grepl("per_inf",name) & par_id!=min(par_id),NA,par_id))
# SI FIGURE 5: number of hospitalisations from simulations vs LIT estimate 
# from (literature estimate)*(hospit probab per infection)
simul_lit_hosp_comp <- simul_hosp %>% 
        mutate(line_size=ifelse(grepl("per_inf",name),1/5,2),
          name=ifelse(grepl("SIMUL",name),"simulated hospitalisations","literature estimate"))
# plot
ggplot() + 
  geom_jitter(data=simul_lit_hosp_comp %>% filter(grepl("simul",name)),
                       aes(x=agegroup_name,y=ifelse(value>0,value/1e3,NA),color=name),alpha=1/2,size=1/2) + 
  geom_hpline(data=simul_lit_hosp_comp %>% filter(!grepl("simul",name)),
      aes(x=agegroup_name,y=ifelse(value>0,value/1e3,NA),color=name),
      size=1.2,width=0.95,position=position_dodge(width=1)) + scale_color_manual(values=c("black","blue")) +
  geom_vline(xintercept=(0:11)+1/2,linetype="dashed",size=1/2) + 
  xlab("age group (years)") + ylab("hospitalisations (thousands) in epi-year") + labs(fill="",color="") + 
  scale_x_discrete(expand=expansion(0,0)) +
  scale_y_log10(breaks=round(10^seq(-2,2,by=1/4),2),expand=expansion(0.02,0)) + 
  theme_bw() + standard_theme + 
  theme(legend.position="top",legend.text=element_text(size=15),legend.title=element_text(size=15),
                         axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
                         axis.title.x=element_text(size=15),axis.title.y=element_text(size=15))
# save
# ggsave(here::here(foldername,"hospit_per_age_group_simul_lit_compare.png"),width=30,height=25,units="cm")

# SI FIGURE 4: simulated attack rates vs LIT estimates?
attack_rates_simul_LIT <- left_join(
        results_summ_all %>% 
                        filter(par_id %in% parsets_regular_dyn$par_id & epi_year==2019) %>%
                        mutate(inf_tot=inf_tot/rsv_age_groups$stationary_popul[agegroup]) %>% 
                        select(agegroup_name,par_id,inf_tot), 
        data.frame(agegroup_name=rsv_age_groups$agegroup_name,
                  cumul_inf_LIT_ESTIM_median=estim_attack_rates$median_est/100,
                  cumul_inf_LIT_ESTIM_min=estim_attack_rates$min_est/100,
                  cumul_inf_LIT_ESTIM_max=estim_attack_rates$max_est/100),by="agegroup_name") %>% 
      mutate(agegroup_name=factor(agegroup_name,levels=unique(agegroup_name))) %>% 
      pivot_longer(!c(agegroup_name,par_id)) %>% 
      mutate(par_id=ifelse(grepl("LIT_ESTIM",name) & par_id!=min(par_id),NA,par_id),
         categ=ifelse(grepl("inf_tot|inf_in_seas",name),"SIMUL","LIT_estim"),
         name=ifelse(grepl("inf_tot|inf_in_seas",name),paste0(name,"_SIMUL"),"LIT_estim")) %>% 
      filter(!is.na(par_id)) %>%
      mutate(name=case_when(name %in% "inf_tot_SIMUL" ~ "cumulative incidence (SIMULATION)",
                        name %in% "LIT_estim" ~ "cumulative incidence (LITERATURE ESTIMATE)"))
# SI FIGURE 4
ggplot(attack_rates_simul_LIT %>% 
         filter(grepl("SIMUL",name)), 
  aes(x=agegroup_name,y=ifelse(value>0,value,NA)*1e2,color=name)) + 
  geom_hpline(size=1/5,width=0.47,color="red") + 
  geom_hpline(data=attack_rates_simul_LIT %>% 
                filter(grepl("LIT",name)) %>% group_by(agegroup_name) %>% 
                arrange(value) %>% 
                mutate(min_med_max=row_number()),
        aes(linetype=ifelse(min_med_max==2,"solid","dashed")),size=1,width=1,color="black") +
  geom_hpline(data=attack_rates_simul_LIT %>% 
                group_by(agegroup_name) %>% 
                summarise(value=median(value)),
              size=1.5,width=1,position=position_dodge(width=1),color="orange") +
  geom_vline(xintercept=(0:11)+1/2,linetype="dashed",size=1/2) + 
  xlab("Age Group") + ylab("attack rate (%)") +
  labs(color="") + theme_bw() + standard_theme + 
  theme(legend.position="none",axis.text.x=element_text(size=13),
                axis.text.y=element_text(size=13),legend.text=element_text(size=16))
# SAVE
# ggsave(here::here(foldername,"attack_rate_comparison_with_lit.png"),width=25,height=20,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# reduce the two parameters exp_dep and age_dep to their value along their 1st principal component
partable_regular_dyn <- left_join(parsets_regular_dyn %>% 
                                    select(!score_reg_dyn),partable,by="par_id")
pred_pca <- data.frame(predict(object=prcomp(partable_regular_dyn %>% select(c(exp_dep,age_dep))),
                newdata=partable_regular_dyn %>% select(c(exp_dep,age_dep))),
      K_exp=partable_regular_dyn$exp_dep,
      K_age=partable_regular_dyn$age_dep,
      par_id=partable_regular_dyn$par_id)

# SI Figure 2: linear relationship between K_age and K_exp for selected param sets
ggplot(pred_pca %>% 
         pivot_longer(!c(PC1,PC2,par_id)),
       aes(x=PC1,y=value)) + geom_point(aes(color=name)) + 
  geom_smooth(aes(group=name,color=name),fill=NA,method='lm') + 
  scale_x_continuous(breaks=(-(2*3):(2*2))/4,limits=c(-1,1)) +
  xlab(expression(kappa)) + ylab(expression(paste(kappa[exp],", ",kappa[age]))) + labs(color="",fill="") +
  theme_bw() + standard_theme + 
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        legend.title=element_text(size=16),legend.text=element_text(size=16))
# ggsave(here::here(foldername,"exp_age_PCA.png"),width=25,height=20,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# order by 1st principal component of k_age and k_exp
suscept_sel_parsets <- left_join(
  pred_pca %>% rename(exp_dep=K_exp,age_dep=K_age),age_exp_dep_uniqvals,by=c("exp_dep","age_dep")) %>%
  rename(`exposure-dependence`=exp_dep,`age-dependence`=age_dep) %>% 
  mutate(age=factor(rsv_age_groups$agegroup_name[age],
            levels=unique(rsv_age_groups$agegroup_name)),
         PC1_grouped=findInterval(PC1,seq(-1,1,by=1/5)) ) %>%
  group_by(PC1_grouped) %>% 
  mutate(PC1_grouped=round(mean(PC1),1)) %>% 
  group_by(PC1_grouped,age,exp) %>% 
  summarise(mean_val=mean(susc_scaled),
            med_val=median(susc_scaled),
            ci50_low=quantile(susc_scaled,c(0.25,0.75))[1],ci50_up=quantile(susc_scaled,c(0.25,0.75))[2],
            ci95_low=quantile(susc_scaled,c(0.025,0.975))[1],
            ci95_up=quantile(susc_scaled,c(0.025,0.975))[2]) %>% rename(exposure=exp)
# color palette
colorpal <- colorRampPalette(colors=c("blue","grey","red"))(length(unique(suscept_sel_parsets$PC1_grouped)))
##################################################
# facet by 'kappa' parameter: SI FIGURE 3
# library(plyr)
label_parseall <- function(variable, value) {
  plyr::llply(value, function(x) parse(text=paste(variable,x,sep = "==")))}
# plot
ggplot(suscept_sel_parsets %>% 
         rename(kappa=PC1_grouped) %>% ungroup() %>% 
         mutate(kappa_num=as.numeric(factor(as.character(kappa)))) %>% 
          filter(kappa_num %in% c(1,3,5,7,9,10)),
  aes(x=age,color=factor(exposure),group=exposure,fill=factor(exposure))) +
  geom_line(aes(y=med_val),size=1.06) + geom_ribbon(aes(ymin=ci50_low,ymax=ci50_up),color=NA,alpha=0.2) +
  facet_wrap(~kappa,labeller=label_parseall)+labs(color="exposure",fill="exposure") +
  scale_y_log10(breaks=unlist(lapply(seq(-4,0,1/2), function(x) round(10^x,abs(x))))) + 
  theme_bw() + standard_theme + 
  theme(strip.text=element_text(size=18),legend.title=element_text(size=16),
      legend.text=element_text(size=15),legend.position="top",axis.text.x=element_text(size=14),
      axis.text.y=element_text(size=12),axis.title.x=element_text(size=16),
      axis.title.y=element_text(size=18)) + xlab("age group") + ylab(expression(delta[exp]^(age)))
# save
# ggsave(here::here(foldername,"suscept_by_dep_level.png"),width=25,height=20,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT: FIGURE 2 in main text (first panel): entire range of outputs by each parameter

partable_full_linear_kage_kexp <- read_csv(here::here(foldername,"partable_full_linear_kage_kexp.csv"))
# PCA on full param table
pred_pca_FULL <- data.frame(
  predict(object=prcomp(partable_full_linear_kage_kexp %>% 
              select(c(exp_dep,age_dep))),newdata=partable_full_linear_kage_kexp %>% 
              select(c(exp_dep,age_dep))),
            exp_dep=partable_full_linear_kage_kexp$exp_dep,
            age_dep=partable_full_linear_kage_kexp$age_dep,
            par_id=partable_full_linear_kage_kexp$par_id) %>% 
  rename(kappa=PC1) %>%
  group_by(age_dep,exp_dep) %>% 
  summarise(kappa=unique(kappa))
mean_ages_broad_age_group <- data.frame(agegroup_broad=1:3,mean_age=c(0.5,1.5,3.5))
# concatenate with PC1 of age-exp parameter, bin values
results_summ_all_fullscan <- left_join(
  read_csv(here::here(foldername,"results_summ_all_fullscan_relevant.csv")) %>%
    mutate(hosp_tot=inf_tot*hosp_probabilities$prob_hosp_per_infection[agegroup],
           peak_hosp=peak_inf*hosp_probabilities$prob_hosp_per_infection[agegroup]) %>% 
    mutate(agegroup_broad=c("<1y","1-2y","2-5y","5+y")[findInterval(agegroup,c(2,4,7)+1)+1]) %>%
    group_by(agegroup_broad,par_id,epi_year,exp_dep,age_dep,
             seasforc_width_wks,seasforce_peak,R0,omega) %>% 
    summarise(hosp_tot=sum(hosp_tot),
              peak_hosp=sum(peak_hosp),
              inf_tot=sum(inf_tot)) %>% 
    group_by(agegroup_broad,par_id) %>%
    mutate(hosp_tot_norm=hosp_tot/hosp_tot[epi_year==2019],
           peak_hosp_norm=peak_hosp/peak_hosp[epi_year==2019]),
  pred_pca_FULL,by=c("age_dep","exp_dep")) %>%
  mutate(kappa_grouped=findInterval(kappa,seq(-1,1,by=1/5)) ) %>% group_by(kappa_grouped) %>% 
  mutate(kappa_grouped=round(mean(kappa),1)) %>% relocate(c(kappa,kappa_grouped),.after=age_dep) %>%
  group_by(par_id,epi_year) %>% 
  mutate(mean_age=
          sum(mean_ages_broad_age_group$mean_age[
            as.numeric(factor(unique(agegroup_broad)))]*inf_tot/sum(inf_tot))) %>%
  group_by(par_id) %>% 
  mutate(mean_age_shift_month=(mean_age-mean_age[epi_year==2019])*12)

# find median parameter set
median_parset <- results_summ_all_fullscan %>% ungroup() %>% 
  select(c(kappa_grouped,seasforc_width_wks,seasforce_peak,R0,omega)) %>% 
  group_by(kappa_grouped,seasforc_width_wks,seasforce_peak,R0,omega) %>% 
  summarise_all(unique) %>% ungroup() %>%
  summarise_all(quantile,p=0.5,type=1)

# calculate ranges: fix 4 params to the median, 
# calculate the values at the min and max of the parameter that is varied
# take median parameter set, calculate the range in ONE parameter, 
# for ALL param sets and for the selected parsets
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
for (k in 1:ncol(median_parset)){ if (k==1) {
  output_ranges_full_scan<-data.frame(); mean_age_shift_ranges<-data.frame()}
  # all parameterisations
  sel_all <- right_join(
    results_summ_all_fullscan, median_parset[,-k],by=colnames(median_parset)[-k]) %>% 
    group_by(agegroup_broad,epi_year) %>% 
    filter(epi_year==2021)
  #
  x_full_scan <- bind_rows(sel_all %>% 
                            summarise(norm_median=median(hosp_tot_norm),norm_min=min(hosp_tot_norm),
                            norm_max=max(hosp_tot_norm)) %>% 
                            mutate(scan_param=colnames(median_parset)[k],range="full",vartype="cumulative"),
                            sel_all %>% 
                            summarise(norm_median=median(peak_hosp_norm),norm_min=min(peak_hosp_norm),
                            norm_max=max(peak_hosp_norm)) %>% 
                            mutate(scan_param=colnames(median_parset)[k],range="full",vartype="peak"))
  # accepted parameterisations
  sel_sel_parsets <- right_join(
    results_summ_all_fullscan %>% 
      filter(par_id %in% parsets_regular_dyn$par_id),
    median_parset[,-k],by=colnames(median_parset)[-k]) %>% 
  group_by(agegroup_broad,epi_year) %>% 
    filter(epi_year==2021)
  # 
  x_parset_sel <- bind_rows(sel_sel_parsets %>% 
                              summarise(norm_median=median(hosp_tot_norm),
                                  norm_min=min(hosp_tot_norm),
                                  norm_max=max(hosp_tot_norm)) %>% 
                              mutate(scan_param=colnames(median_parset)[k],range="sel",vartype="cumulative"),
                            sel_sel_parsets %>% 
                              summarise(norm_median=median(peak_hosp_norm),
                                    norm_min=min(peak_hosp_norm),norm_max=max(peak_hosp_norm)) %>% 
                              mutate(scan_param=colnames(median_parset)[k],range="sel",vartype="peak"))
  # collect all results (cumul, peak)
  output_ranges_full_scan <- bind_rows(output_ranges_full_scan,x_full_scan,x_parset_sel)
  # mean age shift 
  x_mean_age <- bind_rows(sel_all %>% 
                              summarise(norm_median=median(mean_age_shift_month),
                              norm_min=min(mean_age_shift_month),
                              norm_max=max(mean_age_shift_month)) %>%
                              mutate(scan_param=colnames(median_parset)[k],range="full"),
                        sel_sel_parsets %>% summarise(norm_median=median(mean_age_shift_month),
                            norm_min=min(mean_age_shift_month),
                            norm_max=max(mean_age_shift_month)) %>%
                    mutate(scan_param=colnames(median_parset)[k],range="sel")) %>% 
    group_by(scan_param,range) %>% 
    select(!agegroup_broad) %>% 
    summarise_all(unique)
  # 
  mean_age_shift_ranges <- bind_rows(x_mean_age,mean_age_shift_ranges)
  
  if (k==ncol(median_parset)) {
    output_ranges_full_scan <- output_ranges_full_scan %>% 
                                  relocate(c(scan_param,range,vartype),.after=epi_year) }
}
### END OF LOOP
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT
# load tidybayes library
if (!any(grepl("tidybayes",row.names(installed.packages())))) {
  install.packages("tidybayes")}; library(tidybayes)
# 
df_plot_fullrange <- output_ranges_full_scan %>%
                        mutate(
                          scan_param=case_when(
                              grepl("kappa_grouped",scan_param) ~ "exposure vs age",
                              grepl("seasforce_peak",scan_param) ~ "seasonal forcing (strength)",
                              grepl("seasforc_width_wks",scan_param) ~ "seasonal forcing (width)",
                              grepl("R0",scan_param) ~ "baseline R0",
                              grepl("omega",scan_param) ~ "waning rate"),
         epi_year=case_when(epi_year %in% "2021" ~ "2021-2022")) %>% filter(epi_year %in% "2021-2022")
# plot change in CUMUL and PEAK HOSP
dodge_val=1
ggplot() + 
  geom_interval(data=df_plot_fullrange %>% filter(range %in% "full"),
          aes(x=norm_median,y=agegroup_broad,group=scan_param,xmin=norm_min,xmax=norm_max,
                    color=factor(scan_param)),position=position_dodge(width=dodge_val),alpha=1/3,size=5) +
  geom_interval(data=df_plot_fullrange %>% filter(range %in% "sel"),
          aes(x=norm_median,y=agegroup_broad,group=scan_param,xmin=norm_min,xmax=norm_max,
                    color=factor(scan_param)),position=position_dodge(width=dodge_val),size=10) +
  geom_vpline(data=df_plot_fullrange %>% filter(range %in% "sel"),
          aes(x=norm_median,y=agegroup_broad,group=scan_param),
              position=position_dodge(width=dodge_val),color="black",size=1.5,height=0.185) + # 
  facet_grid(~vartype,scales="free_x") + 
  geom_hline(yintercept=(0:4)+1/2,size=1/2) + 
  geom_vline(xintercept=1,size=1/2,linetype="dashed") + 
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  scale_x_log10(breaks=c(0.3,0.5,0.75,1,1.5,2,3,5,10)) + scale_y_discrete(expand=expansion(0,0)) + 
  xlab("relative hospitalisation risk compared to pre-pandemic years") + ylab("") + labs(color="") +
  theme_bw() + standard_theme + theme(strip.text=element_text(size=18),
        axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),
        legend.text=element_text(size=16),legend.position="top",axis.title.x=element_text(size=16))
  
# save
# subfldr_name<-here::here(foldername,"median_interquant_by_param_value/summary_range/")
# if (!dir.exists(subfldr_name)) {dir.create(subfldr_name)}
# ggsave(paste0(subfldr_name,"cumul_peak_hosp_summary_plot.png"),width=28,height=24,units="cm")

# check entire range across all params:
df_plot_fullrange %>% 
  group_by(agegroup_broad,epi_year,range,vartype) %>% 
   summarise(norm_median=round(100*median(norm_median)-100),norm_min=round(100*min(norm_min)-100),
   norm_max=round(100*max(norm_max)-100)) %>% 
  filter(range %in% "sel" & vartype %in% "cumulative")

# shift in mean age
df_plot_mean_age <- mean_age_shift_ranges %>%
  mutate(scan_param=case_when(grepl("kappa_grouped",scan_param) ~ "exposure vs age",
                              grepl("seasforce_peak",scan_param) ~ "seasonal forcing (strength)",
                              grepl("seasforc_width_wks",scan_param) ~ "seasonal forcing (width)",
                              grepl("R0",scan_param) ~ "baseline R0",
                              grepl("omega",scan_param) ~ "waning rate"),
         epi_year=case_when(epi_year %in% "2021" ~ "2021-2022")) %>% 
  filter(epi_year %in% "2021-2022")
# plot
ggplot() +
  geom_interval(data=df_plot_mean_age %>% filter(range %in% "full"),
        aes(x=norm_median,y=1,group=scan_param,xmin=norm_min,xmax=norm_max,
                    color=factor(scan_param)),position=position_dodge(width=dodge_val),alpha=1/3,size=4) + #
  geom_interval(data=df_plot_mean_age %>% filter(range %in% "sel"),
        aes(x=norm_median,y=1,group=scan_param,xmin=norm_min,xmax=norm_max,
                    color=factor(scan_param)),position=position_dodge(width=dodge_val),size=8) +
  geom_vpline(data=df_plot_mean_age %>% filter(range %in% "sel"),
        aes(x=norm_median,y=1,group=scan_param),
              position=position_dodge(width=dodge_val),color="black",size=1.5,height=0.185) + # 
  geom_hline(yintercept=(0:4)+1/2,size=1/2) + geom_vline(xintercept=0,size=1/2,linetype="dashed") + 
  guides(color=guide_legend(nrow=2,byrow=TRUE)) + 
  scale_x_continuous(breaks=(-1:5)) + scale_y_discrete(expand=expansion(0,0)) + 
  xlab("shift in average age (months)") + ylab("") + labs(color="") +
  theme_bw() + standard_theme + theme(strip.text=element_text(size=18),axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),legend.text=element_text(size=16),legend.position="null",
        axis.title.x=element_text(size=16))
# save
# ggsave(here::here(subfldr_name,"aver_age_hosp_summary_plot_all_par.png"),width=28,height=6,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Analyse accepted simulations by epidemiologic parameters

# cumulative infections per epi_year --> add hospitalisations
results_summ_all_hosp <- left_join(
  left_join(
      left_join(results_summ_all %>% 
                  filter(par_id %in% partable_regular_dyn$par_id), 
                  hosp_probabilities %>% 
                  select(c(agegroup_name,prob_hosp_per_infection)),by="agegroup_name"),
      pred_pca %>% select(par_id,PC1),by="par_id"),
  partable_regular_dyn %>% 
    select(par_id,omega) %>% 
    mutate(omega=1/omega) %>% 
    rename(waning=omega), by="par_id") %>% 
  mutate(hosp_tot=prob_hosp_per_infection*inf_tot,
         hosp_seas=prob_hosp_per_infection*inf_in_seas) %>% 
  select(!prob_hosp_per_infection) %>% 
  relocate(agegroup_name,.after=agegroup) %>% 
  mutate(agegroup_broad=c("<1y","1-2y","2-5y","5+y")[findInterval(agegroup,c(2,4,7)+1)+1]) %>% 
  relocate(agegroup_broad,.after=agegroup_name) %>% 
  relocate(c(inf_tot,inf_in_seas),.before=hosp_tot) %>% 
  relocate(par_id,.after=epi_year) %>% 
  relocate(c(exp_dep,age_dep,PC1,seasforce_peak,R0,waning,seasforc_width_wks),.after=par_id)
# write_csv(results_summ_all_hosp,here::here(foldername,"results_summ_all_hosp.csv"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# outcomes by individual parsets # plot sums: pre-NPI, NPI+1 (2021/22), NPI+2 (2022/23)
# LOAD
# parsets_broad_age_groups <- read_csv(here::here(foldername,"parsets_broad_age_groups.csv"))
parsets_broad_age_groups <- results_summ_all_hosp %>% 
  group_by(par_id,epi_year,agegroup_broad) %>% 
  summarise(PC1=unique(PC1),
            seasforce_peak=unique(seasforce_peak),
            waning=round(unique(waning)),
            seasforc_width_wks=unique(seasforc_width_wks),
            R0=unique(R0),
            inf_tot=sum(inf_tot),
            inf_in_seas=sum(inf_in_seas),
            hosp_tot=sum(hosp_tot),
            hosp_seas=sum(hosp_seas),
            seas_share=mean(seas_share),
            attack_rate_perc=mean(attack_rate_perc),
            max_incid_week=mean(max_incid_week)) %>% 
  pivot_longer(!c(par_id,epi_year,agegroup_broad,PC1,seasforce_peak,waning,seasforc_width_wks,R0)) %>%
  group_by(par_id,agegroup_broad,name) %>% 
  summarise(epi_year,PC1,seasforce_peak,waning,seasforc_width_wks,R0,value,
    value_norm=ifelse(epi_year==2019,1,ifelse(name %in% c("max_incid_week","seas_share","attack_rate_perc"),
    value-value[epi_year==2019],value/value[epi_year==2019]))) %>% 
  relocate(name,.after=R0) %>%
  mutate(age_exp_par_bins=findInterval(PC1,seq(-1,1,by=1/5))) %>% 
  group_by(age_exp_par_bins) %>% 
  mutate(age_exp_par_bins=round(mean(PC1),1)) %>% 
  relocate(age_exp_par_bins,.after=PC1)
# write_csv(parsets_broad_age_groups,here::here(foldername,"parsets_broad_age_groups.csv"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# summary plot (median, interquartile range)
summ_broad_age_groups <- parsets_broad_age_groups %>% 
  mutate(value_norm=ifelse(name %in% c("seas_share"),
         value_norm*100,value_norm)) %>% 
  group_by(agegroup_broad,epi_year,name) %>% 
  summarise(mean=mean(value_norm),median=median(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
            ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% 
  filter(epi_year>2019)
# write_csv(summ_broad_age_groups,here::here(foldername,"summ_broad_age_groups.csv"))
# summ_broad_age_groups <- read_csv(here::here(foldername,"summ_broad_age_groups.csv"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# segment age_exp_dep into x values, calculate summary statistics for each param value
summ_broad_age_groups_byvalue <- parsets_broad_age_groups %>% 
  select(!c(PC1,value)) %>% 
  mutate(waning=round(waning)) %>%
  relocate(age_exp_par_bins,.after=R0) %>% 
  rename(varname=name) %>% 
  pivot_longer(!c(par_id,agegroup_broad,epi_year,varname,value_norm)) %>% 
  rename(parname=name,parvalue=value) %>% 
  relocate(c(varname,value_norm),.after=parvalue) %>%
  mutate(value_norm=ifelse(varname %in% c("seas_share"),
        value_norm*100,value_norm)) %>%
  group_by(agegroup_broad,epi_year,parname,parvalue,varname) %>% 
  summarise(mean=mean(value_norm),median=median(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2],
            ci95_low=quantile(value_norm,c(0.025,0.975))[1],
            ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% 
  filter(epi_year>2019)
# write_csv(summ_broad_age_groups_byvalue,here::here(foldername,"summ_broad_age_groups_byvalue.csv"))
# summ_broad_age_groups_byvalue <- read_csv(here::here(foldername,"summ_broad_age_groups_byvalue.csv"))
##############################################################
# PLOTS of summary statistics across ALL parameter values
# the plots below include: Figure 3,4,5 and SI Figure 5
sel_vars <- c("attack rate","in-season hospitalisations","cumulative hospitalisations",
              "in-season infections","cumulative infections","% cases in-season")
for (k_plot in c(3,5)) {
  for (k_age_excl in 1:2) {
    dodge_val=0.9
    ylab_tag <- ifelse(grepl("season peak|% cases|attack",sel_vars[k_plot]),
                       " (change from 2019 level)"," (relative to pre-NPI)")
    # subset data to plot
    df_plot <- summ_broad_age_groups %>% 
      mutate(name=case_when(
                  grepl("attack",name) ~ "attack rate",
                  grepl("hosp_seas",name) ~ "in-season hospitalisations", 
                  grepl("hosp_tot",name) ~ "cumulative hospitalisations", 
                  grepl("inf_in_seas",name) ~ "in-season infections",
                  grepl("inf_tot",name) ~ "cumulative infections",
                  grepl("max_incid_week",name) ~ "season peak (calendar week)",
                  grepl("seas_share",name) ~ "% cases in-season")) %>%
      filter(epi_year>2020 & (name %in% sel_vars[k_plot]))
    if (k_age_excl>1) {
      df_plot <- df_plot %>% 
        filter(!agegroup_broad %in% "5+y"); size_adj<-1.33 } else {
          size_adj<-1
          }
    p <- ggplot(df_plot, 
      aes(x=factor(epi_year),color=factor(agegroup_broad),group=agegroup_broad)) + 
      geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),
                     position=position_dodge(width=dodge_val),
                     alpha=0.3,size=12*size_adj,show.legend=FALSE) +
      geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),
                     position=position_dodge(width=dodge_val),
                     alpha=0.6,size=12*size_adj)+
      geom_hpline(aes(y=median),position=position_dodge(width=dodge_val),
                  width=size_adj/6,size=0.8,color="black") +
      geom_vline(xintercept=(0:4)+1/2,size=1/2) + 
      scale_x_discrete(expand=expansion(0,0)) + 
      xlab("") + ylab(paste0(sel_vars[k_plot],ylab_tag)) + labs(color="age groups") + 
      theme_bw() + standard_theme + theme(strip.text=element_text(size=15),
          legend.text=element_text(size=15),
          axis.text.x=element_text(size=13),axis.text.y=element_text(size=13),
          legend.title=element_text(size=15),panel.grid.major.x=element_blank())
    # 
    if (grepl("season peak|% cases|attack",sel_vars[k_plot])) {
      if (grepl("season peak",sel_vars[k_plot])) {
        break_vals <- (-5:15)*10} else {break_vals <- (-10:10)*10}
      p <- p + geom_hline(yintercept=0,linetype="dashed",size=1/2) + scale_y_continuous(breaks=break_vals) 
                    } else {
        p <- p + scale_y_log10(breaks=c(0.1,1/4,1/2,3/4,1,5/4,3/2,7/4,2,2.5,3,5)) + 
          geom_hline(yintercept=1,linetype="dashed",size=1/2)
        }
    p
    # save
    sel_var_filename <- gsub("%","share",gsub("_calendar_week","",
                              gsub("\\(|\\)","",gsub("-|\\s","_",sel_vars[k_plot]))))
    subfldr_name<-"median_interquant_all_collapsed/"
    if (!dir.exists(here::here(foldername,subfldr_name))) {
      dir.create(here::here(foldername,subfldr_name)) }
    if (!dir.exists(here::here(foldername,subfldr_name,"above_5y/"))){
      dir.create(here::here(foldername,subfldr_name,"above_5y/")) }
    ggsave(here::here(foldername,subfldr_name,
                      paste0(ifelse(k_age_excl==1,"above_5y/",""),"summ_stats_relative_2019_",
                      paste0(sel_var_filename,collapse="_"),".png") ),
            width=25,height=20,units="cm")
    print(sel_vars[k_plot])
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT summary stats disaggregated by parameter VALUES
# the figures produced by this loop include: Figure 3A-B (main text) and SI Fig 6-11

sel_vars <- c("attack rate","cumulative hospitalisations", 
              "cumulative infections","% cases in-season")[c(2,3)]
sel_pars <- c("age_exp_par_bins","R0","seasforc_width_wks","seasforce_peak","waning")
start_year<-2021
for (k_plot_var in 1:length(sel_vars)) {
  for (k_plot_par in 1:length(sel_pars)) {
    sel_par <- sel_pars[k_plot_par]; dodge_val=1
    ylab_tag <- ifelse(grepl("season peak|cases in-season|attack",sel_vars[k_plot_var]),
                       " (change from 2019 level)"," (relative to pre-NPI)")
    df_plot <- summ_broad_age_groups_byvalue %>% 
      mutate(varname=case_when(grepl("attack",varname) ~ "attack rate",
                     grepl("hosp_tot",varname) ~ "cumulative hospitalisations",
                     grepl("inf_tot",varname) ~ "cumulative infections", 
                     grepl("max_incid_week",varname) ~ "season peak (calendar week)",
                     grepl("seas_share",varname) ~ "% cases in-season")) %>% 
      filter(epi_year>=start_year & epi_year<2024 & (varname %in% sel_vars[k_plot_var]) & 
           (parname %in% sel_par) & (!agegroup_broad %in% "5+y") ) %>%
      mutate(parname=case_when(grepl("age_exp_par_bins",parname) ~ "exposure (-1) <-> age (1)", 
                               grepl("seasforc_width_wks",parname) ~ "season width (weeks)", 
                               grepl("seasforce_peak",parname) ~ "seasonal forcing (above baseline)",
                               grepl("R0",parname) ~ "R0 (baseline)",
                               grepl("waning",parname) ~ "waning (days)"))
    n_par_value <- length(unique(df_plot$parvalue))
    # colour palette
    if (!grepl("age_exp_par_bins",sel_par)){ 
      colorpal=colorRampPalette(colors=c("orange","red"))(n_par_value)} else  {
      colorpal=colorRampPalette(colors=c("blue","grey","red"))(n_par_value) 
      }
    p <- ggplot(df_plot,
                aes(x=factor(epi_year),color=factor(parvalue),group=parvalue)) + 
      facet_wrap(~agegroup_broad,scales="free_y") + 
      geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge(width=dodge_val),
                     alpha=0.3,size=28/n_par_value) +
      geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),position=position_dodge(width=dodge_val),
                     alpha=0.6,size=24/n_par_value) +
      geom_hpline(aes(y=median),position=position_dodge(width=dodge_val),width=(1/n_par_value)*0.75,
                  size=0.8,color="black") + 
      geom_vline(xintercept=(0:4)+1/2,size=1/5) + labs(color=unique(df_plot$parname)) +
      scale_x_discrete(expand=expansion(0.02,0)) + 
      xlab("") + ylab(paste0(sel_vars[k_plot_var],ylab_tag)) + 
      theme_bw() + standard_theme + theme(strip.text=element_text(size=15),axis.text.x=element_text(size=13),
          axis.text.y=element_text(size=12),legend.text=element_text(size=11),legend.title=element_text(size=12),
          legend.position=ifelse(grepl("expos|forcing",unique(df_plot$parname)),"bottom","right")) + 
    scale_color_manual(values=colorpal)
      
    if (grepl("season peak|cases in-season|attack",sel_vars[k_plot_var])) { #  # + scale_y_log10() 
      if (grepl("season peak",sel_vars[k_plot_var])) {
        break_vals <- (-5:15)*10} else { break_vals <- (-10:10)*10 }
      p <- p + geom_hline(yintercept=0,linetype="dashed",size=1/2) + scale_y_continuous(breaks=break_vals) } else {
        p <- p + geom_hline(yintercept=1,linetype="dashed",size=1/2)}; p
    # create filename and folders
    sel_var_filename <- gsub("%","share",gsub("_calendar_week","",
                                            gsub("\\(|\\)","",gsub("-|\\s","_",sel_vars[k_plot_var]))))
    if (!dir.exists("median_interquant_by_param_value")) {
          dir.create(here::here(foldername,"median_interquant_by_param_value"))}
    subfldr_name <- "median_interquant_by_param_value/cumul/"
    if (!dir.exists(here::here(foldername,subfldr_name))) {
          dir.create(here::here(foldername,subfldr_name))}
    # save
    ggsave(here::here(foldername,paste0(subfldr_name,"summ_stats_relative_2019_",
                          ifelse(start_year==2020,"incl2020_",""),
                          paste0(sel_var_filename,collapse="_"),"_",sel_par,".png")),
           width=28,height=16,units="cm")
    print(paste0(c(sel_vars[k_plot_var],sel_pars[k_plot_par]),collapse=", "))
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Calculate the average age of infection (by indiv parsets)

parsets_mean_age_inf <- left_join(
  results_summ_all_hosp,rsv_age_groups %>% 
          select(agegroup_name,mean_age_weighted),by="agegroup_name") %>% 
  group_by(par_id,epi_year) %>% 
  summarise(PC1=unique(PC1),
            seasforce_peak=unique(seasforce_peak),
            waning=round(unique(waning)),
            seasforc_width_wks=unique(seasforc_width_wks),
            R0=unique(R0), 
    mean_age_inf_tot_under_5=sum(12*inf_tot[agegroup<=7]*mean_age_weighted[agegroup<=7]/sum(inf_tot[agegroup<=7])),
    mean_age_inf_seas_under_5=
      sum(12*inf_in_seas[agegroup<=7]*mean_age_weighted[agegroup<=7]/sum(inf_in_seas[agegroup<=7])),
    mean_age_hosp_tot_under_5=
      sum(12*hosp_tot[agegroup<=7]*mean_age_weighted[agegroup<=7]/sum(hosp_tot[agegroup<=7])),
    mean_age_hosp_seas_under_5=
      sum(12*hosp_seas[agegroup<=7]*mean_age_weighted[agegroup<=7]/sum(hosp_seas[agegroup<=7]))) %>%
  pivot_longer(!c(par_id,epi_year,PC1,seasforce_peak,waning,seasforc_width_wks,R0)) %>%
  group_by(par_id,name) %>% 
  summarise(epi_year,PC1,seasforce_peak,waning,seasforc_width_wks,R0,value,
      value_norm=ifelse(epi_year==2019,1,value-value[epi_year==2019])) %>% 
  relocate(name,.after=R0) %>%
  mutate(age_exp_par_bins=findInterval(PC1,seq(-1,1,by=1/5))) %>% 
  group_by(age_exp_par_bins) %>% 
  mutate(age_exp_par_bins=round(mean(PC1),1)) %>% 
  relocate(age_exp_par_bins,.after=PC1)
###########################################################
# summary plot (median, interquartile range)
summ_mean_age_infs <- parsets_mean_age_inf %>% 
  group_by(epi_year,name) %>% 
  summarise(mean=mean(value_norm),
            median=median(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2],
            ci95_low=quantile(value_norm,c(0.025,0.975))[1],
            ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% filter(epi_year>2019)
###########################################################
# segment age_exp_dep into x values, calculate summary statistics for each param value
summ_mean_age_inf_byvalue <- parsets_mean_age_inf %>% 
  select(!c(PC1,value)) %>% 
  mutate(waning=round(waning)) %>%
  relocate(age_exp_par_bins,.after=R0) %>% rename(varname=name) %>% 
  pivot_longer(!c(par_id,epi_year,varname,value_norm)) %>%
  rename(parname=name,parvalue=value) %>% relocate(c(varname,value_norm),.after=parvalue) %>%
  mutate(value_norm=ifelse(varname %in% c("seas_share"),
         value_norm*100,value_norm)) %>%
  group_by(epi_year,parname,parvalue,varname) %>% 
  summarise(mean=mean(value_norm),
            median=median(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2],
            ci95_low=quantile(value_norm,c(0.025,0.975))[1],
            ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% 
  filter(epi_year>2019)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PLOT summary stats (mean age inf)
# Figure 3C
sel_vars <- c("mean_age_hosp_seas_under_5","mean_age_hosp_tot_under_5")
dodge_val=0.9
for (k_plot in 1:length(sel_vars)) {
  df_plot <- summ_mean_age_infs %>% 
    filter(epi_year>2020 & (name %in% sel_vars[k_plot])) %>% 
    mutate(name=case_when(grepl("mean_age_hosp_seas_under_5",name) ~ "mean age of in-season hospitalisations", 
                          grepl("mean_age_hosp_tot_under_5",name) ~ "mean age of cumulative hospitalisations"))
  # plot
  ggplot(df_plot,
    aes(x=factor(epi_year))) + 
    geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge(width=dodge_val),
                   alpha=0.3,size=12,show.legend=FALSE) +
    geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),position=position_dodge(width=dodge_val),alpha=0.6,size=12) + 
    geom_hpline(aes(y=median),position=position_dodge(width=dodge_val),width=1/6,size=1.25,color="black") + # 
    geom_vline(xintercept=(0:4)+1/2,size=1/2) + 
    geom_hline(yintercept=0,linetype="dashed",size=1/2) +
    scale_x_discrete(expand=expansion(0,0)) + scale_y_continuous(breaks=(-10:15)) +
    xlab("") + ylab(paste0(unique(df_plot$name)," (change from 2019 in months)")) + 
    theme_bw() + standard_theme + 
    theme(strip.text=element_text(size=15),legend.text=element_text(size=15),
        axis.text.x=element_text(size=13),axis.text.y=element_text(size=13),legend.title=element_text(size=15))
  # create folders
  subfldr_name <- "median_interquant_all_collapsed/mean_age_all_param/"
  if (!dir.exists(here::here(foldername,subfldr_name))) {
    dir.create(here::here(foldername,subfldr_name))}
  # save
  ggsave(here::here(foldername,paste0(subfldr_name,gsub("-|\\s","_",sel_vars[k_plot]),".png") ),
         width=28,height=22,units="cm")
  print(sel_vars[k_plot])
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Plot of shift in average age of inf toggled by param values
# Figure 8, SI Fig 10, 11
sel_vars <- c("mean_age_hosp_tot_under_5") # "mean_age_hosp_seas_under_5",
sel_pars <- c("age_exp_par_bins","R0","seasforc_width_wks","seasforce_peak","waning")
dodge_val=1
for (k_plot_var in 1:length(sel_vars)) {
  for (k_plot_par in 1:length(sel_pars)) {
    sel_par <- sel_pars[k_plot_par];
    df_plot <- summ_mean_age_inf_byvalue %>% 
      filter(epi_year>2020 & epi_year<2024 & (varname %in% sel_vars[k_plot_var]) & (parname %in% sel_par)) %>%
      mutate(varname=case_when(
              grepl("mean_age_hosp_seas_under_5",varname) ~ "mean age of in-season hospitalisations",
              grepl("mean_age_hosp_tot_under_5",varname) ~ "change in mean age of hosp.")) %>% 
      mutate(parname=case_when(
              grepl("age_exp_par_bins",parname) ~ "exposure (-1) <-> age (1)", 
              grepl("seasforc_width_wks",parname) ~ "season width (weeks)",
              grepl("seasforce_peak",parname) ~ "seasonal forcing (above baseline)",
              grepl("R0",parname) ~ "R0 (baseline)",
              grepl("waning",parname) ~ "waning (days)"))
    n_par_value <- length(unique(df_plot$parvalue))
    # colour palette
    if (!grepl("age_exp_par_bins",sel_par)){
      colorpal=colorRampPalette(colors=c("orange","red"))(n_par_value)} else  {
      colorpal=colorRampPalette(colors=c("blue","grey","red"))(n_par_value) }
    # plot
    ggplot(df_plot,
      aes(x=factor(epi_year),color=factor(parvalue),group=parvalue)) + 
      geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),position=position_dodge(width=dodge_val),
                     alpha=0.6,size=85/n_par_value) + #
      geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge(width=dodge_val),
                     alpha=0.3,size=65/n_par_value) + #
      geom_hpline(aes(y=median),position=position_dodge(width=dodge_val),
                  width=(1/n_par_value)*0.85,size=1.25,color="black") + 
      geom_vline(xintercept=(0:4)+1/2,size=1/5) + labs(color=unique(df_plot$parname)) + 
      geom_hline(yintercept=0,linetype="dashed",size=1/2) +
      xlab("") + ylab(paste0(unique(df_plot$varname)," (from 2019)")) + 
      scale_x_discrete(expand=expansion(0.02,0)) + scale_y_continuous(breaks=-10:10) +
      theme_bw() + standard_theme + theme(strip.text=element_text(size=15),axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),axis.title.y=element_text(size=18),legend.text=element_text(size=15),
        legend.title=element_text(size=17),legend.position=ifelse(grepl("expos|forcing",
        unique(df_plot$parname)),"bottom","right")) + scale_color_manual(values=colorpal)
    # file name, folders
    subfldr_name <- "median_interquant_by_param_value/mean_age_by_paramval/"
    if (!dir.exists(here::here(foldername,subfldr_name))) {
      dir.create(here::here(foldername,subfldr_name))}
    # save
    ggsave(here::here(foldername,paste0(subfldr_name,sel_vars[k_plot_var],"_",sel_par,".png")),
           width=35,height=18,units="cm")
    print(paste0(c(sel_vars[k_plot_var],sel_pars[k_plot_par]),collapse=", "))
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Calculate peak cases/hospitalisations and length of season (defined as cases > some threshold value)

# unzip dynamics file
unzip(zipfile=here::here(foldername,"dyn_all_parsets_broad_age.zip"),exdir=here::here(foldername))
# LOAD dynamics of simulations
dyn_all_parsets_broad_age <- read_csv(here::here(foldername,"dyn_all_parsets_broad_age.csv"))
# delete csv: 
# unlink(here::here(foldername,"dyn_all_parsets_broad_age.csv"))

##############################################################
# plot pre-pandemic dynamics for all parameter sets (Figure 1): all age groups <5y, compare to data

# medians across simulated (accepted) parameter sets
median_weekly_pred <- dyn_all_parsets_broad_age %>% 
  filter(agegroup_broad %in% c("<1y","1-2y","2-5y")) %>%
  mutate(sel_par=ifelse(par_id %in% parsets_regular_dyn$par_id,"selected","discarded")) %>%
  filter(sel_par %in% "selected") %>% 
  group_by(date,par_id) %>% 
  summarise(incid_hosp=sum(incid_hosp)) %>% 
  mutate(week_number=week(date),
         incid_hosp_per100k=1e5*incid_hosp/sum(rsv_age_groups$stationary_popul[1:7])) %>% 
  group_by(date) %>% 
  summarise(incid_hosp_per100k=median(incid_hosp_per100k),
            incid_hosp_med_val=median(incid_hosp))
######
# scaling: take mean value per week for years with data in SARIwatch, take the maximal weekly value -> 
# divide the maximal value of the median simulated weekly peak by this value
# data from SARI_watch
SARI_watch_all_hosp <- read_csv(here::here(foldername,"SARI_watch_all_hosp.csv"))
# this data is under-reported (not everyone hospitalised is tested for RSV), so we scale by the ratio 
scale_fact <- max((median_weekly_pred %>% 
                     filter(date<as.Date("2020-04-01")) %>% 
                     mutate(week_number=week(date)) %>% 
                     group_by(week_number) %>% 
                     summarise(incid_hosp_per100k=median(incid_hosp_per100k)))$incid_hosp_per100k)/
                                max((SARI_watch_all_hosp %>% 
                                    filter(epi_year<2021) %>% 
                                    group_by(week_number) %>% 
                                    summarise(rate_per_100000=median(rate_per_100000)))$rate_per_100000)
# scale hospit admissions in data by constant
SARI_watch_all_hosp_scaled <- left_join(
  SARI_watch_all_hosp,median_weekly_pred %>% 
    mutate(week_number=week(date),epi_year=ifelse(week_number>26,year(date),year(date)-1)),
    by=c("week_number","epi_year")) %>% 
    mutate(rate_per_100000_scaled=rate_per_100000*scale_fact) %>% 
    select(!c(incid_hosp_per100k,incid_hosp_med_val))
# overlay plot
date_limits <- as.Date(c("2017-09-15","2020-04-01"))
ggplot() + 
  geom_line(data=dyn_all_parsets_broad_age %>% 
                    filter(date>date_limits[1] & date<date_limits[2] & 
                             agegroup_broad %in% c("<1y","1-2y","2-5y")) %>%
                    group_by(date,par_id) %>% 
                    summarise(incid_hosp=sum(incid_hosp)) %>% 
                    mutate(week_number=week(date),
                           incid_hosp_per100k=1e5*incid_hosp/sum(rsv_age_groups$stationary_popul[1:7]),
                           sel_par=ifelse(par_id %in% parsets_regular_dyn$par_id,"selected","discarded")),
      aes(x=date,y=incid_hosp_per100k,group=par_id,alpha=sel_par,color=sel_par)) +
  scale_color_manual(values=c("grey","blue")) + scale_alpha_manual(values=c(1/10,1/10)) +
  geom_point(data=SARI_watch_all_hosp_scaled %>% 
               filter(date>date_limits[1] & date<date_limits[2]),
      aes(x=date,y=rate_per_100000_scaled)) + # overlay data # ,shape=21,fill=NA,size=2
  geom_line(data=median_weekly_pred %>% 
              filter(date>date_limits[1] & date<date_limits[2]),
      aes(x=date,y=incid_hosp_per100k),color="red",linetype="dashed",size=1.02) + # median simulation
  theme_bw() + standard_theme + 
  xlab("") + ylab("weekly hospitalisations <5y per 100.000 persons") + labs(alpha="",color="") +
  scale_x_date(expand=expansion(1/100,0),date_breaks="2 months") + 
  scale_y_continuous(expand=expansion(0.01,0)) +
  theme(strip.text=element_text(size=15),axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=12),legend.text=element_text(size=11),
        legend.title=element_text(size=12),legend.position="null")
# SAVE
# ggsave(here::here(foldername,
#           "prepandemic_dynamics_all_sel_pars_hosp_per_popul_SARIwatch_under5_median_simul.png"),
#           width=28,height=16,units="cm")

##############################################################
# filter out discarded parameterisations 

# summary by year
epi_year_week_start<-week(as.Date("2020-06-01")); comp_year<-2018
# summary by year
summ_dyn_max_incid_seas_length <- dyn_all_parsets_broad_age %>% 
  filter(par_id %in% parsets_regular_dyn$par_id & 
           date>as.Date("2018-06-01") & date<as.Date("2024-05-01")) %>% 
  mutate(epi_year=ifelse(week(date)>=epi_year_week_start,year(date),year(date)-1)) %>% 
  rename(incid_case=value) %>% 
  pivot_longer(!c(agegroup_broad,date,par_id,epi_year)) %>% 
  group_by(epi_year,par_id,agegroup_broad,name) %>%
  mutate(seas_tot_2018=ifelse(epi_year==2018,sum(value),NA)) %>% 
  group_by(par_id,agegroup_broad,name) %>%
  mutate(seas_tot_2018=min(seas_tot_2018,na.rm=T)) %>% ungroup() %>% 
  mutate(above_baseline=value>seas_tot_2018/52) %>%
  group_by(epi_year,par_id,agegroup_broad,name) %>% 
  summarise(max_value=max(value),sum_value=sum(value),
     seas_length_wk=sum(above_baseline)) %>% filter(epi_year>2017)
# left-join with partable to have input parameters
parsets_max_incid_seas_length <- left_join(
  left_join(summ_dyn_max_incid_seas_length,pred_pca %>% 
              select(par_id,PC1),by="par_id"),
  partable_regular_dyn %>% 
        select(par_id,omega,seasforc_width_wks,R0,seasforce_peak) %>%
        mutate(omega=1/omega) %>% 
        rename(waning=omega), by="par_id") %>% 
  relocate(c(name,max_value,sum_value,seas_length_wk),.after=seasforce_peak) %>% 
  rename(varname=name) %>% 
  pivot_longer(!c(epi_year,par_id,agegroup_broad,PC1,waning,seasforc_width_wks,R0,seasforce_peak,varname)) %>%
  rename(vartype=name) %>% 
  group_by(par_id,agegroup_broad,PC1,waning,seasforc_width_wks,R0,seasforce_peak,varname,vartype) %>% 
  mutate(value_norm=ifelse(grepl("seas_length_wk",vartype),
                           value-value[epi_year==comp_year],
                           value/value[epi_year==comp_year])) %>%
  mutate(age_exp_par_bins=findInterval(PC1,seq(-1,1,by=1/5))) %>% 
  relocate(age_exp_par_bins,.after=PC1) %>%
  group_by(age_exp_par_bins) %>% 
  mutate(age_exp_par_bins=round(mean(PC1),1))
##############################################################
# summary statistics collapse across all param values
summ_max_incid_seas_length <- parsets_max_incid_seas_length %>% 
  group_by(agegroup_broad,epi_year,varname,vartype) %>% 
  summarise(mean=mean(value_norm),median=median(value_norm),
    ci50_low=quantile(value_norm,c(0.25,0.75))[1],ci50_up=quantile(value_norm,c(0.25,0.75))[2],
    ci95_low=quantile(value_norm,c(0.025,0.975))[1],ci95_up=quantile(value_norm,c(0.025,0.975))[2])
##############################################################
# calculate summary statistics for each param value
summ_max_incid_seas_length_byvalue <- parsets_max_incid_seas_length %>% mutate(waning=round(waning)) %>% 
  select(!c(PC1,value)) %>% 
  relocate(age_exp_par_bins,.after=R0) %>% 
  pivot_longer(!c(par_id,agegroup_broad,epi_year,varname,vartype,value_norm)) %>% 
  rename(parname=name,parvalue=value) %>% relocate(c(varname,value_norm),.after=parvalue) %>%
  group_by(agegroup_broad,epi_year,parname,parvalue,varname,vartype) %>% 
  summarise(mean=mean(value_norm),median=median(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2],
            ci95_low=quantile(value_norm,c(0.025,0.975))[1],
            ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% filter(epi_year>2019)

##############################################################
# plot PEAK VALUE of cases/hospitalisations (+ season length), disaggregated by parameter values

# Figure 3B, SI 7-10
start_year<-2021; end_year<-2023
sel_vars <- c("incid_case","incid_hosp")
sel_pars <- c("age_exp_par_bins","R0","seasforc_width_wks","seasforce_peak","waning")
sel_vartypes <- c("max_value","seas_length_wk")[1]
for (k_plot_var in 1:length(sel_vars)) {
  for (k_plot_par in 1:length(sel_pars)) {
    for (k_plot_vartype in 1:length(sel_vartypes)) {
      sel_vartype <- sel_vartypes[k_plot_vartype]; sel_par <- sel_pars[k_plot_par]; dodge_val=1
      # sel data to plot
      df_plot <- summ_max_incid_seas_length_byvalue %>% 
                    filter(epi_year>=start_year & epi_year<=end_year & 
                             (varname %in% sel_vars[k_plot_var]) & (parname %in% sel_par) &
                             (vartype %in% sel_vartype) & (!agegroup_broad %in% "5+y")) %>%
        mutate(varname=case_when(grepl("incid_case",varname) ~ "cases",
                                 grepl("incid_hosp",varname) ~ "hospitalisations"),
               vartype=case_when(grepl("max_value",vartype) ~ "peak", 
                                 grepl("seas_length_wk",vartype) ~ "above baseline (weeks)")) %>% 
        mutate(parname=case_when(grepl("age_exp_par_bins",parname) ~ "exposure (-1) <-> age (1)", 
                                 grepl("seasforc_width_wks",parname) ~ "season width (weeks)", 
                                 grepl("seasforce_peak",parname) ~ "seasonal forcing (above baseline)",
                                 grepl("R0",parname) ~ "R0 (baseline)",
                                 grepl("waning",parname) ~ "waning (days)"))
      ylab_tag <- paste0(paste0(unique(df_plot$vartype)," ",unique(df_plot$varname)),
                              ifelse(grepl("above",unique(df_plot$vartype)),
                                     " (change from 2019 level)"," (relative to pre-NPI)"))
      n_par_value <- length(unique(df_plot$parvalue))
      # colour palette
      if (!grepl("age_exp_par_bins",sel_par)){
              colorpal=colorRampPalette(colors=c("orange","red"))(n_par_value)} else  {
              colorpal=colorRampPalette(colors=c("blue","grey","red"))(n_par_value) }
      p <- ggplot(df_plot,
                  aes(x=factor(epi_year),color=factor(parvalue),group=parvalue)) + 
        facet_wrap(~agegroup_broad,scales="free_y") + 
        geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),
                       position=position_dodge(width=dodge_val),alpha=0.6,size=24/n_par_value) +
        geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),
                       position=position_dodge(width=dodge_val),alpha=0.3,size=24/n_par_value) +
        geom_hpline(aes(y=median),
                    position=position_dodge(width=dodge_val),
                    width=(1/n_par_value)*0.75,size=0.8,color="black") + 
        geom_vline(xintercept=(0:4)+1/2,size=1/5) + 
        scale_x_discrete(expand=expansion(0.02,0)) + 
        xlab("") + ylab(ylab_tag) + labs(color=unique(df_plot$parname)) +
        theme_bw() + standard_theme + 
        theme(strip.text=element_text(size=15),axis.text.x=element_text(size=13),
              axis.text.y=element_text(size=12),legend.text=element_text(size=11),
              legend.title=element_text(size=12),
              legend.position=ifelse(grepl("expos|forcing",unique(df_plot$parname)),"bottom","right")) + 
        scale_color_manual(values=colorpal)
      if (grepl("season peak|cases in-season|attack",sel_vars[k_plot_var])) {
        if (grepl("season peak",sel_vars[k_plot_var])) {
          break_vals <- (-5:15)*10} else {break_vals <- (-10:10)*10 }
        p <- p + geom_hline(yintercept=0,linetype="dashed",size=1/2) + 
          scale_y_continuous(breaks=break_vals) 
        } else { 
          p <- p + geom_hline(yintercept=1,linetype="dashed",size=1/2)}; p
      # save
      # folders
      subfldr_name<-"median_interquant_by_param_value/peak_duration/"
      if (!dir.exists(here::here(foldername,subfldr_name))) {
        dir.create(here::here(foldername,subfldr_name))}
      # filename
      plot_filename <- here::here(foldername,paste0(subfldr_name,"summ_stats_",
                paste0(c(sel_vars[k_plot_var],sel_vartype),collapse="_"),"_",
                sel_par,ifelse(start_year==2020,"_incl2020",""),".png"))
      # save
      ggsave(plot_filename,width=28,height=16,units="cm")
      # print progress of loop
      print(paste0(c(sel_vars[k_plot_var],sel_pars[k_plot_par],sel_vartype),collapse=", "))
    }
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Evaluating DYNAMICS as a function of parameters
comp_year <- c(2018,2019)[2] # 
dyn_all_parsets_broad_age_params <- left_join(
  dyn_all_parsets_broad_age %>% 
    filter(par_id %in% parsets_regular_dyn$par_id & 
             date>as.Date(paste0(comp_year,"-06-01")) & date<as.Date("2024-05-01")) %>% ungroup() %>% 
    mutate(epi_year=ifelse(week(date)>=epi_year_week_start,year(date),year(date)-1)) %>% 
    rename(incid_case=value) %>% 
    pivot_longer(c(incid_case,incid_hosp)) %>% 
    rename(varname=name,varvalue=value) %>% rename(value=varvalue) %>%
    ungroup() %>% group_by(agegroup_broad,par_id,varname) %>%
    mutate(peak_value_sel_year=max(value[epi_year %in% comp_year])) %>% # peak level of selected year
    group_by(epi_year,par_id,agegroup_broad,varname) %>% mutate(epi_week=row_number()) %>%
    group_by(agegroup_broad,par_id,varname) %>% 
    mutate(peak_week_sel_year=min(epi_week[value==peak_value_sel_year # peak week of year of reference
                                           & (epi_year %in% comp_year)])) %>% 
    ungroup() %>% 
    mutate(peak_week_distance=epi_week-peak_week_sel_year) %>% # distance from peak week of sel year
    mutate(value_norm=value/peak_value_sel_year) %>% 
    select(!peak_value_sel_year),
  left_join(partable_regular_dyn %>% 
              select(par_id,seasforc_width_wks,R0,seasforce_peak,omega),
                      pred_pca %>% select(par_id,PC1),by="par_id"), by="par_id") %>% 
  mutate(age_exp_par_bins=findInterval(PC1,seq(-1,1,by=1/5))) %>% 
        group_by(age_exp_par_bins) %>% 
        mutate(age_exp_par_bins=round(mean(PC1),1)) %>% rename(value_abs=value) %>%
        select(!c(seasforc_width_wks,R0,seasforce_peak,PC1)) %>% 
  pivot_longer(c(omega,age_exp_par_bins)) %>%
  rename(parname=name,parvalue=value,value=value_abs) %>% 
  relocate(c(varname,value,value_norm),.after=parvalue)

#######
# calculate statistics at time points relative to peak week of selected reference year
# (this command can take 10-20 secs)
summ_dyn_all_parsets_broad_age_relat_time <- dyn_all_parsets_broad_age_params %>% 
  group_by(agegroup_broad,epi_year,peak_week_distance,parname,parvalue,varname) %>%
  summarise(mean=mean(value_norm),median=median(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2],
            ci95_low=quantile(value_norm,c(0.025,0.975))[1],
            ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>% 
  pivot_longer(c(mean,median,ci50_low,ci50_up,ci95_low,ci95_up)) %>% 
  rename(metric=name) %>% 
  group_by(agegroup_broad,parname,parvalue,varname,metric) %>%
  mutate(parname=ifelse(parname=="omega","waning","exposure (-1) <-> age (1)"),
         parvalue=ifelse(parname=="waning",1/parvalue,parvalue),value=round(value,3))
# dyn_all_parsets_broad_age_params is a large dataframe, free up memory by deleting it
rm(dyn_all_parsets_broad_age_params)

##############################################################
# Plot dynamics faceted by the age vs immunity dependence of susceptibility (Figure 4)

# FIGURE 4 and SI FIG : plot faceted by years, color-coded by parameter values
sel_years <- c("2019","2021","2022","2023")
for (k_par in c("exposure (-1) <-> age (1)","waning")){
  for (k_age in 2) {
    for (var_name in c("incid_case","incid_hosp")[2]) {
          sel_agegr<-c("<1y","1-2y","2-5y")[1:ifelse(k_age==1,2,3)]
  # subset data to plot
          df_plot <- summ_dyn_all_parsets_broad_age_relat_time %>% 
                filter(agegroup_broad %in% sel_agegr & epi_year %in% sel_years & 
                       (metric %in% c("median","ci50_low","ci50_up")) & 
                       (parname %in% k_par) & varname %in% var_name) %>% 
                pivot_wider(names_from=metric,values_from=value) %>% 
                rename(`epi-year`=epi_year)
  # color palette
  n_par<-length(unique(df_plot$parvalue))
  colorpal <- colorRampPalette(colors=c("blue","grey","red"))(n_par)
  p <- ggplot(df_plot,
    aes(x=peak_week_distance,group=parvalue,color=factor(parvalue))) + 
    geom_line(aes(y=median),size=1.1) +
    geom_ribbon(aes(ymin=ci50_low,ymax=ci50_up,fill=factor(parvalue)),color=NA,alpha=0.05) + 
    facet_grid(agegroup_broad~`epi-year`,
               labeller=labeller(`epi-year`=label_both),scales="free_y") +
    scale_x_continuous(limits=c(-23,15)) + theme_bw() + standard_theme + 
    theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
          axis.title.y=element_text(size=16),legend.position="top",
          legend.text=element_text(size=15),legend.title=element_text(size=16),
          strip.text=element_text(size=18)) +
    xlab("distance in weeks from (pre-pandemic) peak week") +  
    ylab(paste0("weekly hospitalisations (1=pre-pandemic peak)")) + labs(color=k_par,fill=k_par) +
    geom_hline(yintercept=1,linetype="dashed",size=1/3) + 
    geom_vline(xintercept=c(-9,9),linetype="dashed",size=1/3)
if (n_par>3) {p <- p + scale_color_manual(values=colorpal) + 
              scale_fill_manual(values=colorpal)}; p 
  # save
  # folder to save
  if (!dir.exists(here::here(foldername,"dynamics"))) {
    dir.create(here::here(foldername,"dynamics"))}
  # filename
  plot_fn <- here::here(foldername, paste0("dynamics/weekly_hosp_by_",
                               ifelse(grepl("exp",k_par),"age_exp",k_par),
                            "_norm_",comp_year,"_peak_until",max(sel_years),
                            ifelse(k_age>1,"_1_5y",""),".png")) 
  # save
  ggsave(plot_fn,width=30,height=18,units="cm"); print(gsub(foldername,"",plot_fn))
     }
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# quantify total variation in hospitalisations (figures quoted in Discussion)
# comp_year<-2019
hosp_tot_by_parsets <- results_summ_all_hosp %>% 
  filter(epi_year %in% c(comp_year,2021) & agegroup<=7) %>% 
  select(!c(hosp_seas,max_incid_week,attack_rate_perc,seas_share,inf_in_seas,
            hosp_seas,agegroup_name,final,omega,exp_dep,age_dep)) %>% 
  group_by(epi_year,par_id) %>% 
  summarise(hosp_tot=sum(hosp_tot),inf_tot=sum(inf_tot),PC1=unique(PC1)) %>%
  pivot_longer(c(hosp_tot,inf_tot)) %>% 
  arrange(par_id,name) %>% group_by(par_id,name) %>% 
  mutate(value_norm=value/value[epi_year==comp_year])

# statistics on relative change from 2019 to 2021, ALL parsets
results_summ_all_hosp %>% 
  filter(epi_year %in% c(comp_year,2021) & agegroup<=7) %>% 
  select(!c(hosp_seas,max_incid_week,attack_rate_perc,seas_share,inf_in_seas,
            hosp_seas,agegroup_name,final,omega,exp_dep,age_dep)) %>% 
  pivot_longer(c(hosp_tot,inf_tot)) %>% 
  arrange(par_id,name) %>% 
  group_by(agegroup_broad,par_id,name) %>% 
  mutate(value_norm=value/value[epi_year==comp_year]) %>% 
  filter(epi_year==2021) %>% ungroup() %>% 
  summarise(mean=mean(value_norm),median=median(value_norm),
            min=min(value_norm),max=max(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2])

# CUMUL burden increase
cumul_hosp_by_age <- results_summ_all_hosp %>% 
  filter(epi_year %in% c(comp_year,2021) & agegroup<=7) %>% 
  select(!c(hosp_seas,max_incid_week,attack_rate_perc,seas_share,
            inf_in_seas,hosp_seas,agegroup_name,final,omega,exp_dep,age_dep)) %>% 
  pivot_longer(c(hosp_tot,inf_tot)) %>% 
  arrange(par_id,name) %>% 
  group_by(agegroup_broad,par_id,name) %>% 
  mutate(value_norm=value/value[epi_year==comp_year]) %>% 
  mutate(age_exp_par_bins=findInterval(PC1,seq(-1,1,by=1/5))) %>% 
  group_by(age_exp_par_bins) %>% 
  mutate(age_exp_par_bins=round(mean(PC1),1)) %>% ungroup() 

# CUMUL burden SUM relative increase
cumul_hosp_by_age %>% filter(name %in% "hosp_tot") %>% 
  group_by(epi_year,par_id) %>% 
  summarise(value=sum(value)) %>% group_by(par_id) %>% 
  mutate(value_norm=value/value[epi_year==comp_year]) %>% 
  filter(epi_year==2021) %>% ungroup() %>%
  summarise(mean=mean(value_norm),median=median(value_norm),
            min=min(value_norm),max=max(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2])

# CUMUL burden increase each age group, all params
cumul_hosp_by_age %>% group_by(agegroup_broad) %>%
  filter(epi_year==2021) %>%
  summarise(mean=mean(value_norm),median=median(value_norm),
            min=min(value_norm),max=max(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2])


# CUMUL BURDEN range at maximal exposure dependence, by age groups
cumul_hosp_by_age %>% group_by(agegroup_broad) %>%
  filter(age_exp_par_bins==min(age_exp_par_bins) & epi_year==2021) %>%
          summarise(mean=mean(value_norm),median=median(value_norm),
              min=min(value_norm),max=max(value_norm),
              ci50_low=quantile(value_norm,c(0.25,0.75))[1],
              ci50_up=quantile(value_norm,c(0.25,0.75))[2])

# CUMUL BURDEN range at minimal exposure dependence
cumul_hosp_by_age %>% group_by(agegroup_broad) %>%
filter(age_exp_par_bins==max(age_exp_par_bins) & epi_year==2021) %>%
  summarise(mean=mean(value_norm),
            median=median(value_norm),
            min=min(value_norm),max=max(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2])

# complete range of increase in (weekly) peaks in ALL AGE GROUPS
parsets_max_incid_seas_length %>% 
  filter(varname %in% "incid_hosp" & vartype %in% "max_value" & 
           epi_year %in% c(comp_year,2021) & 
           agegroup_broad %in% c("<1y","1-2y","2-5y")) %>%
  group_by(epi_year,par_id) %>% 
  summarise(value=sum(value),kappa=unique(age_exp_par_bins)) %>%
  group_by(par_id) %>% 
  mutate(value_norm=value/value[epi_year==comp_year]) %>% 
  arrange(par_id) %>% ungroup() %>% 
  filter(epi_year==2021 & kappa==min(kappa)) %>% 
  summarise(mean=mean(value_norm),median=median(value_norm),
            min=min(value_norm),max=max(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2])

# complete range of increase in (weekly) peaks in AGE GROUPS >1y
parsets_max_incid_seas_length %>% 
  filter(varname %in% "incid_hosp" & vartype %in% "max_value" & 
           epi_year %in% c(comp_year,2021) & 
           agegroup_broad %in% c("1-2y","2-5y")) %>%
  group_by(epi_year,par_id) %>% 
  summarise(value=sum(value),kappa=unique(age_exp_par_bins)) %>%
  group_by(par_id) %>% 
  mutate(value_norm=value/value[epi_year==comp_year]) %>% 
  arrange(par_id) %>% ungroup() %>% 
  filter(epi_year==2021 & kappa==min(kappa)) %>% 
  summarise(mean=mean(value_norm),median=median(value_norm),
            min=min(value_norm),max=max(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2])


# ranges of increase in peak values by agegroup AT MAX EXP DEP
parsets_max_incid_seas_length %>% 
  filter(varname %in% "incid_hosp" & vartype %in% "max_value" & 
          epi_year %in% c(comp_year,2021) & 
          agegroup_broad %in% c("<1y","1-2y","2-5y")) %>%
  group_by(par_id,agegroup_broad) %>% 
  mutate(value_norm=value/value[epi_year==comp_year]) %>% 
  arrange(par_id) %>% ungroup() %>% 
  # strong exposure dependence
  filter(epi_year==2021 & age_exp_par_bins==min(age_exp_par_bins)) %>% 
  group_by(agegroup_broad) %>%
  summarise(mean=mean(value_norm),median=median(value_norm),
            min=min(value_norm),max=max(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2])

# ranges of increase in peak values by agegroup AT MAX AGE DEP
parsets_max_incid_seas_length %>% 
  filter(varname %in% "incid_hosp" & vartype %in% "max_value" & 
           epi_year %in% c(comp_year,2021) & 
           agegroup_broad %in% c("<1y","1-2y","2-5y")) %>%
  group_by(par_id,agegroup_broad) %>% 
  mutate(value_norm=value/value[epi_year==comp_year]) %>% 
  arrange(par_id) %>% ungroup() %>% 
  # strong exposure dependence
  filter(epi_year==2021 & age_exp_par_bins==max(age_exp_par_bins)) %>% 
  group_by(agegroup_broad) %>%
  summarise(mean=mean(value_norm),median=median(value_norm),
            min=min(value_norm),max=max(value_norm),
            ci50_low=quantile(value_norm,c(0.25,0.75))[1],
            ci50_up=quantile(value_norm,c(0.25,0.75))[2])

# % simuls where peak goes down if MAX AGE DEP
parsets_max_incid_seas_length %>% 
  filter(varname %in% "incid_hosp" & vartype %in% "max_value" & 
           epi_year %in% c(comp_year,2021) & 
           agegroup_broad %in% c("<1y","1-2y","2-5y")) %>%
  group_by(par_id,agegroup_broad) %>% 
  mutate(value_norm=value/value[epi_year==comp_year]) %>% 
  arrange(par_id) %>% ungroup() %>% 
  # strong exposure dependence
  filter(epi_year==2021 & age_exp_par_bins==max(age_exp_par_bins)) %>% 
  group_by(agegroup_broad) %>%
  summarise(perc_decrease=sum(value_norm<1)/n(),perc_increase=sum(value_norm<1)/n())

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Plot UK RSV case data, calculate seasonal concentration

season_weeks=c(13,40)
resp_detects_weekly_all_age <- read_csv(here(foldername,
                  "Respiratory_viral_detections_by_any_method_UK.csv")) %>% 
  mutate(year_week=factor(paste0(Year,"-",Week),unique(paste0(Year,"-",Week))), 
         RSV_rolling_av=rollmean(RSV,k=7,align="center",fill=NA) ) %>% 
  select(-(contains("virus")|contains("flu"))) %>% 
  mutate(epi_year=ifelse(Week>=season_weeks[2],Year+1,Year)-min(Year)+1) %>% 
  group_by(epi_year) %>% mutate(perc_yearly=RSV/sum(RSV)) %>% group_by(epi_year) %>% 
  mutate(season_share=sum(perc_yearly[Week>=season_weeks[2] | Week<=season_weeks[1]]),
         on_off_season=ifelse(findInterval(Week,season_weeks+c(1,0))==1,"off","on"))

# SI Table 5
resp_detects_weekly_all_age_means_shares <- left_join(
  resp_detects_weekly_all_age %>% 
          group_by(epi_year,on_off_season) %>% 
          summarise(mean_on_off=mean(RSV,na.rm=T),cal_year=paste0(unique(Year),collapse="_")) %>% 
            filter(epi_year>1&epi_year<8) %>% 
            mutate(cal_year=ifelse(on_off_season %in% "off",paste0(as.numeric(cal_year)-1,"_",
                                      as.numeric(cal_year)),cal_year)) %>% 
            pivot_wider(names_from=on_off_season,values_from=mean_on_off,names_prefix="mean_"),  
  resp_detects_weekly_all_age %>% 
            group_by(epi_year) %>% 
            summarise(season_share=unique(season_share)),by="epi_year" ) %>%
  mutate(on_off_ratio=round(mean_on/mean_off,1),
  mean_on=round(mean_on,1),mean_off=round(mean_off,1),season_share=round(season_share,3))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SI Figure 1 (RSV case data with age)
resp_virus_data_uk_tidy <- read_csv(here::here(foldername,
                            "Respiratory_viral_detections_by_any_method_UK_Ages.csv")) %>% 
  pivot_longer(!c("Year","startweek","Age")) %>% 
  mutate(Age=factor(gsub(" Y","Y",Age),levels=unique(gsub(" Y","Y",Age))),
         date=as.Date(paste(Year,startweek,1,sep="-"),"%Y-%U-%u"))
# PLOT all years
ggplot(subset(resp_virus_data_uk_tidy,name %in% "RSV"),
       aes(x=date,y=value,group=Age)) + 
  geom_area(aes(fill=Age),position=position_stack(reverse=T),color="black",size=1/4) +
  scale_x_date(breaks="2 month",expand=expansion(0.01,0)) + 
  scale_y_continuous(expand=expansion(0.01,0)) + 
  xlab("")  +  ylab("number of reported RSV cases") + labs(fill="") + 
  theme_bw() + standard_theme  +
  theme(axis.text.x=element_text(size=13),axis.text.y=element_text(size=14),
        legend.position="bottom",legend.text=element_text(size=15),
        axis.title.y=element_text(size=16))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### THE END :)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PS: 
# To run individual simulations and plot results go to "indiv_simul.R" 
#