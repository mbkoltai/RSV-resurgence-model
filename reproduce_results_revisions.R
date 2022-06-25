# This script is for reproducing the results and figures of the manuscript at [...]
# Mihaly Koltai, May/2022
####

# clear workspace
rm(list=ls())
# To set the path we need the "here" package
if (!any(row.names(installed.packages()) %in% "here")) {install.packages("here")}; library(here)
# load constant parameters and functions for simulations, specify folder where inputs are stored
source(here("load_params.R"))
foldername <- "repo_data/"; figs_folder <- "repo_data/FIGS"
# NPI dates
npi_dates<-as.Date(c("2020-03-26","2021-05-17"))
# set up the table of parameter vectors by Latin Hypercube Sampling (LHS)
new_partable=F
if (new_partable){
source("fcns/create_lhs_partable.R")
} else {
  partable <- read_csv("repo_data/partable_full_lhs.csv")
}

# check the size of objects (>x Mb) in the workspace by: fcn_objs_mem_use(min_size=1)

# start date of saved simulations
start_date_dyn_save <- "2016-09-01"

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

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
# contact_red<-0.9
# agegroup_res<-"broad_age" # this if summary stats should be saved in 0-1, 1-2, 2-5, 65+ groups only
# # these are parameters selected by criteria of 1) attack rates 2) seasonal concentration of cases
# partable_filename <- "repo_data/partable_full_lhs.csv"; n_row <- nrow(read_csv(partable_filename))
# # we split the parameter table into `n_core` batches and run them in parallel, the sh file will launch the jobs
# # write the files launching jobs
# command_print_runs <- paste0(c("Rscript fcns/write_run_file.R",n_core,n_row,
#                                  simul_length_yr,n_post_npi_yr,contact_red,
#                                  partable_filename,"SAVE sep_qsub_files",start_date_dyn_save,
#                                  agegroup_res,memory_max),collapse=" ")
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
# this is a ZIP file in the folder, first unzip: 
# unzip("repo_data/results_summ_all.zip",exdir = "repo_data/")
results_summ_all <- read_csv("repo_data/results_summ_all.csv")

# param sets filtered out bc of attack rates or seasonal concentr
AR_crit=8; seas_conc_crit=AR_crit; seas_share_cutoff=85/100
# score: how many age groups satisfy
fullscan_score_AR_seasconc <- left_join(
            results_summ_all %>%
                        filter(epi_year %in% c(2017,2018)) %>% 
                        group_by(par_id,agegroup) %>%
                        summarise(attack_rate_perc=mean(attack_rate_perc),
                                  seas_share=mean(seas_share)),
            estim_attack_rates %>% 
                select(!c(attack_rate,sympt_attack_rate,n_test,RSV_posit,RSV_sympt_posit))  %>%
                mutate(agegroup=as.numeric(factor(agegroup_name,
                                levels=agegroup_name))),by="agegroup") %>%
                mutate(attack_rate_check=(attack_rate_perc>=min_est & attack_rate_perc<=max_est),
                       seas_share_check=(seas_share>=seas_share_cutoff)) %>%
                group_by(par_id) %>% summarise(n_attack_rate_check=sum(attack_rate_check), 
                                                 n_seas_share_check=sum(seas_share_check))
# summary stats of filtering (how many accepted/rejected and by what criteria)
filtering_summ_stats <- fullscan_score_AR_seasconc %>% 
                    mutate(attack_rate_fail=n_attack_rate_check<AR_crit,
                           seas_conc_fail=n_seas_share_check<seas_conc_crit,
                           both_fail=attack_rate_fail&seas_conc_fail) %>% 
                    summarise(n_both_fail=sum(both_fail),
                              n_attack_rate_fail=sum(attack_rate_fail)-n_both_fail,
                              n_seas_conc_fail=sum(seas_conc_fail)-n_both_fail,
                              n_accepted=sum(!attack_rate_fail&!seas_conc_fail))
# filter for parameter sets where x/11 age groups satisfy criteria 
# for attack rates and seasonal concentration
sel_yrs<-2019; n_sel_yr=length(sel_yrs)
# parameter sets where x/11 age-groups satisfy the AR and seasonal concentr criteria
parsets_AR_seas_share <- (fullscan_score_AR_seasconc %>% 
                            filter(n_attack_rate_check>=AR_crit & 
                                     n_seas_share_check>=seas_conc_crit))$par_id
all_sum_inf_epiyear_age_filtered <- results_summ_all %>% filter(par_id %in% parsets_AR_seas_share)
# this leads to ~4.5e3 out of 2e4 parameter sets:
length(unique(all_sum_inf_epiyear_age_filtered$par_id))
# select parameter sets matching the first two criteria
create_AR_seas_filtered_partable=T
if (create_AR_seas_filtered_partable){
partable_filtered_AR_seasconc <- partable %>% filter(par_id %in% parsets_AR_seas_share)
write_csv(partable_filtered_AR_seasconc,here(foldername,"partable_filtered_AR_seasconc.csv")) 
} else {
partable_filtered_AR_seasconc <- read_csv(here("repo_data/partable_filtered_AR_seasconc.csv")) 
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SI Figure 14: plot total infections in epi-years 2017,17,18 compared to 2019
# this is to see if seasonality is annual or there are large differences between years for some parameter sets

# Plot cumul incid (relative to 2019) of accepted parsets
results_summ_all %>% 
  filter(epi_year<2020 & 
           par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id)) %>%
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],
                              levels=rsv_age_groups$agegroup_name)) %>%
  group_by(par_id,agegroup_name) %>% 
  mutate(inf_in_seas=inf_in_seas/(inf_in_seas[epi_year==2019])) %>% # *1+(1/26)
  filter(epi_year<2019) %>%
ggplot() + 
  geom_jitter(aes(x=factor(epi_year),y=inf_in_seas,group=par_id),alpha=1/8,size=1/2) +
  facet_wrap(~agegroup_name,scales = "free") + scale_y_log10() +
  xlab("epi-year") + ylab("number of infections compared to 2019 (1=2019 level)") + 
  theme(legend.position="top") + theme_bw() + standard_theme
# save
ggsave(here(figs_folder,"2020_RSV_check_log.png"),width=32,height=20,units="cm")

# Biennial patterns? Number of param sets where cumul incidence is within 15% of 2019 level (or not)
irreg_toler=0.15; scaling_2019=1+1/26
results_summ_all %>% 
  filter(epi_year<2020 & 
           par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id)) %>%
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],
                              levels=rsv_age_groups$agegroup_name)) %>%
  group_by(par_id,agegroup_name) %>% 
  mutate(inf_in_seas=inf_in_seas/(scaling_2019*inf_in_seas[epi_year==2019])) %>%
  filter(epi_year<2019) %>% group_by(par_id,epi_year) %>% 
  summarise(regular_by_agegr=sum(inf_in_seas>1-irreg_toler & 
                                   inf_in_seas<1+irreg_toler)) %>%
  group_by(epi_year) %>% 
  summarise(regular=sum(regular_by_agegr>9),irregular=n()-regular)

# parsets with regular (annual) vs irregular patterns
reg_irreg_parsets <- results_summ_all %>% 
  filter(epi_year<2020 & par_id %in% parsets_AR_seas_share) %>%
  group_by(par_id,agegroup) %>% 
  mutate(inf_in_seas=inf_in_seas/(scaling_2019*inf_in_seas[epi_year==2019])) %>%
  filter(epi_year<2019) %>% group_by(par_id,epi_year) %>% 
  summarise(regular_by_agegr=sum(inf_in_seas>1-irreg_toler & 
                                   inf_in_seas<1+irreg_toler)) %>% 
  group_by(par_id) %>% summarise(regular=all(regular_by_agegr==11))
# number of regular/irregular 
reg_irreg_parsets %>% summarise(regular=sum(regular),irreg=n()-regular)

# CDF of percentage deviation from 2019
results_summ_all %>% 
  filter(epi_year<2020 & 
           par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id)) %>%
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],
                              levels=rsv_age_groups$agegroup_name)) %>%
  group_by(par_id,agegroup_name) %>% 
  mutate(inf_in_seas=inf_in_seas/(scaling_2019*inf_in_seas[epi_year==2019])) %>%
  filter(epi_year<2019) %>% 
  select(par_id,agegroup,inf_in_seas) %>% 
  mutate(inf_in_seas=abs(1-inf_in_seas)) %>%
ggplot(aes(inf_in_seas)) + stat_ecdf(geom="step") +
  facet_wrap(~agegroup,nrow=2,labeller=labeller(agegroup=label_both)) + 
  geom_vline(xintercept=0.15,linetype="dashed",size=1/2,color="red") +
  scale_x_log10(limits=c(0.001,1)) + # ,expand=expansion(0.02,0)
  xlab("relative difference in incidence") + ylab("CDF") + labs(color="agegroups")+
  theme_bw() + standard_theme + 
  theme(legend.position="top",strip.text=element_text(size=14),
    legend.title=element_text(size=15),legend.text=element_text(size=14),
    axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
    axis.title.x=element_text(size=16))
# save
ggsave(here(figs_folder,"interyear_difference_cumul_incid_reg_dyn.png"),
       width=35,height=25,units="cm")

# they haven't settled or stable biennial? read in one of the dynamic batches to check
# unzip first: 
# unzip("repo_data/dyn_parsets_main1_324.zip",exdir="repo_data/")
dyn_parsets_main1_324 <- read_csv("repo_data/dyn_parsets_main1_324.csv") %>%
  filter(par_id %in% unique(all_sum_inf_epiyear_age_filtered$par_id)) %>% 
  mutate(date=as.Date("2016-09-01")+t-min(t)) %>% filter(date<as.Date("2020-04-01"))
# plot
ggplot(dyn_parsets_main1_324 %>% filter(date<as.Date("2020-04-01") & 
                  par_id %in% reg_irreg_parsets$par_id & agegroup==1)) +
  geom_line(aes(x=date,y=value),alpha=1/2) + facet_wrap(~par_id,scales="free_y") + 
  xlab("") + ylab("infections") + theme_bw() + standard_theme
# curves look stabilised, but let's also calculate

# % difference in dynamic curves, calculated as: sum(abs(diff(2016,2018)))/sum(mean(2016,2018)) 
# and same for 2017/2019
interyr_diff_cumul_incid <- left_join(
  dyn_parsets_main1_324 %>% filter(date<as.Date("2020-04-01")),
                                      reg_irreg_parsets) %>%
      mutate(week=isoweek(date)) %>% filter(week>=35|week<=9) %>% 
      mutate(epi_year=ifelse((week>=35 & yday(date)>245) | 
                               yday(date)>245,year(date),year(date)-1)) %>% 
      filter(epi_year>=2016) %>%
      group_by(agegroup,epi_year,par_id,regular) %>% 
      arrange(epi_year,date) %>% 
      mutate(epi_day=row_number()) %>% filter(epi_day<=180) %>% 
      group_by(par_id,epi_day,agegroup,regular) %>%
      summarise(`diff_2016_18`=abs(value[epi_year==2016]-value[epi_year==2018]),
                `diff_2017_19`=abs(value[epi_year==2017]-value[epi_year==2019]),
      mean_2016_18=mean(value[epi_year %in% c(2016,2018)]),
      mean_2017_19=mean(value[epi_year %in% c(2017,2019)])) %>%
      group_by(agegroup,par_id,regular) %>% 
      summarise(rel_diff_16_18=sum(diff_2016_18)/sum(mean_2016_18),
                rel_diff_17_19=sum(diff_2017_19)/sum(mean_2017_19)) %>%
      group_by(par_id,regular) %>% 
      summarise(rel_diff_16_18=mean(rel_diff_16_18)*100,
                rel_diff_17_19=mean(rel_diff_17_19)*100) %>%
      group_by(regular) %>% 
      summarise(perc_diff_16_18=mean(rel_diff_16_18),
                perc_diff_17_19=mean(rel_diff_17_19))
# all under 1%, on average under 0.8% difference, so simulations have stabilised

# Plot curves comparing 2016 to 2018 and 2017 to 2019
left_join(dyn_parsets_main1_324 %>% filter(date<as.Date("2020-04-01")),
          reg_irreg_parsets) %>%
  mutate(week=isoweek(date)) %>% filter(week>=35|week<=9) %>%
  mutate(epi_year=ifelse((week>=35 & yday(date)>245) | 
                           yday(date)>245,year(date),year(date)-1)) %>%
  filter(epi_year>=2016) %>%
  group_by(agegroup,epi_year,par_id,regular) %>% arrange(epi_year,date) %>% 
  mutate(epi_day=row_number()) %>% 
  filter(epi_day<=180 & agegroup==2 & !regular) %>%
  mutate(year_type=epi_year %% 2) %>% group_by(year_type) %>% 
  mutate(year_rank=ifelse(epi_year==min(epi_year),0,1)) %>%
ggplot(aes(x=epi_day,y=value,color=factor(epi_year),
           group=epi_year,linetype=factor(year_rank))) + 
  geom_line() + facet_wrap(~par_id,scales="free_y") + 
  xlab("day of season") + labs(linetype="",color="epi-year") +
  theme_bw() + standard_theme
# save
# ggsave(here("simul_output/interyear_difference_irregular_parsets.png"),width=30,height=25,units="cm")
ggsave(here(figs_folder,"interyear_difference_regular_parsets.png"),width=30,height=25,units="cm")

# Plot again the comparison to 2019 for regular parameter sets only
results_summ_all %>% 
  filter(epi_year<2020 & par_id %in% parsets_AR_seas_share) %>%
  filter(par_id %in% reg_irreg_parsets$par_id[reg_irreg_parsets$regular]) %>%
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],
                              levels=rsv_age_groups$agegroup_name)) %>%
  group_by(par_id,agegroup_name) %>% 
  mutate(inf_in_seas=inf_in_seas/(scaling_2019*inf_in_seas[epi_year==2019])) %>%
  filter(epi_year<2019) %>%
ggplot() + geom_jitter(aes(x=factor(epi_year),y=inf_in_seas,group=par_id),alpha=1/8) +
  facet_wrap(~agegroup_name) + xlab("epi-year") + ylab("number of infections (1=2019 level)") + 
  theme(legend.position="top") + theme_bw() + standard_theme
# save
ggsave(here(figs_folder,"/cumul_incid_compared2019_by_agegrs.png"),width=30,height=25,units="cm")

# keep only parameters with regular annual patterns
regular_dyn_parset_load <- T
if (!regular_dyn_parset_load){
partable_regular_dyn <- partable_filtered_AR_seasconc %>% 
  filter(par_id %in% reg_irreg_parsets$par_id[reg_irreg_parsets$regular])
write_csv(partable_regular_dyn,"simul_output/2e4_parsets/crit_8_11/partable_regular_dyn.csv")
# keep outputs with correct attack rate
results_summ_all_reg_dyn <- results_summ_all %>% filter(par_id %in% partable_regular_dyn$par_id)
write_csv(results_summ_all_reg_dyn,"repo_data/results_summ_all_reg_dyn_FULL.csv")
} else {
  # unzip("results_summ_all_reg_dyn.zip",exdir="repo_data/")
  results_summ_all_reg_dyn <- read_csv(here(foldername,"results_summ_all_reg_dyn.csv"))
  partable_regular_dyn <- read_csv(here(foldername,"partable_regular_dyn.csv"))
}

# seasonal concentration density plot
results_summ_all_reg_dyn %>% filter(epi_year<2019) %>%
ggplot() + stat_ecdf(aes(seas_share,group=epi_year,color=factor(epi_year)),geom="step") +
  # geom_density(aes(x=seas_share,group=epi_year,color=factor(epi_year))) + 
  geom_vline(xintercept=seas_share_cutoff,color="red",linetype="dashed") +
  xlab("seasonal (w40≤x≤w13) share of cumulative incidence") + ylab("CDF") + 
  labs(color="epi-year") + theme_bw() + standard_theme
# save
ggsave(here(figs_folder,"seas_conc_CDF.png"),width=30,height=25,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# After filtering by attack rate and seasonal concentration we need full dynamics of *selected* parsets
# for likelihood calculations

# accepted parsets
create_dynamics_df <- F
if (create_dynamics_df) {
  source("fcns/create_dynamics_df.R")
} else {
  # download these two files from 
  # https://drive.google.com/file/d/1PAO0OjaI-MFugKwZVrD7dLFTktszi-1p/view?usp=sharing
  # and unzip
  all_dynamics_accepted <- read_csv(here(foldername,"all_dynamics_accepted_2e3.csv")) %>% 
                                mutate(date=as.Date(start_date_dyn_save)+t-min(t)) %>%
                                filter(date>=as.Date("2017-09-25") & date<as.Date("2024-10-01"))
  all_dynamics_rejected <- read_csv(here(foldername,"all_dynamics_rejected_2e3.csv")) %>% 
                                mutate(date=as.Date(start_date_dyn_save)+t-min(t)) %>%
                                filter(date>=as.Date("2017-09-25") & date<as.Date("2024-10-01"))
  # remove timepoints after SARI-Watch data
  # all_dynamics_accepted <- all_dynamics_accepted %>% filter(date<=as.Date("2020-05-11"))
  # all_dynamics_rejected <- all_dynamics_rejected %>% filter(date<=as.Date("2020-05-11"))
}


# load RSV hospitalisation counts
source("fcns/load_SARI_data.R")

# dataframes created: 
# SARIwatch_RSVhosp_under5_2018_2020_weekly_counts, SARIwatch_RSVhosp_over65_2018_2020_weekly_counts
# under-reporting: under_report_factor_under5, under_report_factor_over65y

# calculate weekly hospitalisations of accepted
# agegroup_mapping=findInterval(1:11,c(3,5,8,11))+1; names(agegroup_mapping)=paste0("agegroup",1:11)
# calculate hospitalisations
create_weekly_hosp_df = F
if (create_weekly_hosp_df){

# create dataframe of weekly hospitalisations, by broad age groups (0-1, 1-2, 2-5, 65+) and by <5 and 65+:
# simul_hosp_rate_weekly_all_broad_agegroups | simul_hosp_rate_weekly_under5_over65
# this takes ~2 mins for 4e3 param sets
source("fcns/create_weekly_hosp_df.R")

# calculate the (Poisson) likelihood per data point
  simul_hosp_rate_weekly_SARIdates_LLH <- left_join(
          left_join(SARIwatch_RSVhosp_under5_2018_2020_weekly_counts,
                    SARIwatch_RSVhosp_over65_2018_2020_weekly_counts,
                    by=c("wk_n","date","year")), 
          simul_hosp_rate_weekly_under5_over65, by="date") %>%
                mutate(simul_hosp_scaled=ifelse(grepl("65y",broad_age),
                                  simul_hosp_sum_full*under_report_factor_over65y*(
                                      pop_AGE65PLUS/rsv_age_groups$stationary_popul[11]),
                                  simul_hosp_sum_full*under_report_factor_under5*(
                                      pop_AGEUNDER5/sum(rsv_age_groups$stationary_popul[1:7])) ),
          log_lklh_poiss=dpois(x=ifelse(grepl("65y",broad_age),cases65plustotal,casesunder5total),
                              lambda=simul_hosp_scaled,log=T)) %>% 
    relocate(date,.after=wk_n) %>%
    group_by(par_id,broad_age) %>% 
    mutate(sum_neg_llh=-sum(log_lklh_poiss,na.rm=T),
           accepted=ifelse(par_id %in% partable_regular_dyn$par_id,"accepted","rejected"))
# save
write_csv(simul_hosp_rate_weekly_SARIdates_LLH,
          "repo_data/simul_hosp_rate_weekly_SARIdates_LLH.csv")
} else {
  # unzip: unzip("repo_data/simul_hosp_rate_weekly_LLH.zip",exdir="repo_data/")
  simul_hosp_rate_weekly_SARIdates_LLH <- 
    read_csv("repo_data/simul_hosp_rate_weekly_SARIdates_LLH.csv") %>%
        mutate(year_week=factor(year_week,unique(year_week)))
}

# `simul_hosp_scaled` means that simulated hospitalisations were scaled by the under-reporting factor

# calculate likelihoods by parameter set (accepted/rejected)
hosp_dyn_likelihoods <- simul_hosp_rate_weekly_SARIdates_LLH %>% 
  group_by(par_id,accepted,broad_age) %>% 
  summarise(sum_neg_llh=unique(sum_neg_llh)) %>%
  group_by(par_id,accepted) %>% mutate(sum_neg_llh_all=sum(sum_neg_llh)) %>% 
  pivot_longer(!c(par_id,accepted,broad_age)) %>%
  mutate(broad_age=ifelse(grepl("all",name),"all",broad_age )) %>% 
  select(!name) %>% distinct()

# plot likelihoods: <5y, 65+y, both combined
ggplot(hosp_dyn_likelihoods %>% mutate(broad_age=paste0(broad_age," hospitalisations")),
       aes(x=accepted,y=value,color=accepted)) + 
  geom_point(position=position_jitterdodge(seed=1,dodge.width=0.9),
             alpha=1/2,size=1.5,shape=21) +
  geom_boxplot(fill=NA,size=1/2,width=1/2,outlier.colour=NA,color="black") + 
  facet_wrap(~broad_age,scales = "free_y") +
  scale_color_manual(values=c("grey","blue")) + scale_y_log10() +
  xlab("") + ylab("negative log-likelihood") + labs(color="accepted") +
  theme_bw() + standard_theme # theme(legend.position="none")
# save
ggsave(here(figs_folder,"hosp_dyn_LLH_accepted_rejected_boxplot.png"),
       width=30,height=25,units="cm")

# calculate likelihoods for attack rates
create_rejected_summary_stats=F
if (create_rejected_summary_stats){
  # read_csv("simul_output/2e4_parsets/results_summ_all.csv")
  results_summ_rejected <- results_summ_all %>% 
    filter(!(par_id %in% unique(partable_regular_dyn$par_id)) )
  # save
  write_csv(results_summ_rejected,"simul_output/2e4_parsets/crit_8_11/results_summ_rejected.csv")
} else {
  # unzip first
    results_summ_rejected <- read_csv("repo_data/results_summ_rejected.csv")
}

create_LLH=T
# create LLH tables
if (create_LLH) {
  
likelihoods_attackrates <- bind_rows(
  results_summ_all_reg_dyn %>% 
    filter(par_id %in% unique(hosp_dyn_likelihoods$par_id)) %>% mutate(accepted=TRUE),
  results_summ_rejected %>% mutate(accepted=FALSE)) %>%
  filter(epi_year<2019) %>% 
  select(c(epi_year,par_id,agegroup,inf_in_seas_AR,attack_rate_perc,accepted)) %>%
  mutate(AR_log_binom_LLH=dbinom(x=estim_attack_rates$RSV_sympt_posit[agegroup],
         size=estim_attack_rates$n_test[agegroup],
         prob=attack_rate_perc/100,log=T)) %>%
  group_by(par_id,accepted) %>% 
  summarise(AR_negLLH_binom=-sum(AR_log_binom_LLH)) 

# sum of all likelihoods
all_likelihoods = left_join(
  hosp_dyn_likelihoods %>% 
    filter(!broad_age %in% "all") %>% mutate(accepted=accepted %in% "accepted") %>% 
    pivot_wider(names_from="broad_age",
                values_from="value",names_prefix="LLH from hospit. "),
  likelihoods_attackrates) %>% 
  mutate(`complete likelihood`=`LLH from hospit. <5y` + 
           `LLH from hospit. >65y`+AR_negLLH_binom) %>% 
  rename(`LLH from attack rates`=AR_negLLH_binom) %>% 
  pivot_longer(!c(par_id,accepted))
# SAVE
write_csv(all_likelihoods,"repo_data/all_likelihoods.csv") 
} else {
  all_likelihoods <- read_csv("repo_data/all_likelihoods.csv")
}

# plot as jitter plot
ggplot(all_likelihoods,aes(x=accepted,y=value,color=accepted)) + 
  geom_jitter(width=0.4,alpha=1/3) + #  geom_violin(fill=NA,show.legend=F) + 
  geom_boxplot(fill=NA,width=0.88,size=3/4,outlier.colour=NA,color="black") + # 
  facet_wrap(~name,scales="free_y",nrow=1) + scale_y_log10() +
  scale_color_manual(values=c("darkgrey","blue")) + # labs(color="accepted") + # coord_fixed(ratio=3) + 
  xlab("") + ylab("negative log-likelihood") + theme_bw() + standard_theme + # element_text(size=14)
  theme(strip.text=element_text(size=15),
        axis.text.x=element_blank(),axis.text.y=element_text(size=15),
        axis.title.y=element_text(size=18),
        legend.position="top",legend.title=element_text(size=15),
        legend.text=element_text(size=14)) 
  # manuscript_large_font_theme
# save
ggsave(here(figs_folder,"all_LLH_accepted_rejected_boxplot.png"),
       width=35,height=22,units="cm")

# density plot
median_llh = all_likelihoods %>% group_by(accepted,name) %>% 
                summarise(median_llh=median(value,na.rm=T))
ggplot(all_likelihoods,aes(x=value,color=accepted)) + 
  geom_density(aes(y=..scaled..)) + # aes(y=..scaled..) # aes(y=..count..)
  geom_vline(data=median_llh,aes(xintercept=median_llh,color=accepted),
             linetype="dashed",size=1/2,show.legend=F) +
  geom_text(data=median_llh,show.legend=F,
            aes(x=median_llh*0.88,y=1/8,color=accepted,label=round(median_llh))) +
  facet_wrap(~name,scales="free",nrow=3) + scale_x_log10(expand=expansion(0.01,0)) +
  scale_color_manual(values=c("black","blue")) + labs(color="accepted") + 
  xlab("negative log-likelihood") + ylab("density") + 
  theme_bw() + standard_theme + 
  theme(strip.text=element_text(size=15),axis.text.y=element_text(size=12))
# save
ggsave(here(figs_folder,"all_LLH_accepted_rejected_density.png"),width=35,height=22,units="cm")

# cumulative density of parsets up to LLH<x
ggplot(all_likelihoods) + stat_ecdf(aes(x=value,color=accepted),geom="step") + 
  facet_wrap(~name,scales="free",nrow=3) + 
  geom_segment(data=median_llh,
               aes(x=median_llh,xend=median_llh,y=0,yend=1/2,color=accepted),
               linetype="dashed",size=1/2,show.legend=F) +
  geom_text(data=median_llh,show.legend=F,
            aes(x=median_llh*1.12,y=0.44,color=accepted,label=round(median_llh))) +
  scale_x_log10(expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.01,0)) +
  scale_color_manual(values=c("black","blue")) + labs(color="accepted") + 
  xlab("negative log-likelihood") + ylab("cumulative density function") +
  theme_bw() + standard_theme + 
  theme(legend.position="top",strip.text=element_text(size=14),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16))
# SAVE
ggsave(here(figs_folder,"all_LLH_accepted_rejected_CDF.png"),width=35,height=22,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  
# what happens for parsets where fit is good?
# top 1%: array(quantile(all_likelihoods$value[grepl("complete",all_likelihoods$name)],probs=0.01,na.rm=T))

bind_rows(results_summ_rejected,
          results_summ_all_reg_dyn) %>% 
  filter(par_id %in% 
           (all_likelihoods %>% filter((name %in% "complete likelihood") & value<1.5e3))$par_id & 
        epi_year<2019) %>% 
  select(c(par_id,agegroup,attack_rate_perc,seas_share)) %>% 
  mutate(accepted=
           ifelse(par_id %in% unique(all_likelihoods$par_id[all_likelihoods$accepted]),
                  T,F)) %>% 
  pivot_longer(!c(par_id,agegroup,accepted)) %>% 
  mutate(name=ifelse(name %in% "attack_rate_perc",
                     "attack rate (%)", 
                     "seasonal concentration (%)")) %>%
ggplot() + 
  geom_jitter(aes(x=factor(agegroup),y=ifelse(grepl("seas",name),value*100,value),
                  color=accepted,group=accepted),alpha=1/4,
              position=position_jitterdodge(dodge.width=0.9,jitter.width=0.35)) + 
  geom_vline(xintercept=(1:10)+1/2,size=1/2) +
  geom_segment(data=estim_attack_rates %>% 
                 mutate(agegroup=row_number(),sympt_attack_rate=100*sympt_attack_rate) %>% 
                 select(c(agegroup,min_est,max_est,sympt_attack_rate)) %>% 
                 pivot_longer(!agegroup) %>% 
                 mutate(type=name,name="attack rate (%)"),
              aes(x=agegroup-1/2,xend=agegroup+1/2,y=value,yend=value,group=name,
                  linetype=ifelse(type %in% "sympt_attack_rate","solid","dashed")),
                  size=1/3,show.legend=F) +
  geom_hline(aes(yintercept=ifelse(name %in% "seas_share",0.85,NA)),linetype="dashed") +
  facet_wrap(~name,scales="free",nrow=2) + xlab("") + ylab("") +
  scale_x_discrete(expand=expansion(0,0)) + scale_y_log10() +  
  scale_color_manual(values = c("darkgrey","blue")) + theme_bw() + standard_theme + 
  theme(strip.text=element_text(size=15),
        axis.text.x=element_blank(),axis.text.y=element_text(size=15),
        axis.title.y=element_text(size=18),
        legend.position="top",legend.title=element_text(size=15),legend.text=element_text(size=14)) 
# it's because attack rates are out of the range for agegroups > 8
# ggsave(here("simul_output/2e4_parsets/goodfits_negLLH_below1e3_filtered.png"),width=35,height=16,units="cm")
ggsave(here(figs_folder,"goodfits_negLLH_below1500_filtered.png"),width=35,height=20,units="cm")

# dynamics
goodfit_pars <- (all_likelihoods %>% filter((name %in% "complete likelihood") & value<1.5e3))$par_id
sari_hosp_data_joint = simul_hosp_rate_weekly_SARIdates_LLH %>% ungroup() %>% 
  select(c(year_week,casesunder5total,cases65plustotal,broad_age,year)) %>% 
  distinct() %>% pivot_longer(!c(year_week,broad_age,year)) %>%
  filter(!(broad_age %in% "<5y" & name %in% "cases65plustotal")) %>%
  filter(!(broad_age %in% ">65y" & name %in% "casesunder5total"))
# plot
simul_hosp_rate_weekly_SARIdates_LLH %>% 
  select(c(par_id,accepted,broad_age,simul_hosp_scaled,year_week,year)) %>%
  filter(par_id %in% goodfit_pars) %>% 
  mutate(year_week=factor(year_week,levels=unique(year_week))) %>%
ggplot(aes(x=year_week)) + 
  geom_line(aes(y=simul_hosp_scaled,color=accepted,group=par_id),alpha=1/2) + 
  geom_point(data=sari_hosp_data_joint,aes(y=value)) +
  geom_line(data=sari_hosp_data_joint,aes(y=value,group=1),linetype="dashed",size=1/2) +
  scale_color_manual(values=c("grey","blue")) + 
  scale_x_discrete(breaks=show_every_nth(n=2)) + 
  facet_grid(broad_age~year,scales="free") + 
  xlab("") + ylab("hospitalisations (count)") + 
  theme_bw() + standard_theme
# save
# ggsave(here("simul_output/2e4_parsets/goodfits_negLLH_below2e3_filtered.png"),width=35,height=16,units="cm")
# ggsave(here("simul_output/2e4_parsets/goodfits_negLLH_below1.5e3_filtered.png"),width=35,height=16,units="cm")
ggsave(here(figs_folder,"goodfits_negLLH_below1e3_filtered.png"),width=35,height=16,units="cm")

### ### ### ### ### ### ### ### ### ### ###
# Plot showing why 'good' (low negLLH) parsets are filtered out
x_dodge=1/5
all_crit_summarised = 
left_join(
  left_join(
    left_join(  # %>% filter(par_id %in% unique(all_dynamics_accepted$par_id))
            bind_rows(results_summ_all_reg_dyn %>% mutate(accepted=T), 
                    results_summ_rejected  %>% mutate(accepted=F)) %>% 
                    filter(epi_year<2020) %>% 
                    select(c(par_id,agegroup,attack_rate_perc,seas_share,accepted)), 
            reg_irreg_parsets),
    fullscan_score_AR_seasconc) %>% group_by(par_id) %>% 
    mutate(seas_share=mean(seas_share)) %>%
    select(!c(agegroup,attack_rate_perc)) %>% distinct(),
  all_likelihoods %>% 
    filter(name %in% "complete likelihood") %>% select(!name) %>% rename(negLLH=value)) %>%
    mutate(x_pos_AR=ifelse(accepted,
                         n_attack_rate_check+x_dodge+runif(1,min=-1/7,max=1/7),
                         n_attack_rate_check-x_dodge+runif(1,min=-1/7,max=1/7)))
# n_attack_rate_check=paste0(n_attack_rate_check,"/11"),
# n_attack_rate_check=factor(n_attack_rate_check,levels=unique(n_attack_rate_check))

# if sampling 1000-1000 of all parsets
# sample_pars <- c(sample(all_crit_summarised$par_id[all_crit_summarised$accepted],size=1000),
#                  sample(all_crit_summarised$par_id[!all_crit_summarised$accepted],size=1000))

# plot  
ggplot(all_crit_summarised) + #  %>% filter(par_id %in% sample_pars)
  geom_point(aes(x=x_pos_AR,y=negLLH,fill=ifelse(accepted,"accepted","rejected")),
             alpha=0.4,shape=21,stroke=0,size=3) +
  # put a cross (x) if attack rates <80% correct
  geom_point(data=all_crit_summarised %>% filter(n_seas_share_check<8),color="black",
            aes(x=x_pos_AR,y=negLLH,shape="seasonal concentration\n <8/11 correct")) +
  # red circle if seasons irregular
  geom_point(data=all_crit_summarised %>% filter(!regular), 
             aes(x=x_pos_AR,y=negLLH,color="biennial/irregular seasons"),
             shape=21,fill=NA,stroke=2/3,size=2.5) +
  scale_fill_manual(values=c("accepted"="blue","rejected"="darkgrey"),
                    guide=guide_legend(nrow=2,byrow=TRUE,
                    override.aes=list(color=NA,stroke=NA))) +
  scale_color_manual(values=c("biennial/irregular seasons"="red"),
                     guide=guide_legend(override.aes=list(size=4,stroke=1.2))) +
  scale_shape_manual(values=c("seasonal concentration\n <8/11 correct"=4),
                     guide=guide_legend(override.aes=list(size=3))) +
  geom_vline(xintercept=c(5:6,8:10)+1/2,size=1/2,linetype="dashed") + 
  geom_vline(xintercept=7.5,size=1/2) +  
  scale_x_continuous(breaks=1:11,expand=expansion(0.01,0),limits=c(4.6,11.36)) + 
  scale_y_log10(limits=c(720,1e4)) +
  xlab("age groups with correct attack rates") + ylab("Negative log-likelihood") +
  labs(fill="",color="",shape="",size="correct \nattack rates")+
  theme_bw() + standard_theme + 
  theme(strip.text=element_text(size=15),
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        legend.position="top",legend.title=element_text(size=15),
        legend.text=element_text(size=15)) 
# save
ggsave(here(figs_folder,"parsets_phase_diagram_x_attackrate.png"),width=35,height=16,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# SI FIGURE 4: simulated attack rates vs LIT estimates
attack_rate_limits <- estim_attack_rates %>% 
  mutate(agegroup=row_number(),sympt_attack_rate=100*sympt_attack_rate) %>% 
  select(c(agegroup,min_est,max_est,sympt_attack_rate)) %>% 
  pivot_longer(!agegroup) %>% 
  mutate(type=name,name="attack_rate_perc")
# plot
bind_rows(results_summ_all_reg_dyn %>% 
            filter(par_id %in% sample(unique(par_id),size=5e2)) %>% mutate(accepted=T),
          results_summ_rejected %>% 
            filter(par_id %in% sample(unique(par_id),size=5e2)) %>% mutate(accepted=F)) %>% 
  filter(epi_year<2020) %>% 
  select(c(par_id,agegroup,attack_rate_perc,accepted)) %>%
ggplot() + 
  geom_jitter(aes(x=factor(rsv_age_groups$agegroup_name[agegroup],
                           levels=unique(rsv_age_groups$agegroup_name)),
                  y=attack_rate_perc,color=accepted,group=accepted),
              position=position_jitterdodge(dodge.width=0.9,jitter.width=0.4),
              alpha=1/4,size=2/3) + 
  geom_vline(xintercept=(0:11)+1/2,size=1/2,color="grey") +
  geom_segment(data=attack_rate_limits,
               aes(x=agegroup-1/2,xend=agegroup+1/2,y=value,yend=value,group=name,
               linetype=ifelse(type %in% "sympt_attack_rate","solid","dashed")),
               show.legend=F,size=1/2,color="red") +
  scale_y_log10(limits=c(0.1,110),expand=expansion(0.02,0)) + 
  scale_x_discrete(expand=expansion(0,0)) +
  scale_color_manual(values=c("black","blue"),guide=guide_legend(override.aes=list(size=3))) + 
  xlab("age group") + ylab("attack rate (%)") +
  theme_bw() + standard_theme + 
  theme(legend.position="top",axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13))
# SAVE
ggsave(here(figs_folder,"attack_rate_comparison_with_lit.png"),width=25,height=20,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# distribution of all parameters separated by accepted/rejected
# partable <- read_csv("repo_data/partable_full_lhs.csv")

# jitter with boxplot
plot_partable_histogram = partable %>% 
  mutate(accepted=par_id %in% partable_regular_dyn$par_id) %>% select(!const_delta) %>%
  mutate(`waning (days)`=1/omega,`peak R0`=R0*(1+seasforce_peak),
         `maximal forcing (% above baseline)`=1e2*seasforce_peak) %>%
  select(!c(omega,seasforce_peak,peak_week)) %>% 
  rename(`R0 (baseline)`=R0, 
         `age-dependence`=age_dep,`exposure-dependence`=exp_dep,
         `season width (weeks)`=seasforc_width_wks) %>% 
  group_by(accepted) %>% 
  slice_sample(n=2e3) %>% pivot_longer(!c(par_id,accepted))
# median values
median_parvals_accept = plot_partable_histogram %>% 
  group_by(accepted,name) %>% summarise(median_parval=median(value))

# KS significance test
KS_test_accept = plot_partable_histogram %>% filter(!name %in% "peak forcing (week)") %>% 
  # select(!const_delta) %>% pivot_longer(!c(par_id,`early off season`)) %>%
  group_by(name) %>% 
  summarise(min_val=min(value),max_val=max(value),
            p_val=ks.test(x=value[accepted],y=value[!accepted])$p.value,
            signif=p_val<0.01)

# plot distrib by jitter
ggplot(plot_partable_histogram %>% 
         filter(!grepl("maximal",name) & !grepl("baseline",name)),
       aes(x=accepted,y=value,color=accepted)) + 
  geom_jitter(width=0.4,alpha=1/4) + # geom_violin(fill=NA,show.legend=F) + 
  geom_boxplot(fill=NA,width=0.88,size=3/4,outlier.colour=NA,color="black") + # 
  facet_wrap(~name,scales="free_x",nrow = 4) + # scale_y_log10() +
  scale_color_manual(values=c("grey","blue"),
                     guide=guide_legend(override.aes=list(size=3))) + 
  geom_text(data=KS_test_accept %>% 
                  filter(!grepl("maximal",name) & !grepl("baseline",name)),color="black",
            aes(x=2/3,y=max_val*0.95,label=paste0("p=",signif(p_val,3),ifelse(signif,"**","")))) +
  xlab("") + ylab("parameter values") + theme_bw() + standard_theme +
  theme(strip.text=element_text(size=13),
        axis.text.x=element_text(size=13),axis.text.y=element_text(size=13),
        legend.position="top",legend.text=element_text(size=13),
        legend.title=element_text(size=13)) + coord_flip()
# save
ggsave(figs_folder,"param_distrib_accept_jitter.png",width=28,height=18,units="cm")

# CDF plot
ggplot(plot_partable_histogram %>% filter(!grepl("maximal",name) & !grepl("baseline",name)), 
       aes(x=value,color=accepted)) + 
  stat_ecdf(geom="step") + 
  geom_vline(data=median_parvals_accept %>% 
               filter(!grepl("maximal",name) & !grepl("baseline",name)),
             aes(xintercept=median_parval,color=accepted),
             linetype="dashed",size=1/2,show.legend=F) +
  geom_text(data=KS_test_accept %>% filter(!grepl("maximal",name) & !grepl("baseline",name)),
            aes(x=max_val*0.9,y=0.1,label=paste0("p=",signif(p_val,3),ifelse(signif,"**","")) ),
            color="black") +
  facet_wrap(~name,scales="free",nrow=3) +
  scale_color_manual(values=c("grey","blue")) + labs(color="accepted") + 
  xlab("") + ylab("density") + theme_bw() + standard_theme +
  theme(strip.text=element_text(size=15),axis.text.y=element_text(size=12),
        legend.position="top")
#
# ggsave("simul_output/2e4_parsets/crit_8_11/param_distrib_accept_CDF.png",width=28,height=18,units="cm")
ggsave("simul_output/2e4_parsets/crit_8_11/param_distrib_accept_CDF_peakR0.png",width=28,height=18,units="cm")

# plot correlations btwn params
library(GGally)
ggpairs(plot_partable_histogram %>% mutate(name=gsub("dependence","dep.",name),
              name=gsub("width","width\n",name), name=gsub("forcing","forcing\n",name)) %>%
          # filter(!grepl("peak",name)) %>%
          filter(!grepl("maximal",name) & !grepl("baseline",name)) %>% 
          pivot_wider(names_from=name) %>%
          select(!c(par_id)), # `R0 peak`,`peak forcing (week)`,`season width (weeks)`
      columns=2:(length(unique(plot_partable_histogram$name))-1), 
      aes(color=accepted,alpha=1/4)) + 
  scale_color_manual(values=c("grey","blue")) + scale_fill_manual(values=c("grey","blue")) +  
  theme_bw() + standard_theme + theme(strip.text=element_text(size=9))
# save
ggsave("simul_output/2e4_parsets/crit_8_11/param_correlations.png",width=28,height=18,units="cm")
# ggsave("simul_output/2e4_parsets/crit_8_11/param_correlations_peakR0.png",width=28,height=18,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# reduce the two parameters exp_dep and age_dep to their value along their 1st principal component

# distribution of accepted parameter sets in terms of age_dep-exp_dep
ggplot(partable_regular_dyn %>% 
         filter(par_id %in% sample(partable_regular_dyn$par_id,size=2e3))) + 
  geom_point(aes(x=age_dep,y=exp_dep),color="blue",alpha=1/2) +
  scale_x_continuous(limits=c(0,0.4)) + scale_y_continuous(limits=c(0,5/4)) + 
  theme_bw() + standard_theme

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Plot susceptibility as a function of age and exposure (SI FIG 3)

age_exp_dep_uniqvals <- list(exp_dep=seq(3/10,1.25,0.19),
                             age_dep=seq(1/15,1/3,1/15),age=1:11,exp=1:3) %>% 
  expand.grid %>% bind_rows %>% 
  mutate(suscept_unscaled=exp(-(exp_dep*exp+age_dep*age)))

age_exp_dep_uniqvals <- age_exp_dep_uniqvals %>% 
  mutate(const_delta=1/unlist(lapply(lapply(1:nrow(age_exp_dep_uniqvals), 
          function(n_p) { sapply(1:n_age,function(x) {
                  (1*exp(-age_exp_dep_uniqvals$exp_dep[n_p]*(1:3)))/(exp(age_exp_dep_uniqvals$age_dep[n_p]*x))})}),
          function(x) R0_calc_SIRS(C_m,x,rho,n_inf))),susc_scaled=suscept_unscaled*const_delta)

# PLOT
ggplot(age_exp_dep_uniqvals %>% 
         filter(exp_dep %in% seq(3/10,1.25,0.19)[c(1,3,5)] 
                & age_dep %in% seq(1/15,1/3,1/15)[c(1,3,5)]) %>%
         mutate(age_dep=round(age_dep,2)) %>%
         rename(`exposure-dependence`=exp_dep,`age-dependence`=age_dep) %>% 
         mutate(age=factor(rsv_age_groups$agegroup_name[age],
                           levels=unique(rsv_age_groups$agegroup_name))) )  + 
  geom_line(aes(x=age,color=factor(exp),group=exp,y=susc_scaled),size=1.06) + 
  facet_grid(`exposure-dependence`~`age-dependence`,
             labeller=labeller(`exposure-dependence`=label_both,
                               `age-dependence`=label_both)) + 
  scale_y_log10() + 
  labs(color="exposure") + xlab("age group") + ylab(expression(delta[exp]^(age))) +
  theme_bw() + standard_theme
  
# ggsave
ggsave(here(figs_folder,"age_exp_dep.png"),width=22,height=18,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Entire range of outputs by each parameter

source("fcns/create_results_fullscan_hosp.R")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# calculate ranges: fix 4 params to the median, 
# calculate the values at the min and max of the parameter that is varied
# take median parameter set, calculate the range in ONE parameter, 
# for ALL param sets and for the accepted parsets
source("fcns/create_df_fullrange_plots.R")
# creates the dataframes: output_ranges_full_scan | mean_age_shift_ranges

# create plots  (not included in manuscript)
source("fcns/create_fullrange_plots.R")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Sensitivity analysis: PRCC
# library(epiR)

# create dataframe of PRCC values
source("fcns/create_df_prcc.R")

# plot cumul age and peak
bind_rows(list_prcc_output) %>% 
  filter((!name %in% "kappa") & 
           !grepl("mean_age",output) & (!agegroup_broad %in% "5+y")) %>% 
  rename(`age`=agegroup_broad) %>% 
  mutate(name=factor(name,levels=unique(name))) %>%
ggplot() + facet_wrap(~output,scales = "free_x") + 
  geom_col(aes(y=est,x=name,group=age,fill=age),
           color="black",size=1/3,position=position_dodge(width=0.85),width=4/5) + 
  coord_flip() + xlab("") + ylab("PRCC") + 
  geom_vline(xintercept=1/2+1:5,size=1/3) + geom_hline(yintercept=0) +
  theme_bw() + standard_theme + manuscript_large_font_theme + 
  theme(panel.grid.major.y = element_blank())
# save
ggsave(here(figs_folder,"PRCC_cumul_peak_hosp_under5.png"),width=28,height=20,units="cm")

# plot mean age shift
bind_rows(list_prcc_output) %>% 
  filter((!name %in% "kappa") & 
            grepl("mean_age",output)) %>% rename(`age`=agegroup_broad) %>% 
  mutate(output=paste0(gsub("mean_age_shift_","average age of hospitalisation (",output),")"),
         output=ifelse(grepl("2020",output),
                       gsub("2020","2020/21",output),gsub("2021","2021/22",output)),
         output=ifelse(grepl("22",output),
                        gsub("\\)","",
                        gsub("average age of hospitalisation \\(","",output)), output),
         output=factor(output,levels=(unique(output))),
         name=factor(name,levels=unique(name)) ) %>%
ggplot() + 
  geom_col(aes(y=est,x=name,group=output,fill=output),
           color="black",size=1/3,position=position_dodge(width=0.85),width=4/5) +
  coord_flip() + xlab("") + ylab("PRCC") + labs(fill="") +
  geom_vline(xintercept=1/2+1:5,size=1/3) + geom_hline(yintercept=0) +
  theme_bw() + standard_theme + 
  manuscript_large_font_theme + theme(panel.grid.major.y=element_blank())
# save
ggsave(here(figs_folder,"PRCC_mean_age_hosp_under5.png"),width=28,height=20,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# outcomes by the value of parameters
# LOAD
# summ_broad_age_groups_byvalue <- read_csv(here(foldername,"summ_broad_age_groups_byvalue.csv"))
summ_broad_age_groups_byvalue <- results_fullscan_hosp %>% 
  filter(!agegroup_broad %in% "5+y" & 
           par_id %in% partable_regular_dyn$par_id & 
           epi_year>2019) %>%
  select(agegroup_broad,par_id,epi_year,exp_dep,age_dep,kappa,
         seasforc_width_wks,seasforce_peak,R0,omega,
         hosp_tot_norm,peak_hosp_norm,mean_age_cumul_hosp_shift) %>%
  mutate(age_exp_ratio=age_dep/exp_dep) %>%
  pivot_longer(!c(agegroup_broad,par_id,epi_year,hosp_tot_norm,
                  peak_hosp_norm,mean_age_cumul_hosp_shift)) %>%
  rename(parname=name,parvalue=value) %>% 
  relocate(c(parname,parvalue),.before=hosp_tot_norm) %>%
  mutate(epi_year_20_21_merged=
           ifelse(epi_year %in% c(2020,2021),"2020-21",epi_year)) %>% 
  # off-season outbreaks fall into epi-year 2020 -> 
  # take sum of 2020+2021 for cumulative and 20-21 epiyear maximum for peak
  group_by(agegroup_broad,par_id,epi_year_20_21_merged,parname) %>%
  summarise(parvalue=unique(parvalue),hosp_tot_norm=sum(hosp_tot_norm),
            peak_hosp_norm=max(peak_hosp_norm),
            mean_age_cumul_hosp_shift=max(mean_age_cumul_hosp_shift)) %>%
  group_by(parname) %>%
  mutate(par_bin=ntile(parvalue,10)) %>% relocate(par_bin,.after=parvalue) %>%
  pivot_longer(!c(agegroup_broad,par_id,epi_year_20_21_merged,
                  parname,parvalue,par_bin)) %>% rename(varname=name) %>%
  # calculate summary stats (median, ci50, ci95)
  group_by(agegroup_broad,epi_year_20_21_merged,parname,par_bin,varname) %>%
  summarise(parvalue=median(parvalue),mean_parvalue=mean(parvalue),
            median=median(value,na.rm=T),mean=mean(value,na.rm=T),
            ci50_l=quantile(value,probs=0.25,na.rm=T),
            ci50_u=quantile(value,probs=0.75,na.rm=T),
            ci95_l=quantile(value,probs=0.025,na.rm=T),
            ci95_u=quantile(value,probs=0.975,na.rm=T)) %>%
  group_by(parname,par_bin) %>% 
  mutate(parvalue=median(parvalue),mean_parvalue=mean(parvalue)) %>%
  filter(!(varname %in% "mean_age_cumul_hosp_shift" & !agegroup_broad %in% "<1y")) %>%
  mutate(agegroup_broad=ifelse(varname %in% "mean_age_cumul_hosp_shift",
                               "<5y",agegroup_broad)) %>%
  rename(epi_year=epi_year_20_21_merged)

# save
write_csv(summ_broad_age_groups_byvalue,here("repo_data/summ_broad_age_groups_byvalue.csv"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT summary stats disaggregated by parameter VALUES
# the figures produced by this loop include: 

# peak values are not generated from this dataframe because of the agegroup resolution being different,
# instead they are calculated from the dynamics below
source("fcns/create_summary_plots_narrow_age_groups.R")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Evaluating the DYNAMICS 

# plot pre-pandemic dynamics for all parameter sets (Figure 1): all age groups <5y, compare to data

# overlay plot (SARI-Watch data on all simulations, Figure 1B)
date_limits <- as.Date(c("2017-09-15","2020-04-01"))
# take paramsets with best LLH
subsample_par = (all_likelihoods %>% 
                   filter(name %in% "complete likelihood" & accepted & value<1500))$par_id
# simul_hosp_rate_weekly_under5_over65: contains both accepted and rejected parsets
# simul_hosp_rate_weekly_under5_over65_grad_relax: only accepted parsets, with gradual relaxation
hosp_plot_df <- simul_hosp_rate_weekly_under5_over65 %>%
  select(broad_age,par_id,date,simul_hosp_rate_100k) %>% ungroup() %>%
  filter(broad_age %in% "<5y") %>% 
  mutate(sel_par=ifelse(par_id %in% partable_regular_dyn$par_id,"accepted","rejected")) %>%
  # select best LLH params?
  filter((par_id %in% subsample_par) | sel_par %in% "rejected") %>%
  group_by(sel_par) %>% filter(par_id %in% sample(unique(par_id),size=100)) %>%
  filter(date>=date_limits[1] & date<=date_limits[2])

# medians across simulated (accepted) parameter sets
# simul_hosp_rate_weekly_under5_over65 OR simul_hosp_rate_weekly_under5_over65_grad_relax
median_weekly_pred <- simul_hosp_rate_weekly_under5_over65 %>% ungroup() %>%
  filter(par_id %in% hosp_plot_df$par_id[hosp_plot_df$sel_par %in% "accepted"]) %>%
  ungroup() %>% group_by(date,year_week,broad_age) %>%
  summarise(incid_hosp_med_val=median(simul_hosp_sum_full),
            incid_hosp_per100k=median(simul_hosp_rate_100k),
            incid_hosp_per100k_with_underrep=incid_hosp_per100k*ifelse(broad_age %in% "<5y",
                    under_report_factor_under5,under_report_factor_over65y) ) %>% distinct()

# data from SARI_watch
SARI_watch_under5y_hosp_rate <- left_join(
  read_csv(here(foldername,"SARI_watch_under5y_hosp.csv")) %>%
    mutate(year_week=gsub("-0","-",year_week)),
  median_weekly_pred %>% select(date,year_week) %>% ungroup() %>% distinct()) %>% 
  relocate(date,.before=year_week) 
# this data is under-reported, so simulations are scaled by under-reporting factor

# plot
ggplot() +
  geom_line(data=hosp_plot_df,
      aes(x=date,y=simul_hosp_rate_100k*under_report_factor_under5,
          group=par_id,alpha=sel_par,color=sel_par)) +
  scale_color_manual(values=c("blue","black")) + scale_alpha_manual(values=c(1/5,1/8)) +
  geom_point(data=SARI_watch_under5y_hosp_rate %>% 
               filter(date>date_limits[1] & date<date_limits[2]),
      aes(x=date,y=rate_under5yrs)) + # overlay data # ,shape=21,fill=NA,size=2
  geom_line(data=median_weekly_pred %>% 
              filter(date>date_limits[1] & date<date_limits[2] & broad_age %in% "<5y"),
      aes(x=date,y=incid_hosp_per100k_with_underrep),
      color="red",linetype="dashed",size=1.02) + # median simulation
  theme_bw() + standard_theme + 
  xlab("") + ylab("weekly hospitalisations <5y per 100.000 persons") + labs(alpha="",color="") +
  scale_x_date(expand=expansion(1/100,0),date_breaks="2 months") +
  scale_y_continuous(expand=expansion(0.01,0)) +
  theme(strip.text=element_text(size=15),axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=12),legend.text=element_text(size=11),
        legend.title=element_text(size=12),legend.position="none")
# SAVE
# ggsave(here(foldername,"prepandemic_dyn_compare_SARIwatch_under5_median_simul.png"),
#           width=28,height=16,units="cm")
ggsave(here(foldername,"prepandemic_dyn_compare_SARIwatch_under5_median_simul_LLH1500.png"),
       width=28,height=16,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# select params where off-season outbreak is earlier (as what was observed)

# calculate distance for 2020-2021
subsample_par = (all_likelihoods %>% 
                   filter(name %in% "complete likelihood" & accepted & value<2e3))$par_id
hosp_plot_df = simul_hosp_rate_weekly_under5_over65_grad_relax %>%
  select(broad_age,par_id,date,simul_hosp_rate_100k) %>% ungroup() %>%
  filter(broad_age %in% "<5y") %>% 
  mutate(sel_par=ifelse(par_id %in% partable_regular_dyn$par_id,"accepted","rejected")) %>%
  # select best LLH params?
  filter((par_id %in% subsample_par) | sel_par %in% "rejected")

# we calculate euclidean distance from observed hosp rate in 2021-22
dist_hosp_2021_22 = left_join(
    hosp_plot_df %>% filter(date>=as.Date("2021-01-01") & 
                              date<=max(SARI_watch_under5y_hosp_rate$date)),
    SARI_watch_under5y_hosp_rate,by=c("date","year_week","broad_age")) %>%
  mutate(abs_dist=abs(rate_under5yrs-simul_hosp_rate_100k),sqrd_dist=abs_dist^2) %>%
  group_by(par_id,broad_age) %>% summarise(mean_abs_dist=mean(abs_dist),mean_sqrd_dist=mean(sqrd_dist))

early_off_season = (dist_hosp_2021_22 %>% 
          filter(mean_sqrd_dist<=quantile(dist_hosp_2021_22$mean_sqrd_dist,probs=0.1)))$par_id

# plot param sets where off-season outbreak in 2021 is earlier (as it was in reality)
ggplot() +
  geom_line(data=hosp_plot_df %>% 
              filter(par_id %in% early_off_season | sel_par %in% "rejected") %>%
              group_by(sel_par) %>% 
              filter(par_id %in% sample(unique(par_id),size=length(early_off_season)) & 
                       date<=as.Date("2022-04-01")),
            aes(x=date,y=simul_hosp_rate_100k*under_report_factor_under5,
                group=par_id,alpha=sel_par,color=sel_par)) +
  scale_color_manual(values=c("blue","black")) + 
  scale_alpha_manual(values=c(1/5,1/8)) +
  geom_point(data=SARI_watch_under5y_hosp_rate %>% 
               filter(date>date_limits[1] & date<date_limits[2]),
             aes(x=date,y=rate_under5yrs)) + # overlay data # ,shape=21,fill=NA,size=2
  theme_bw() + standard_theme + 
  xlab("") + ylab("weekly hospitalisations <5y per 100.000 persons") + labs(alpha="",color="") +
  scale_x_date(expand=expansion(1/100,0),date_breaks="2 months") + 
  scale_y_continuous(expand=expansion(0.01,0)) +
  theme(strip.text=element_text(size=15),axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=12),legend.text=element_text(size=11),
        legend.title=element_text(size=12),legend.position="none")
# SAVE
ggsave(here(figs_folder,
            "prepandemic_dyn_compare_SARIwatch_under5_median_simul_early_offseason_gradrelax.png"),
        width=28,height=16,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# summary by year
epi_year_week_start <- isoweek(as.Date("2020-06-01"))
# normalise by pre-pandemic epi-years 2017 and 2018
summ_dyn_peak_cumul_meanage <- simul_hosp_rate_weekly_all_broad_agegroups %>%
  filter(par_id %in% partable_regular_dyn$par_id) %>%
  mutate(epi_year=ifelse(week(date)>=40,year(date),year(date)-1),
         epi_year_june=ifelse(week(date)>=epi_year_week_start,year(date),year(date)-1)) %>% ungroup() %>%
  filter(epi_year>2016) %>%
  pivot_longer(cols=c(epi_year,epi_year_june)) %>% rename(epi_year_wk_start=name,epi_year=value) %>%
  mutate(epi_year_wk_start=ifelse(grepl("june",epi_year_wk_start),epi_year_week_start,40)) %>%
  group_by(epi_year,epi_year_wk_start,par_id,agegroup) %>%
  summarise(peak_hosp=max(simul_hosp_sum_full), cumul_hosp=sum(simul_hosp_sum_full),
            peak_yday=mean(yday(date[simul_hosp_sum_full %in% max(simul_hosp_sum_full)]))) %>%
  pivot_longer(!c(epi_year,epi_year_wk_start,par_id,agegroup)) %>%
  group_by(par_id,agegroup,name,epi_year_wk_start) %>% # calculate pre-pandemic mean
  # because of some biennial patterns need to normalise 
  # by the right pre-pandemic year (2022/2018, 2021/2017)
  mutate(norm_value=case_when(
    grepl("hosp",name) ~ value/value[epi_year==ifelse(epi_year %% 2==0,2018,2017)],
    !grepl("hosp",name) ~ value-value[epi_year==ifelse(epi_year %% 2==0,2018,2017)])) %>%
  relocate(par_id,.before=epi_year)
# add mean age
if (!any(grepl("age",unique(summ_dyn_peak_cumul_meanage$name)))) {
  summ_dyn_peak_cumul_meanage = 
    bind_rows(summ_dyn_peak_cumul_meanage,
              summ_dyn_peak_cumul_meanage %>% 
                group_by(par_id,epi_year_wk_start,name,epi_year) %>%
                    mutate(mean_age_under5y=ifelse(agegroup==1 & 
                              name %in% "cumul_hosp",sum(value[agegroup==1]*1/2+
                              value[agegroup==1]*1.5+value[agegroup==1]*3.5)/sum(value),NA)) %>% 
                ungroup() %>% select(epi_year,epi_year_wk_start,par_id,mean_age_under5y) %>% 
                filter(!is.na(mean_age_under5y)) %>% distinct() %>% 
                group_by(par_id,epi_year_wk_start) %>%
                mutate(mean_age_shift=
                         mean_age_under5y-
                         mean_age_under5y[epi_year==ifelse(epi_year %% 2==0,2018,2017)]) %>%
                rename(value=mean_age_under5y,
                       norm_value=mean_age_shift) %>% mutate(name="mean_age_under5y") )
  }

# left-join with partable to have input parameters, calculate summary stats by quantiles of parameters
# summ_max_incid_seas_length_byvalue
summ_dyn_peak_cumul_meanage_byparvalue <- left_join(
  summ_dyn_peak_cumul_meanage %>% filter(epi_year>2019),
  partable_regular_dyn %>% 
    select(!const_delta) %>% 
    rename(waning=omega,age_exp_PC1=PC1) %>% 
    mutate(age_exp_ratio=age_dep/exp_dep) %>% 
    pivot_longer(!par_id) %>% rename(parname=name,parvalue=value), by="par_id") %>% 
  relocate(c(name,value,norm_value),.after=last_col()) %>% 
  relocate(par_id,.before=epi_year) %>% rename(varname=name) %>%
  group_by(parname) %>% mutate(par_bin=ntile(parvalue,10)) %>% 
  relocate(par_bin,.after=parvalue) %>% 
  # calculate summary stats (median, ci50, ci95)
  group_by(agegroup,epi_year,epi_year_wk_start,parname,par_bin,varname) %>%
  summarise(parvalue=median(parvalue),mean_parvalue=mean(parvalue),
            median=median(norm_value),mean=mean(norm_value),
            ci50_l=quantile(norm_value,probs=0.25),
            ci50_u=quantile(norm_value,probs=0.75),
            ci95_l=quantile(norm_value,probs=0.025),
            ci95_u=quantile(norm_value,probs=0.975))

##############################################################
# plot PEAK VALUE of cases/hospitalisations (+ season length), disaggregated by parameter values

# Figure 3B, SI 7-10
sel_pars <- unique(summ_dyn_peak_cumul_meanage_byparvalue$parname)
sel_varnames <- unique(summ_dyn_peak_cumul_meanage_byparvalue$varname)
start_year=2020
# this file also saves the plots!
source("fcns/create_summ_plots_from_dyn.R")

##############################################################
# calculate PRCC values wrt parameters

source("fcns/create_df_prcc_from_dynamics.R")
# creates the list: list_prcc_dyn

# plot PRCC on cumulative/peak hospitalisations, 
# timing of peak (yday relative to pre-NPI), mean age of hosp
source("fcns/create_prcc_plot_from_dyn.R")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Evaluating DYNAMICS as a function of parameters (relative timing of peak and duration/shape of season)

# all_dynamics_accepted <- read_csv(here(foldername,"all_dynamics_accepted_2e3.csv")) %>%
#  mutate(date=as.Date(start_date_dyn_save)+t-min(t)) %>%
#   filter(date>=as.Date("2018-06-01") & date<as.Date("2024-10-01"))

epi_year_week_start=23
# if dataframe starts with Autumn 2017!!!!
startweek_startyr=isoweek(min(simul_hosp_rate_weekly_all_broad_agegroups_grad_relax$date))-23
start_yr=isoyear(min(simul_hosp_rate_weekly_all_broad_agegroups_grad_relax$date))

dyn_hosp_weekly_norm_to_peak <- simul_hosp_rate_weekly_all_broad_agegroups %>% ungroup() %>% 
      filter(par_id %in% partable_regular_dyn$par_id & agegroup %in% 1:3) %>% 
      rename(value=simul_hosp_sum_full) %>%
      mutate(epi_year=ifelse(week(date)>=epi_year_week_start,
                             year(date),year(date)-1)) %>% 
      group_by(epi_year,par_id,agegroup) %>% 
      mutate(epi_week=ifelse(epi_year==start_yr,
                             row_number()+startweek_startyr,row_number())) %>%    
      ungroup() %>% group_by(agegroup,par_id) %>%
      mutate(peak_value_2017=max(value[epi_year %in% 2017]),
             peak_value_2018=max(value[epi_year %in% 2018])) %>% # peak level of selected year
      group_by(agegroup,par_id) %>% 
      # peak week of year of reference  
      mutate(peak_week_2017=min(epi_week[value==peak_value_2017 & (epi_year %in% 2017)]),
             peak_week_2018=min(epi_week[value==peak_value_2018 & (epi_year %in% 2018)])) %>% 
      # distance from peak week of sel year    
      ungroup() %>% 
  mutate(peak_week_distance=ifelse(epi_year %% 2==0,
                               epi_week-peak_week_2018,epi_week-peak_week_2017)) %>% 
      mutate(value_norm=ifelse(epi_year %% 2==0,
                               value/peak_value_2018,value/peak_value_2017))

# bin parameter values 
binned_partable_regular_dyn <- partable_regular_dyn %>% 
  select(!c(const_delta,seasforc_width_wks,peak_week)) %>%
  filter(par_id %in% unique(dyn_hosp_weekly_norm_to_peak$par_id)) %>%
  mutate(age_exp_ratio=age_dep/exp_dep) %>%
  pivot_longer(!par_id) %>% group_by(name) %>% 
  mutate(par_bin=ntile(value,8)) %>% 
  group_by(name,par_bin) %>% mutate(median_parbin=median(value))

#######
# calculate statistics at time points relative to peak week of selected reference year,
# for each bin of the selected parameters
# this takes about 1 min
source("fcns/create_dyn_wk_hosp_norm_time_norm_value.R")

# this creates the dataframe `summ_dyn_all_parsets_broad_age_relat_time`

##############################################################
# Plot dynamics faceted by the age vs immunity dependence of susceptibility (Figure 4)

# FIGURE 4 and SI FIG : plot faceted by years, color-coded by parameter values
source("fcns/create_norm_dyn_plots.R")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# DATA PLOTS (for SI)
# Plot UK RSV case data, calculate seasonal concentration

season_weeks=c(13,40)
resp_detects_weekly_all_age <- 
  read_csv(here(foldername,"Respiratory_viral_detections_by_any_method_UK.csv")) %>% 
  mutate(year_week=factor(paste0(Year,"-",Week),unique(paste0(Year,"-",Week))), 
         RSV_rolling_av=roll_mean(RSV,n=7,align="center",fill=NA) ) %>%
  select(-(contains("virus")|contains("flu"))) %>% 
  mutate(epi_year=ifelse(Week>=season_weeks[2],Year+1,Year)-min(Year)+1) %>% 
  group_by(epi_year) %>% mutate(perc_yearly=RSV/sum(RSV)) %>% group_by(epi_year) %>% 
  mutate(season_share=sum(perc_yearly[Week>=season_weeks[2] | Week<=season_weeks[1]]),
         on_off_season=ifelse(findInterval(Week,season_weeks+c(1,0))==1,"off","on"))

# SI Table 5
resp_detects_weekly_all_age_means_shares <- left_join(
  resp_detects_weekly_all_age %>% 
    group_by(epi_year,on_off_season) %>% 
    summarise(mean_on_off=mean(RSV,na.rm=T),
              cal_year=paste0(unique(Year),collapse="_")) %>% 
    filter(epi_year>1&epi_year<8) %>%
    mutate(cal_year=ifelse(on_off_season %in% "off",
                           paste0(as.numeric(cal_year)-1,"_",
                           as.numeric(cal_year)),cal_year)) %>% 
    pivot_wider(names_from=on_off_season,values_from=mean_on_off,names_prefix="mean_"),  
  resp_detects_weekly_all_age %>% 
    group_by(epi_year) %>% 
    summarise(season_share=unique(season_share)),by="epi_year" ) %>%
  mutate(on_off_ratio=round(mean_on/mean_off,1),
         mean_on=round(mean_on,1),mean_off=round(mean_off,1),
         season_share=round(season_share,3))

### ### ### ### ### ### ### ### ### ### ###
# SI Figure 1 (RSV case data with age)
resp_virus_data_uk_tidy <- 
  read_csv(here(foldername,"Respiratory_viral_detections_by_any_method_UK_Ages.csv")) %>% 
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

# DATAMART data
datamart_0_5y_2012_2022 = read_csv("data/datamart_0_5y_2012_2022.csv")
# seasonal concentration (compare to `resp_detects_weekly_all_age_means_shares`)
datamart_seas_conc = datamart_0_5y_2012_2022 %>% #
  mutate(epi_year=ifelse(week_number>=40,year(date),year(date)-1)) %>%
  filter(epi_year>2012 & epi_year<2020 & !is.na(epi_year)) %>% group_by(epi_year) %>% 
  summarise(in_season_num=sum(num_positives_0_5y[week_number>=40 | week_number<=13]),
            off_season_num=sum(num_positives_0_5y[week_number<40 & week_number>13],na.rm=T),
            in_season_posit=sum(positivity_0_5y[week_number>=40 | week_number<=13]),
            off_season_posit=sum(positivity_0_5y[week_number<40 & week_number>13],na.rm=T),
            in_season_share_num=in_season_num/(in_season_num+off_season_num),
            in_season_share_posit=in_season_posit/(in_season_posit+off_season_posit))

# SARI-Watch seasonal concentration
SARI_watch_all_hosp = read_csv("repo_data/SARI_watch_all_hosp.csv") 
SARI_watch_all_hosp %>% filter(epi_year<2020) %>% group_by(epi_year) %>% 
  summarise(in_season_hosp=sum(rate_per_100000[week_number>=40 | week_number<=13],na.rm=T),
            off_season_hosp=sum(rate_per_100000[week_number<40 & week_number>13],na.rm=T),
            in_season_share=in_season_hosp/(in_season_hosp+off_season_hosp))

# end of plots
# below are calculations for the range of outputs (cumul and peak hospitalisations) by parameters

# concentration of cases
subset(resp_virus_data_uk_tidy,name %in% "RSV")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# annual burden of hospitalisation

# data
SARI_watch_under5y_hosp_rate %>% 
  mutate(week=isoweek(date),epi_year=ifelse(week>23,isoyear(date),isoyear(date)-1)) %>% 
  group_by(broad_age,epi_year) %>% summarise(sum(rate_under5yrs))
# simul
simul_hosp_rate_weekly_under5_over65_from2017 %>% 
  mutate(week=isoweek(date),epi_year=ifelse(week>23,isoyear(date),isoyear(date)-1)) %>% 
  ungroup() %>% group_by(broad_age,epi_year,par_id) %>% 
  summarise(annual_burden=sum(simul_hosp_rate_100k)) %>%
  group_by(broad_age,epi_year) %>% 
  summarise(median(annual_burden),min(annual_burden),max(annual_burden))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
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
  summarise(mean=mean(value_norm), median=median(value_norm),
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



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### THE END ¯\_(ツ)_/¯
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PS: 
# To run individual simulations and plot results go to "indiv_simul.R" 
#

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# first principal component of age_dep - exp_dep
# pred_pca <- data.frame(
#   predict(object=prcomp(partable_regular_dyn %>% select(c(exp_dep,age_dep))),
#           newdata=partable_regular_dyn %>% select(c(exp_dep,age_dep))),
#   K_exp=partable_regular_dyn$exp_dep,K_age=partable_regular_dyn$age_dep,
#   par_id=partable_regular_dyn$par_id)


# PC1 of age_dep exp_dep
# 
# binned_pc1 = pred_pca %>% 
#   mutate(PC1_bin=findInterval(PC1,vec=seq(-0.65,0.35,by=0.1))) %>% 
#   group_by(PC1_bin) %>% # 
#   summarise(PC1=mean(PC1),K_exp=mean(K_exp),K_age=mean(K_age)) %>% 
#   pivot_longer(!c(PC1,PC1_bin))
# 
# # SI Figure 2: linear relationship between K_age and K_exp for selected param sets
# ggplot(pred_pca %>% pivot_longer(!c(PC1,PC2,par_id)), aes(x=PC1,y=value)) + 
#   geom_point(aes(color=name),alpha=1/2) + 
#   # geom_smooth(aes(group=name,color=name),fill=NA,method='lm',color="black") + 
#   geom_line(data=binned_pc1) + facet_wrap(~name,scales = "free_y") +
#   scale_x_continuous(breaks=((-3):2)/4,limits=c(-3/4,1/2)) +
#   xlab(expression(kappa)) + ylab(expression(paste(kappa[exp],", ",kappa[age]))) + 
#   labs(color="",fill="") + theme_bw() + standard_theme + 
#   theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
#         legend.title=element_text(size=16),legend.text=element_text(size=16))


# PC1 capture ~all variation in k_exp, 
# but there's residual variation in k_Age, how different are these parsets?
# sel_pars_PC1_high = pred_pca %>% filter(PC1>0.25 & (K_age<0.1 | K_age>0.3)) %>% 
#   mutate(age_dep_low_high=ifelse(K_age<0.1,"low","high")) %>% 
#   group_by(age_dep_low_high) %>%
#   filter(par_id %in% sample(par_id,size=11))
# # difference in dynamics? yes, dynamics quite 
# different at high PC1 value but different k_age values
# left_join(sel_pars_PC1_high %>% 
#             select(c(par_id,age_dep_low_high,K_age)), all_dynamics_accepted) %>% 
#   filter(!is.na(agegroup) & 
#            date>as.Date("2018-09-01") & date<as.Date("2020-04-01")) %>%
#   ggplot(aes(x=date,y=value,group=par_id,color=factor(age_dep_low_high))) + 
#   scale_x_date(expand=expansion(0.01,0)) +
#   geom_line(alpha=1/2) + facet_wrap(~agegroup) + theme_bw() + standard_theme

# partable <- read_csv(here(foldername,"partable.csv"))
# PCA on full param table
# pred_pca_FULL <- data.frame(
#   predict(object=prcomp(partable %>% select(c(exp_dep,age_dep))), 
#           newdata=partable %>% select(c(exp_dep,age_dep))),
#   exp_dep=partable$exp_dep, age_dep=partable$age_dep, par_id=partable$par_id) %>% 
#   rename(kappa=PC1) %>%
#   group_by(age_dep,exp_dep) %>% 
#   summarise(kappa=unique(kappa))
# 

# suscept_sel_parsets <- left_join(bind_rows(lapply(sample(1:nrow(pred_pca),size=2e3), 
#           function(x) age_exp_dep_uniqvals %>% 
#                 select(age,exp) %>% distinct() %>%
#                 mutate(suscept_unscaled=exp(-(pred_pca$K_exp[x]*exp+pred_pca$K_age[x]*age)),
#                       `exposure-dependence`=pred_pca$K_exp[x],
#                       `age-dependence`=pred_pca$K_age[x],PC1=pred_pca$PC1[x],
#                       par_id=pred_pca$par_id[x]) )), 
#           partable_regular_dyn %>% select(par_id,const_delta) ) %>% 
#   mutate(susc_scaled=suscept_unscaled*const_delta) %>%
#   mutate(age=factor(rsv_age_groups$agegroup_name[age],
#                     levels=unique(rsv_age_groups$agegroup_name)),
#          PC1_grouped=findInterval(PC1,seq(min(summary(pred_pca$PC1)),
#                                           max(summary(pred_pca$PC1)),by=1/10)) ) %>%
#   group_by(PC1_grouped) %>% 
#   mutate(PC1_grouped=round(mean(PC1),1)) %>% 
#   group_by(PC1_grouped,age,exp) %>% 
#   summarise(mean_val=mean(susc_scaled),med_val=median(susc_scaled),
#             ci50_low=quantile(susc_scaled,c(0.25,0.75))[1],
#             ci50_up=quantile(susc_scaled,c(0.25,0.75))[2],
#             ci95_low=quantile(susc_scaled,c(0.025,0.975))[1],
#             ci95_up=quantile(susc_scaled,c(0.025,0.975))[2]) %>%
#   rename(exposure=exp)
# # color palette
# colorpal <- colorRampPalette(colors=c("blue","grey","red"))(
#               length(unique(suscept_sel_parsets$PC1_grouped)))
# ##################################################
# # library(plyr)
# label_parseall <- function(variable, value) {
#   plyr::llply(value, function(x) parse(text=paste(variable,x,sep = "==")))}
# 
# # plot
# ggplot(suscept_sel_parsets %>%
#          rename(kappa=PC1_grouped) %>% ungroup() %>% 
#          mutate(kappa_num=as.numeric(factor(as.character(kappa)))) %>% 
#          filter(kappa_num %in% c(1,3,5,7,9,10)),
#        aes(x=age,color=factor(exposure),group=exposure,fill=factor(exposure))) +
#   geom_line(aes(y=med_val),size=1.06) + 
# geom_ribbon(aes(ymin=ci50_low,ymax=ci50_up),color=NA,alpha=0.2) +
#   facet_wrap(~kappa,labeller=label_parseall) + labs(color="exposure",fill="exposure") +
#   scale_y_log10(breaks=unlist(lapply(seq(-4,0,1/2), function(x) round(10^x,abs(x))))) + 
#   theme_bw() + standard_theme + manuscript_large_font_theme + 
# xlab("age group") + ylab(expression(delta[exp]^(age)))
# # save
# # ggsave(here(foldername,"suscept_by_dep_level.png"),width=25,height=20,units="cm")
