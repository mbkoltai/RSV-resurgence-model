# install.packages("devtools"); library("devtools")
# install_github("SineadMorris/shinySIR")

# ode solving, maximum likelihood, rcpp
# contact data from https://bisaloo.github.io/contactdata/index.html (Prem 2017)
# functions
rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# library(contactdata); library(fitdistrplus);  library(bbmle); library(Rcpp); library(GillespieSSA)
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","plotly","Rcpp","zoo","lubridate","tsibble","qs"),library,character.only=TRUE)
source('fcns/RSV_model_functions.R')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SET PARAMETERS --------------------------------------------------------
### ### ### ### ### ### ### ###
# constant parameters
# selected country
country_sel="United Kingdom"
# time resolution (in days)
elem_time_step=0.5
# population data
standard_age_groups <- fun_cntr_agestr(country_sel,"2020",seq(0,75,5),c(seq(4,74,5),99))
# RSV age groups (population data from wpp2019)
rsv_age_groups_low=c(0,0.5,1,1.5, 2,3,4, 5,10,15, 20); rsv_age_group_sizes=c(rep(0.4,4),rep(0.9,3),rep(4,3),79)
rsv_age_groups=fun_rsv_agegroups(standard_age_groups,rsv_age_groups_low,rsv_age_group_sizes)
# number of age groups, reinfections and variables (S,I,R)
n_age=nrow(rsv_age_groups); varname_list=c('S','I','R'); n_compartment=length(varname_list); n_inf=3
dim_sys=n_age*n_compartment*n_inf; n_days_year=365
# query variables: fun_sub2ind(1:3,11,"R",varname_list,n_age,n_inf)
# force of infection terms
# linear indices of the I & S variables
l_inf_susc=fun_inf_susc_index_lists(n_age,n_inf,varname_list); inf_vars_inds=l_inf_susc[[1]]; susc_vars_inds=l_inf_susc[[2]]
# CONTACT MATRIX
# contact matrix from covidm ("home","work","school","other")
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
list_contmatrs=fun_covidm_contactmatrix(country_sel="Germany",currentdir_path=currentdir_path,cm_path=cm_path)
C_m_full=Reduce('+',list_contmatrs)
# make matrix symmetric
C_m_fullrecipr=fun_recipr_contmatr(C_m_full,age_group_sizes=standard_age_groups$values)
# create for our age groups
C_m=fun_create_red_C_m(C_m_fullrecipr,rsv_age_groups,
      orig_age_groups_duration=standard_age_groups$duration,orig_age_groups_sizes=standard_age_groups$values)
# make it reciprocal for the larger group
C_m=fun_recipr_contmatr(C_m,age_group_sizes=rsv_age_groups$value)
# bc of reinfections we need to input contact matrix repeatedly
contmatr_rowvector=t(do.call(cbind, lapply(1:nrow(C_m), function(x){diag(C_m[x,]) %*% matrix(1,n_age,n_inf)})))
# build kinetic matrix
# WANING (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j
omega=1/150 # 1/runif(1,60,200)
# RECOVERY
rho=1/7 # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,rsv_age_groups)
# BIRTH RATE into S_1_1 (Germany 2019: 778e3 births)
birth_rate=713e3/365; birth_term=matrix(c(birth_rate,rep(0,dim_sys-1)),dim_sys,1) # 0.01
### ### ### ### ### ### ### ###
# variable parameters  --------------------------------------------------------
# SUSCEPTIBILITY (normalised by age group sizes, for infection terms ~ delta*(I1+I2+...+In)*S_i/N_i)
# agedep_fact determines strength of age dependence, if agedep_fact>1, decreasing susceptibility with age
agedep_fact=1; delta_primary=4*c(0.1,0.02,0.004) # rep(0.15,3)
delta_susc=sapply(1:n_age, function(x) {delta_primary/((agedep_fact^(x-1))*rsv_age_groups$value[x])})
delta_susc_prop=delta_susc*matrix(rep(rsv_age_groups$value,3),nrow=3,byrow=T)
### PLOT susceptibility: 
fcn_plot_suscept_table(fcn_suscept_agedeptable(rsv_age_groups,delta_susc,n_inf)) # ggsave("simul_output/suscept.png",width=32,height=22,units="cm")
# calculate R0 (at max seasonal forcing=1)
R0_calc_SIRS(C_m,delta_susc_prop,rho,n_inf) # scale_f=3; delta_susc_prop=delta_susc_prop*scale_f; delta_susc=delta_susc*scale_f
####
# DURATION of SIMULATION
n_years=15.25; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=elem_time_step)
# seasonal forcing (above baseline level of 1) | npi_reduc_strength: reduction from baseline 
# shutdown season (if x --> (x+1)th season shut down) | preseas_npi_on/postseas_npi_off: on/off NPI week before/after season onset
forcing_strength=3; npi_reduc_strength=0.3; npi_year=round(n_years-5); preseas_npi_on=0; postseas_npi_off=6
g(forcing_vector_npi,shutdwn_lims,seas_force,seas_lims) %=% fun_shutdown_seasforc(timesteps,elem_time_step=0.5,
      forcing_strength,npi_reduc_strength,npi_year,peak_week=46,season_width=4,preseas_npi_on,postseas_npi_off,n_prec=0.01,n_sd=2) 
# set seas lims from UK data: peak is weeks 49/50, on/off is 41,11
seas_lims$on=floor(seas_lims$on)+41*7/365; seas_lims$off=round(seas_lims$off)+11*7/365
# PLOT seasonal forcing with NPI
fcn_plot_seas_forc(timesteps,seas_force,forcing_vector_npi,shutdwn_lims,seas_lims) 
# SAVE: ggsave(paste0("simul_output/NPI_y",npi_year,"_on",preseas_npi_on,"w_off",postseas_npi_off,"w.png"),units="cm",height=10,width=20)

# INITIAL CONDITIONS. # introduce stationary state as init state? | init_set: "previous" or anything else
prev_or_new=c("previous","fromscratch")[1]; filename="simul_output/df_ode_solution_UK_long.RDS"
# init_cond_src: "file" or "output"? 
initvals_sirs_model=fcn_set_initconds(init_set=prev_or_new,init_cond_src="output",ode_solution,init_seed=10,seed_vars="all",filename)
# manually choose a timepoint from prev simul
# initvals_sirs_model=matrix(as.numeric(round(ode_solution[min(which(round(timesteps/365,3)==19.5)),2:ncol(ode_solution)])))

### integrate ODE --------------------------------------------------------
params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,elem_time_step,delta_susc)
# interpolation fcns for seas forcing & extern introds
approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*10))
tm<-proc.time(); ode_solution<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params); round(proc.time()-tm,2)
# reshape data | # check size of objs: fcn_objs_mem_use(1)
g(ode_solution,df_ode_solution_tidy) %=% fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups,neg_thresh=-1e-3)
####
# PLOT how age group totals change (initial vs final popul sizes: fun_agegroup_init_final_pop(df_ode_solution_tidy))
fcn_plotagegroup_totals(df_ode_solution_tidy,scale_val=c('free_y','fixed')[1])
# ggsave(paste("simul_output/agegroup_totals_",scale_val,"yscale.png",sep=""),width=28,height=16,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Plot time course --------------------------------------------------------
xval_lims=c(npi_year-1.27,npi_year+4.5); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2]); agegr_lim=8
xval_breaks=seq(0,max_time/365,by=1/4)
# PLOT
# tags: y-axis fixed/free | facet by AGE/(age&INFECTION) | abs values/fraction
# k=7 --> (free scale, infects separate, absval). k=8 --> (free, infects separate, fractional)
g(scale_val,facet2tag,value_type,y_axis_tag,ncol_val,facet_formula,foldername,caption_txt,subtitle_str,timecourse_filename) %=% 
    fun_tcourse_plottags(k=8,nval=2,rval=3,n_inf,n_age,colvar="age",agegr_lim,delta_susc_prop,delta_primary,
          preseas_npi_on,postseas_npi_off,npi_reduc_strength,forcing_strength)
# PLOT time course absolute values  --------------------------------------------------------
fcn_plot_allcases_absval_stackedinfs(df_ode_solution_tidy,value_type,x_lims=c(0,13),t_subset=7,agegrlim=11,ncolval=3,y_axis_tag,
                                     scale_val,seas_lims_plot,shutdwn_lims,xval_breaks,subtitle_str,caption_txt)
## SAVE
ggsave(timecourse_filename, width=32,height=28,units="cm")

# PLOT time course normalised by max per age group  --------------------------------------------------------
xval_lims=c(npi_year-1.27,npi_year+3.25); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2])
fcn_plot_norm_max_byage(df_ode_solution_tidy,xval_lims,seas_lims_plot,shutdwn_lims,xval_breaks,subtitle_str,caption_txt)
# SAVE
ggsave(gsub(".png","_byage_maxnorm.png",timecourse_filename), width=30,height=20,units="cm")

### PLOT time course faceted by infection --------------------------------------------------------
xval_lims=c(npi_year-1.27,npi_year+3.25); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2])
fcn_plot_norm_max_byage(df_ode_solution_tidy,xval_lims,seas_lims_plot,shutdwn_lims,xval_breaks,subtitle_str,caption_txt)
# SAVE
ggsave(gsub(".png","_byinf_maxnorm.png",timecourse_filename), width=30,height=20,units="cm")

### share of infections --------------------------------------------------------
xval_lims=c(npi_year-2.27,npi_year+3.25); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2])
fcn_plot_share_infs_agefaceted(df_ode_solution_tidy,xval_lims,seas_lims_plot,n_aver=28,shutdwn_lims,xval_breaks,subtitle_str,caption_txt)
# SAVE
ggsave(gsub(".png","_byinf_share.png",timecourse_filename), width=32,height=18,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# dependence on age (>1) | dependence on prev exposure
agedep_fact=1; clinfract_expos=rep(0.4,n_inf) # c(0.8,0.4,0.2) # rep(0.75,3) # c(0.5,0.125,0.05) # rep(0.5,3)
clin_fract_age_exp=fcn_clin_fract_age_exp(agedep_fact,clinfract_expos,rsv_age_groups,n_age,n_inf)
# plot: fcn_plot_clinfract_table(clin_fract_age_exp)
### CALCULATE SYMPTOMATIC CASES & sum of 1,2,3rd infections -------------------------------
df_ode_sol_cases_sum=fun_symptomcases_table(df_ode_solution_tidy,clin_fract_age_exp,c("infection","agegroup"),seas_lims)

# PLOT SYMPTOMATIC CASES -----------------------------------
# k_plot: free/fixed | fraction/cases
g(scale_val,y_axis_tag,plotvar,foldername,symptcases_filename,caption_txt,subtitle_str) %=% fun_sumcase_plot_tags(n_val=2,r_val=2,k_plot=2,
        clin_fract_age_exp,delta_susc_prop,preseas_npi_on,postseas_npi_off,npi_reduc_strength,forcing_strength)
xval_lims=c(npi_year-1.5,npi_year+3.5); seas_lims_plot=subset(seas_lims,on>xval_lims[1]&off<xval_lims[2]) %>% pivot_longer(col=!season)
# plot
fcn_plot_symptomcases_agefacet(df_ode_sol_cases_sum,xval_lims,y_plotvar="symptom_cases",y_axis_tag,scale_val,seas_lims_plot,
                               shutdwn_lims,xval_breaks,n_col=3,n_agelim=9,subtitle_str,caption_txt)
# SAVE
# full_filename=gsub('.png','_reducedresurg.png',full_filename)
ggsave(symptcases_filename,width=31,height=22,units="cm")

### plot symptom cases stacked as area plot -----------
xval_lims=c(npi_year-0.5,npi_year+4.5); seas_lims_plot=subset(seas_lims,on>xval_lims[1]&off<xval_lims[2]) %>% pivot_longer(col=!season)
fcn_plot_symptomcases_agestacked(df_ode_sol_cases_sum,xval_lims,plotvar,y_axis_tag,seas_lims_plot,shutdwn_lims,
                                 xval_breaks,subtitle_str,caption_txt)
# SAVE
ggsave(gsub("symptomcases","symptomcases_areaplot",symptcases_filename),width=31,height=22,units="cm")

### Age distrib before and after shutdown ----------------------------
season_peaks_AUC=fcn_seas_agedistrib(df_ode_sol_cases_sum,max_time,timesteps,seas_case_threshold=1e3) %>%
  pivot_longer(cols=!c(agegroup,agegroup_name,season,agegr_size,mean_age,pre_post_shtdn))
# PLOT
plot_season_agedist=fcn_agedistrib_calc_season(season_peaks_AUC,seaslims=c(npi_year-1,npi_year+4),selvar="auc_case",
                                                 agegr_merge_min=which(rsv_age_groups$age_low==5)-1,rsv_age_groups)
if (length(unique(plot_season_agedist$pre_post_shtdn))==3) {colorvals=c(2,1,3)} else {colorvals=c(2,3)}
# bar plot
fcn_plot_agedistrib_perseas(plot_season_agedist,nrowval=2,yexpval=c(0,0.02),xval_lims,subtitle_str,caption_txt,textsize=3.5)
# SAVE
agedistr_filename=fcn_agedistrib_plot_tags(delta_susc_prop,delta_primary,plot_season_agedist,clin_fract_age_exp,
                                           preseas_npi_on,postseas_npi_off,npi_reduc_strength,seas_case_threshold=1e2)
ggsave(agedistr_filename,width=32,height=16,units="cm")

### ### ### ### ### 
### mean age of infections ----------------------------
# mean_age_perseason_exp_dep # =mean_age_perseason
mean_age_perseason=fcn_calc_mean_age(season_peaks_AUC,season_min=4)
# segment plot
var_sel="value"; if(grepl("norm",var_sel)) {ylab_tag=" (normalised)"} else {ylab_tag=" (year)"}
fcn_plot_mean_age(mean_age_perseason,seas_lims,npi_year,shutdwn_lims,yexpand=c(0.2,0),ylab_tag)
# SAVE
filenm=gsub("pct","pct_geomsegment",gsub("age_distrib","mean_age",agedistr_filename))
ggsave(filenm,width=22,height=14,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### UK RSV data  --------------------------------------------------------
resp_virus_data_uk=read_csv("data/Respiratory viral detections by any method UK Ages.csv")
resp_virus_data_uk_tidy = resp_virus_data_uk %>% pivot_longer(!c("Year","startweek","Age")) %>% 
  mutate(Age=factor(gsub(" Y","Y",resp_virus_data_uk_tidyAge),levels=unique(gsub(" Y","Y",resp_virus_data_uk_tidyAge))))
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
ggplot(subset(resp_virus_data_uk_tidy,name %in% "RSV" & grepl("indiv",type)),aes(x=year_week,y=value,group=Age)) + 
  geom_area(aes(fill=Age),position=position_stack(reverse=T)) + 
  # geom_line(aes(color=Age)) + geom_point(size=0.4) + facet_wrap(~Age,ncol=2,scales="free") + # 
  xlab("Year, week") + ylab("number of reported RSV cases (4-week period)") + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(size=12,vjust=0.5),axis.text.y=element_text(size=13),legend.position="bottom")
#  labs(caption="source: gov.uk/health-and-social-care/health-protection-infectious-diseases")
# ggsave("simul_output/uk_rsv_data2014_2020.png",width=32,height=16,units="cm")
# ggsave("simul_output/uk_rsv_data2014_2020_age_facet.png",width=32,height=16,units="cm")
# plot MEAN vs 2020
ggplot(resp_virus_data_uk_tidy[truthvals_rsv,],aes(x=startweek,y=value,group=Year,color=Year)) + # ,linetype=factor(type)
  geom_line(size=1.25) + geom_point(aes(shape=factor(type),size=factor(width))) + facet_wrap(~Age,scales="free") + 
  scale_x_continuous(breaks=(0:10)*5) + scale_linetype_manual(values=c("solid","dashed")) + theme_bw() + standard_theme + ylab("") +
  scale_size_manual(values = c(1.5,2),guide=FALSE) + labs(shape="data type") # ,linetype="data type"
# SAVE
ggsave("simul_output/uk_rsv_data2020_multiyearaverage.png",width=32,height=16,units="cm")

### data with weekly resolution (no age resol) ----
x=read_csv("data/Respiratory viral detections by any method UK.csv")
resp_detects_weekly_all_age= x %>% mutate(year_week=factor(paste0(Year,"-",Week),unique(paste0(Year,"-",Week))), 
    RSV_rolling_av=rollmean(RSV,k=7,align="center",fill=NA),section=ceiling((as.numeric(year_week)*0.92)/100) ) %>% 
  mutate(section=ifelse(section>3,3,section)) # year_week=factor(year_week,unique(year_week)),
# plot
ggplot(resp_detects_weekly_all_age,aes(x=year_week)) + geom_point(aes(y=RSV,group=1,color=factor(Year))) + 
  geom_line(aes(y=RSV_rolling_av,group=1)) + labs(color="Year") + scale_y_continuous(expand = expansion(0,0.05)) + 
  geom_vline(data=subset(resp_detects_weekly_all_age, Week %in% c(40,11)),aes(xintercept=year_week),color="red",linetype="dashed",size=0.5) + 
  facet_wrap(~section,nrow=3,scale="free_x") + theme_bw() + standard_theme + theme(axis.text.x = element_text(vjust=0.5,size=6))
# save
ggsave("simul_output/uk_rsv_data2020_allagegroups_weekly.png",width=32,height=20,units="cm")