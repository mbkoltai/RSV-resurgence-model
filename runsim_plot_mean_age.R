rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# library(contactdata); library(fitdistrplus);  library(bbmle); library(Rcpp); library(GillespieSSA)
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","plotly","Rcpp","zoo","lubridate","tsibble","pracma",
         "qs","ungeviz"),library,character.only=TRUE) # 
source('fcns/RSV_model_functions.R')
# selected country
country_sel="United Kingdom"
# time resolution (in days)
elem_time_step=0.5
# population data
standard_age_groups <- fun_cntr_agestr(country_sel,i_year="2020",seq(0,75,5),c(seq(4,74,5),99))
popul_struct=fcn_cntr_fullpop(n_year="2020",country_sel)
# RSV age groups (population data from wpp2019)
rsv_age_groups=fun_rsv_agegroups(standard_age_groups,popul_struct,rsv_age_groups_low=c(0,0.5,1,1.5, 2,3,4, 5,10,15, 20),
                                 rsv_age_group_sizes=c(rep(0.4,4),rep(0.9,3),rep(4,3),79))
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
# if UK -> England's contact matrix # check: cm_parameters_SEI3R(cm_uk_locations("UK", 1))$pop[[1]]$matrices 
list_contmatrs=fun_covidm_contactmatrix(country_sel,currentdir_path,cm_path=cm_path) 
# make matrix reciprocal
C_m_full=Reduce('+',list_contmatrs) # fun_recipr_contmatr(Reduce('+',list_contmatrs),age_group_sizes=standard_age_groups$values)
# create for our age groups
C_m_merged_nonrecipr=fun_create_red_C_m(C_m_full,rsv_age_groups,
                                        orig_age_groups_duration=standard_age_groups$duration,orig_age_groups_sizes=standard_age_groups$values)
# make it reciprocal for the larger group
C_m=fun_recipr_contmatr(C_m_merged_nonrecipr,age_group_sizes=rsv_age_groups$value)
# bc of reinfections we need to input contact matrix repeatedly
contmatr_rowvector=t(do.call(cbind, lapply(1:nrow(C_m), function(x){diag(C_m[x,]) %*% matrix(1,n_age,n_inf)})))
# build kinetic matrix
# WANING (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j
omega=1/350 # 1/runif(1,60,200)
# RECOVERY
rho=1/7 # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,rsv_age_groups)
# BIRTH RATE into S_1_1 (Germany 2019: 778e3 births)
birth_rates=matrix(c(713e3/365,rep(0,dim_sys-1)),dim_sys,1)
# DEATHS (2019: 530841 deaths [England and Wales!])
uk_death_rate=read_csv("data/uk_death_rate_byage.csv")
g(rsv_age_groups,death_rates) %=% fun_death_rates(rsv_age_groups,uk_death_rate,nage=n_age,ninf=n_inf,dimsys=dim_sys)
# variable parameters  --------------------------------------------------------
# SUSCEPTIBILITY (normalised by age group sizes, for infection terms ~ delta*(I1+I2+...+In)*S_i/N_i)
# agedep_fact determines strength of age dependence, if agedep_fact>1, decreasing susceptibility with age
agedep_fact=1.25; delta_primary=rep(0.1,3)
# c(0.27,0.03,0.01)/5 # c(0.21,0.11,0.01)/5 # c(0.09,0.07,0.05)/1.6 # rep(0.15,3) # rep(0.09,3)
delta_susc=sapply(1:n_age, function(x) {delta_primary/((agedep_fact^(x-1))*rsv_age_groups$value[x])})
delta_susc_prop=delta_susc*matrix(rep(rsv_age_groups$value,3),nrow=3,byrow=T)
### PLOT susceptibility:
# fcn_plot_suscept_table(fcn_suscept_agedeptable(rsv_age_groups,delta_susc,n_inf)) 
# ggsave(paste0("simul_output/suscept_age_dep/",agedep_fact,"delta_prim",unique(delta_primary),".png"),width=32,height=22,units="cm")
# calculate R0 (at max seasonal forcing=1)
R0_calc_SIRS(C_m,delta_susc_prop,rho,n_inf)
# DURATION of SIMULATION
n_years=11.25; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=elem_time_step)
# seasonal forcing (baseline level=1, forcing_strength=2 means 200% above baseline) | npi_reduc_strength: reduction from baseline 
# shutdown season (if x --> (x+1)th season shut down) | preseas_npi_on/postseas_npi_off: on/off NPI week before/after season onset
forcing_strength=0.5; npi_reduc_strength=0.75; npi_year=round(n_years-4); preseas_npi_on=-2; postseas_npi_off=4
g(forcing_vector_npi,shutdwn_lims,seas_force,seas_lims) %=% fun_shutdown_seasforc(timesteps,seas_lims_real=c(37,9),
      forcing_strength,npi_reduc_strength,npi_year,peak_week=48,season_width=3,preseas_npi_on,postseas_npi_off,n_prec=0.01,n_sd=2)
# fcn_plot_seas_forc(timesteps,seas_force,forcing_vector_npi,shutdwn_lims,seas_lims)
# cosine fcn
# cosfunc <- function(beta, eta, phi, k_exp, time, k_timescale) { beta*(1+eta*cos(k_timescale*pi*(time-phi)/365)^k_exp) }
# l_time=3*365*2; data.frame(time=timesteps[1:(l_time+1)],
#       cos_fcn=cosfunc(time=(0:l_time)/2,beta=1,eta=1,phi=50,k_exp=2,k_timescale=2.3), norm_fcn=forcing_vector_npi[1:(l_time+1)]) %>%
#   pivot_longer(!time) %>% ggplot(aes(x=time,y=value,color=name)) + geom_line() + theme_bw() + standard_theme + 
#   scale_x_continuous(breaks=(0:40)*30)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# INITIAL CONDITIONS. # introduce stationary state as init state? | init_set: "previous" or anything else
initvals_sirs_model=fcn_set_initconds(init_set=c("previous","fromscratch")[1],init_cond_src=c("output","file")[2],ode_solution,
                                      init_seed=10,seed_vars="all",filename="simul_output/df_ode_solution_UK_long.RDS")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### integrate ODE --------------------------------------------------------
params=list(list(birth_rates,death_rates),K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,elem_time_step,delta_susc)
# interpolation fcns for seas forcing & extern introds
approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
# how many introductions every 30 days?
approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*100))
tm<-proc.time(); ode_solution<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params); round(proc.time()-tm,2)
# reshape data | # check size of objs: fcn_objs_mem_use(1)
g(ode_solution,df_cases_infs) %=% fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups,neg_thresh=-0.01)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Plot time course --------------------------------------------------------
xval_lims=c(npi_year-2.3,npi_year+4.9); xval_breaks=seq(0,max_time/365,by=1/4)
# PLOT
# tags: y-axis fixed/free | facet by AGE/(age&INFECTION) | abs values/fraction
# k=7 --> (free scale, infects separate, absval). k=8 --> (free, infects separate, fraction)
g(scale_val,facet2tag,value_type,y_axis_tag,ncol_val,facet_formula,foldername,caption_txt,subtitle_str,timecourse_filename) %=% 
  fun_tcourse_plottags(k=8,nval=2,rval=3,n_inf,n_age,colvar="age",agegr_lim=7,delta_susc_prop,delta_primary,agedep_fact,
                       preseas_npi_on,postseas_npi_off,npi_reduc_strength,forcing_strength)
# PLOT time course absolute/fract values by age group  --------------------------------------------------------
fcn_plot_allcases_absval_stackedinfs(df_cases_infs,value_type,x_lims=xval_lims,t_subset=7,agegrlim=6,ncolval=3,y_axis_tag,
          scale_val,vertline_x=subset(seas_lims,on>xval_lims[1]&off<xval_lims[2]),shutdwn_lims,xval_breaks,subtitle_str,caption_txt)
## SAVE
ggsave(timecourse_filename, width=32,height=20,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### mean age per season of 1st/2nd/3rd infection --------------------------------------------------------
season_means_byinf=subset(fcn_calc_mean_age_byinfs(rsv_age_groups,df_cases_infs,seas_lims,low_thresh=1e3,n_aver=30), 
                          grepl("in",on_off)) %>% mutate(ageweight_aver=mean_age_at_t*suminfs) %>% 
  group_by(epi_year,infection) %>% summarise(weighted_sum=ifelse(sum(suminfs)>5e4,sum(ageweight_aver)/sum(suminfs),NA))
# plot means per season by infection type
plot_xlims=c(npi_year-2,npi_year+3)
ggplot(subset(season_means_byinf %>% mutate(weighted_sum=ifelse(epi_year==npi_year+1,NA,weighted_sum)),
              epi_year>plot_xlims[1]&epi_year<=plot_xlims[2])) + 
  geom_hpline(aes(x=epi_year,y=weighted_sum,color=ifelse(epi_year==npi_year+2,"red",NA)),width=0.9,show.legend=F) +
  facet_wrap(~infection,scales="free") + theme_bw() + standard_theme + theme(panel.grid.major.x=element_blank()) +
  geom_rect(aes(xmin=npi_year+0.5,xmax=npi_year+1.5,ymin=-Inf,ymax=Inf),fill="pink",alpha=0.2) + xlab("RSV season")+ylab("mean age") +
  geom_vline(xintercept=(plot_xlims[1]:plot_xlims[2])+0.5,size=0.2,linetype="dashed",color="black") + 
  theme(axis.text.x=element_text(vjust=0.5)) + labs(subtitle=subtitle_str,caption=caption_txt)
# save
ggsave(paste0(paste0(strsplit(timecourse_filename,"/")[[1]][1:2],collapse = "/"),
        "/meanage_perseas_byinf",ifelse(agedep_fact>1,paste0("_agedep_",agedep_fact,"_delta",unique(delta_primary)),
                                        paste0("_expdep_",paste0(delta_primary,collapse="_")) ),
        "_peakforc",(1+forcing_strength)*1e2,"_NPIred",npi_reduc_strength*1e2,".png"), width=32,height=18,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### CALCULATE SYMPTOMATIC CASES & sum of 1,2,3rd infections -------------------------------
# dependence of SEVERITY on age (>1) | dependence on prev exposure
sever_agedep_fact=1; clinfract_expos=rep(1,n_inf) # c(0.8,0.4,0.2) # rep(0.75,3)
clin_fract_age_exp=fcn_clin_fract_age_exp(sever_agedep_fact,clinfract_expos,rsv_age_groups,n_age,n_inf)
# plot: fcn_plot_clinfract_table(clin_fract_age_exp)
df_sympt_cases=fun_symptomcases_table(df_cases_infs,clin_fract_age_exp,c("infection","agegroup"))

# PLOT time course SYMPTOMATIC CASES -----------------------------------
# k_plot: free/fixed | fraction/cases
g(scale_val,y_axis_tag,plotvar,foldername,symptcases_filename,caption_txt,subtitle_str) %=% fun_sumcase_plot_tags(n_val=2,r_val=2,
      k_plot=1,clin_fract_age_exp,delta_susc_prop,delta_primary,preseas_npi_on,postseas_npi_off,npi_reduc_strength,forcing_strength)
xval_lims=c(npi_year-2.5,npi_year+3.5); seas_lims_plot=subset(seas_lims,on>xval_lims[1]&off<xval_lims[2]) %>% pivot_longer(col=!season)
# plot
fcn_plot_symptomcases_agefacet(df_sympt_cases %>% filter(agegroup<=6),xval_lims,y_plotvar=plotvar,y_axis_tag,scale_val,seas_lims_plot,
                               shutdwn_lims,xval_breaks,n_col=3,n_agelim=11,subtitle_str,caption_txt)
# SAVE
# full_filename=gsub('.png','_reducedresurg.png',full_filename)
ggsave(symptcases_filename,width=31,height=22,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### mean age per season ---
low_thr=1e3
mean_age_symptcases_at_t <- left_join(rsv_age_groups %>% mutate(agegroup=as.numeric(factor(agegroup_name,levels=agegroup_name))) %>% 
      select(agegroup,mean_age_weighted),subset(df_sympt_cases,compartment %in% "I"),by="agegroup") %>% group_by(t) %>% # t,
  summarise(suminfs=sum(symptom_cases),mean_age_at_t=sum(symptom_cases*mean_age_weighted/sum(symptom_cases))) %>% 
  mutate(mean_age_at_t_smooth=ifelse(suminfs>low_thr,rollmean(mean_age_at_t,k=30,align="center",fill=NA),NA)) %>%
  mutate(on_off=ifelse(mod(findInterval(t/365,array(matrix(t(seas_lims[,c("on","off")])))),2)==1 & t/365>=seas_lims$on[1],
                       "in-season","off-season"),epi_year=findInterval(t/365,seas_lims$on)) %>% group_by(epi_year,on_off) %>% # ,infection
  mutate(season_mean=ifelse(mean(suminfs)>low_thr,mean(mean_age_at_t),NA))
# mean age per season
season_means_symptcases=subset(mean_age_symptcases_at_t, grepl("in",on_off)) %>% mutate(ageweight_aver=mean_age_at_t*suminfs) %>% 
  group_by(epi_year) %>% summarise(weighted_sum=ifelse(sum(suminfs)>5e4,sum(ageweight_aver)/sum(suminfs),NA))
# PLOT
plot_xlims=c(npi_year-2,npi_year+4)
ggplot(subset(season_means_symptcases %>% mutate(weighted_sum=ifelse(epi_year==npi_year+1,NA,weighted_sum)),
              epi_year>plot_xlims[1]&epi_year<=plot_xlims[2])) +
  geom_hpline(aes(x=epi_year,y=weighted_sum,color=ifelse(epi_year==npi_year+2,"red",NA)),width=0.9,show.legend=F) +
  theme_bw() + standard_theme + theme(panel.grid.major.x=element_blank()) + scale_x_continuous(breaks=1:round(n_years)) +
  geom_rect(aes(xmin=npi_year+0.5,xmax=npi_year+1.5,ymin=-Inf,ymax=Inf),fill="pink",alpha=0.2) + xlab("RSV season")+ylab("mean age") +
  geom_vline(xintercept=(plot_xlims[1]:plot_xlims[2])+0.5,size=0.2,linetype="dashed",color="black") + 
  theme(axis.text.x=element_text(vjust=0.5,size=13),axis.text.y=element_text(size=13),axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15)) + labs(subtitle=subtitle_str,caption=caption_txt)
# SAVE
ggsave(paste0(paste0(unlist(strsplit(symptcases_filename,"/"))[1:3],collapse="/"),"/mean_age_symptcases_peak_forc",(1+forcing_strength)*1e2,
    "_NPIred",npi_reduc_strength*1e2,"pct",ifelse(agedep_fact>1,paste0("_agedep",agedep_fact),""),"_delta",
    paste0(unique(delta_primary),collapse="_"),".png"),width=22,height=14,units="cm")
