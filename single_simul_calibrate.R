rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","plotly","Rcpp","zoo","lubridate",
         "tsibble","pracma","qs","ungeviz"),library,character.only=TRUE)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SET PARAMETERS --------------------------------------------------------
source("load_params.R")
# contact matrix from covidm ("home","work","school","other")
cm_path="~/Desktop/research/models/covid_model/lmic_model/covidm/"
# if UK -> England's contact matrix # check: cm_parameters_SEI3R(cm_uk_locations("UK", 1))$pop[[1]]$matrices 
list_contmatrs=fun_covidm_contactmatrix(country_sel,currentdir_path,cm_path=cm_path) 
#
# estimated attack rates
n_error_tol<-2.5; estim_attack_rates <- data.frame(agegroup_name=paste0("age=",rsv_age_groups$agegroup_name,"yr"),
      median_est=c(rep(65,4),rep(40,4),10,8,5)) %>% mutate(min_est=median_est/n_error_tol,max_est=median_est*n_error_tol)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# duration, npi timing 
npi_dates=as.Date(c("2020-03-26","2021-04-01")); seaspeakval=1/2; seasforc_width_wks=8
simul_length_yr=10; post_npi_yr=4
l_seas<-fun_shutdown_seasforc(npi_dates,years_pre_post_npi=c(simul_length_yr-post_npi_yr,post_npi_yr),
  season_width_wks=seasforc_width_wks,init_mt_day="06-01",peak_week=c(45,49,51)[2],
  forcing_above_baseline=seaspeakval, npireduc_strength=0.5); g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% l_seas
approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*10))
# deaths/popul (by age group)
death_rates <- matrix(unlist(lapply(c(rep(1e-5,2)*3,rep(1e-6,5),rep(0,2),1e-6,1.79e-4),
                                    function(x) rep(x,n_inf*n_compartment))))
# maternal immunity
mat_imm_flag=TRUE
mat_imm_inds<-list(fun_sub2ind(i_inf=1,j_age=1,"R",c("S","I","R"),11,3),fun_sub2ind(i_inf=c(1,2,3),j_age=9,"R",
    c("S","I","R"),11,3),fun_sub2ind(i_inf=c(1,2,3),j_age=9,"S",c("S","I","R"),11,3))
# susceptibility
const_delta <- 1/12; exp_dep<-3/4; age_dep<-1/5; 
delta_primary<-const_delta*exp(-exp_dep*(1:3)); delta_susc<-sapply(1:n_age, function(x) {delta_primary/(exp(age_dep*x))})
# R0_calc_SIRS(C_m,delta_susc,rho,n_inf)
###
# DEFINE parameters
if (mat_imm_flag){ params<-list(list(birth_rates,death_rates),K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,delta_susc,mat_imm_inds)
} else { params<-list(list(birth_rates,death_rates),K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,delta_susc) }
# initial conditions
initvals_sirs_model <- fcn_set_initconds(rsv_age_groups$stationary_popul,
    init_set=c("previous","fromscratch")[2],init_cond_src=c("output","file")[1],
    ode_solution[1:(ncol(ode_solution)-1),],init_seed=10,seed_vars="all",filename="")
# SIMULATE
tm<-proc.time();ode_solution<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc_mat_immun,parms=params);(proc.time()-tm)[3]
df_cases_infs <- fcn_process_odesol_incid(ode_solution,n_age,n_inf,n_compartment,simul_start_end)

# simple lineplot
sel_weeks <- df_cases_infs %>% mutate(week=week(date),year=year(date)) %>% filter(week %in% c(42,52,10)) %>%
 group_by(year,agegroup,week) %>% filter(date==min(date)&infection==1); sel_yrs<-as.Date(c("2017-10-01","2022-03-01"))
p <- fcn_plot_timecourse_by_agegr(df_cases_infs %>% filter(t %% 7==0 & agegroup<=9 & date>sel_yrs[1] & date<sel_yrs[2]),
  agegroup_name=rsv_age_groups$agegroup_name,sel_agelim=9,varname="value",npidates=npi_dates,date_break_val="2 month",
  selweeks=sel_weeks,alphaval=0.01,vline_w=c(1/4,1/8)); t_str<-paste0("exp_dep=",exp_dep,", age_dep=",
    age_dep," (susc. const.=",round(const_delta,2),")")
p+theme(legend.position="top")+labs(fill="# infection")+ylab("")+ggtitle(t_str)

# plot
# calculate attack rates by epi-year and age group # 81-97% of cases should be within weeks 41 and 9
sum_inf_epiyear_age <- left_join(df_cases_infs %>% mutate(year=year(date),
    agegroup_name=rsv_age_groups$agegroup_name[agegroup],epi_year=ifelse(date>ymd(paste(year(date),"-07-01")),year(date),year(date)-1),
    in_out_season=ifelse(week(date)<=9 | week(date)>=41,"in","out")) %>% group_by(epi_year,agegroup_name) %>%
    summarise(inf_tot=sum(value),inf_in_seas=sum(value[in_out_season=="in"])) %>% group_by(agegroup_name) %>% 
    filter(epi_year>min(epi_year)),rsv_age_groups,by="agegroup_name") %>% 
    mutate(attack_rate=100*inf_tot/stationary_popul,seas_share=inf_in_seas/inf_tot,
         agegroup_name=factor(agegroup_name,levels = unique(rsv_age_groups$agegroup_name)))
# plot attack rates
ggplot(sum_inf_epiyear_age %>% mutate(attack_rate=ifelse(epi_year==2020,NA,attack_rate)) %>%filter(epi_year>2016)) + 
  geom_hpline(aes(x=factor(epi_year),y=attack_rate),width=0.9,size=1) + facet_wrap(~agegroup_name,scales="free_y") + 
  theme_bw() + standard_theme + xlab("epi-year") + ylab("attack rate % age group") +
  geom_rect(xmin=(2020-2016)-0.5,xmax=(2020-2016)+1/2, ymin=-Inf,ymax=Inf,fill="pink",alpha=0.1) + 
  geom_hline(data=estim_attack_rates %>% pivot_longer(!agegroup_name) %>%
        mutate(agegroup_name=factor(gsub("age=|yr","",agegroup_name),levels=unique(rsv_age_groups$agegroup_name))),
      aes(yintercept=value,linetype=ifelse(name %in% "median_est","solid","dashed")),size=1/3,color="red",show.legend=F) +
  labs(caption=paste0("seas peak=",round(max(forcing_vector_npi),2),"x baseline, NPI contact red.=",
        round((1-min(forcing_vector_npi))*1e2),"%")) + ggtitle(paste0("expos_dep=",exp_dep,", age_dep=",age_dep))

# attack rates
# View(df_cases_infs %>% filter(date>as.Date("2018-10-15")&date<as.Date("2019-03-01")) %>% group_by(agegroup,infection) %>%
#   summarise(cumul_inf=round(sum(value))) %>% mutate(agegr_prop=round(100*cumul_inf/rsv_age_groups$value[agegroup],1)))

# write_csv(data.frame(initvals_sirs_model),
#           paste0("simul_output/stat_sol_expon_dep_agedep",age_dep,"_expos",exp_dep,".csv"))
# initvals_sirs_model <- as.matrix(
#   read_csv(paste0("simul_output/stat_sol_expon_dep_agedep",age_dep,"_expos",exp_dep,".csv")))
