# Run individual simulations
#################################
rm(list=ls())
if (!any(row.names(installed.packages()) %in% "here")) {install.packages("here")}; library(here)
# load constant parameters and functions for simulations, specify folder where inputs are stored
source(here::here("load_params.R"))
foldername <- "repo_data/"

# WANING (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j
omega=1/350 # 1/runif(1,60,200)
# RECOVERY
rho=1/7 # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX: this contains linear terms ie. recovery, aging and waning
# (aging terms need to be scaled by duration of age groups!)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,rsv_age_groups)

birth_rates=matrix(c(713e3/365,rep(0,dim_sys-1)),dim_sys,1)
# DEATHS (2019: 530841 deaths [England and Wales!]) # "uk_death_rate_byage_rsv_agegroups.csv" 
# is for 1000 population!

# SUSCEPTIBILITY (normalised by age group sizes, for infection terms ~ delta*(I1+I2+...+In)*S_i/N_i)
# agedep_fact determines strength of age dependence, if agedep_fact>1, decreasing susceptibility with age
exp_dep <- 1.750; age_dep <- 0.3125 # partable$exp_dep[k_par]
const_delta <- 2.549890; delta_primary <- const_delta*exp(-exp_dep*(1:3)) # partable$const_delta[k_par]
delta_susc <- sapply(1:n_age, function(x) {delta_primary/(exp(age_dep*x))})

# calculate R0
R0_calc_SIRS(C_m,delta_susc,rho,n_inf)

# DURATION of SIMULATION
# seasonal forcing (baseline level=1, forcing_strength=2 means 200% above baseline) 
# npi_reduc_strength: reduction from baseline 
# set seas lims from UK data: peak is weeks 49/50, on/off is 41,11
npi_dates=as.Date(c("2020-03-26","2021-05-17")); seaspeakval=1; seasforc_width_wks=5
g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% fun_shutdown_seasforc(npi_dates,
        years_pre_post_npi=c(15,3),
        season_width_wks=seasforc_width_wks,init_mt_day="06-01",
        peak_week=48,forcing_above_baseline=seaspeakval,npireduc_strength=0.5)
# plot seasonal forcing
fcn_plot_seas_forc(simul_start_end,forcing_vector_npi,seas_lims_wks=c(7,42),
                   npi_dates,date_resol="3 month")

# set initial condition (all susceptible)
initvals_sirs_model <- fcn_set_initconds(rsv_age_groups$stationary_popul,
      init_set=c("previous","fromscratch")[2],init_cond_src=c("output","file")[2],
      NA,init_seed=10,seed_vars="all",filename="")

# all parameters as argument to fcn that'll run ODE solver
params<-list(list(birth_rates,
                  matrix(unlist(lapply(uk_death_rate,function(x) rep(x,n_inf*n_compartment))))),
                  K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,delta_susc)
# interpolation fcns for seas forcing & extern introds
approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
# how many introductions every 30 days?
approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*5))
tm<-proc.time(); ode_sol<-lsoda(initvals_sirs_model,timesteps,
                                func=sirs_seasonal_forc,parms=params); round(proc.time()-tm,2)
# reshape data: get all variables in `ode_solution_tidy`
g(final_pop,ode_solution,ode_solution_tidy) %=% fun_process_simul_output(ode_sol,varname_list,
                              incidvar="newinf",incid_only=F,init_date=simul_start_end[1],n_age,n_inf,
                              rsv_age_groups,neg_thresh=-0.01)
# extract incident infections only
df_cases_infs <- fcn_process_odesol_incid(ode_sol,n_age,n_inf,n_compartment,simul_start_end) %>% 
  mutate(value=round(value,2))

# select weeks of season limits
sel_weeks <- df_cases_infs %>% mutate(week=week(date),year=year(date)) %>% 
  filter(week %in% c(9,42,49)) %>% group_by(year,agegroup,week) %>% 
  filter(date==min(date) & infection==1)

# PLOT incident infections, disaggregated by age AND 1/2/3rd infection
df_cases_infs %>% 
  filter(t %% 7==0 & agegroup<=9 & date>as.Date("2018-09-01") & date<as.Date("2022-04-01")) %>% 
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],
                              unique(rsv_age_groups$agegroup_name))) %>%
ggplot() + geom_area(aes(x=date,y=value,fill=factor(infection)),position=position_stack(reverse=T),
                alpha=0.3,color="black",size=0.25) + 
  facet_wrap(~agegroup_name,scales="free_y") + xlab("") + 
  geom_vline(data=sel_weeks %>% filter(agegroup<=9),
         aes(xintercept=date,linetype=ifelse(week==unique(sel_weeks$week)[2],"solid","dashed"),
           size=ifelse(week==49,"a","b")),show.legend=F) + 
  scale_size_manual(breaks=c("a","b"),values=c(1/4,1/8)) +
  scale_x_date(date_breaks="2 month",expand=expansion(0.01,0)) + 
  scale_y_continuous(expand=expansion(0.01,0)) +
  geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="pink",alpha=1/100) +
  theme_bw() + standard_theme + ylab(ifelse(varname=="value","# cases","% age group"))

# disaggregate by age, summed over infection type (1/2/3)
df_cases_infs %>% 
  filter(agegroup<=9 & date>as.Date("2018-09-01") & date<as.Date("2022-04-01")) %>% # t %% 7==0 & 
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],
  unique(rsv_age_groups$agegroup_name))) %>% 
  group_by(t,date,agegroup_name) %>% summarise(value=sum(value)) %>%
  ggplot() + geom_line(aes(x=date,y=value)) + facet_wrap(~agegroup_name,scales="free_y") + 
  geom_vline(data=sel_weeks %>% filter(agegroup<=9),
             aes(xintercept=date,linetype=ifelse(week==unique(sel_weeks$week)[2],"solid","dashed"),
                 size=ifelse(week==49,"a","b")),show.legend=F) + 
  scale_size_manual(breaks=c("a","b"),values=c(1/4,1/8)) +
  scale_x_date(date_breaks="2 month",expand=expansion(0.01,0)) + 
  scale_y_continuous(expand=expansion(0.01,0)) + 
  xlab("") + ylab(ifelse(varname=="value","# cases","% age group")) +
  geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="pink",alpha=1/100) +
  theme_bw() + standard_theme


# sum across age groups
df_cases_infs %>% filter(agegroup<=9 & date>as.Date("2018-09-01") & # t %% 7==0 & 
      date<as.Date("2023-04-01")) %>% group_by(date,infection) %>% summarise(value=sum(value)) %>%
  ggplot(aes(x=date,y=value,fill=factor(infection),group=infection)) + 
  geom_area(position=position_stack(reverse=T),alpha=0.3,color="black",size=0.25) + 
  scale_x_date(date_breaks="3 months",expand=expansion(0.01,0)) + 
  scale_y_continuous(expand=expansion(0.01,0)) + 
  geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="pink",alpha=0.01) + 
  geom_vline(data=sel_weeks %>% filter(week %in% c(42,9)),
             aes(xintercept=date,linetype="solid"),size=1/4,show.legend=F,linetype="dashed") + 
  ylab("# new infections") + theme_bw() + standard_theme
