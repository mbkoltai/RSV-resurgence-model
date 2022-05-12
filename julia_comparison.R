# Run individual simulations
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
rm(list=ls())
# setwd("Desktop/research/models/RSV_model/transmission_model/RSV-resurgence-model/")
if (!any(row.names(installed.packages()) %in% "here")) {install.packages("here")}; library(here)
here::i_am()
# currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# load libraries
x1<-c("tidyverse","deSolve","lubridate","gtools","ISOweek")
# "rstudioapi","devtools","wpp2019","Rcpp","tsibble","pracma","qs","ungeviz","zoo","RcppRoll"
# if (!any(grepl("ungeviz",row.names(installed.packages())))) { devtools::install_github("wilkelab/ungeviz") }
x2 <- x1 %in% row.names(installed.packages()); if (any(x2 == FALSE)) { install.packages(x1[! x2]) }
# Load all packages    
lapply(x1,library,character.only=TRUE) # as.Date <- zoo::as.Date
# select <- dplyr::select; # row_number <- dplyr::row_number; summarise <- dplyr::summarise
source('fcns/essential_fcns.R')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# set plotting theme (remove Calibri if it throws an error)
standard_theme=theme(plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90,vjust=1/2),
                     axis.text.y=element_text(size=9),axis.title=element_text(size=14), panel.grid.minor=element_blank()) 
# panel.grid=element_line(linetype="solid",colour="black",size=0.1) # , text=element_text(family="Calibri")
every_nth = function(n) { return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]}) }
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SET PARAMETERS --------------------------------------------------------
source("load_params.R")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# set variable/control parameters
# SUSCEPTIBILITY (normalised by age group sizes, for infection terms ~ delta*(I1+I2+...+In)*S_i/N_i)
# agedep_fact determines strength of age dependence, if agedep_fact>1, decreasing susceptibility with age
exp_dep <- 1.750; age_dep <- 0.3125 # partable$exp_dep[k_par]
const_delta <- 2.55; delta_primary <- const_delta*exp(-exp_dep*(1:3)) # partable$const_delta[k_par]
# matrix of susceptibility parameters
delta_susc <- round(sapply(1:n_age, function(x) {delta_primary/(exp(age_dep*x))}),5)

# calculate R0 *at the baseline* (no forcing)
R0_calc_SIRS(C_m,delta_susc,rho,n_inf)

### ### ### ### ### ### ### ### ### ### ### ### 
# DURATION of SIMULATION
# seasonal forcing (baseline level=1, forcing_strength=2 means 200% above baseline) 
# npi_reduc_strength: reduction from baseline 
# set seas lims from UK data: peak is weeks 49/50, on/off is 41,11
npi_dates=as.Date(c("2020-03-26","2021-05-17"))
g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% fun_shutdown_seasforc(npi_dates,
                    years_pre_post_npi=c(9,3),season_width_wks=5,init_mt_day="01-01", # "06-01"
                    peak_week=44,forcing_above_baseline=1,npireduc_strength=0.9)
# plot seasonal forcing
fcn_plot_seas_forc(simul_start_end,forcing_vector_npi,seas_lims_wks=c(7,42),npi_dates,date_resol="3 month")
# interpolation fcns for seas forcing & extern introds
approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
# how many introductions every 30 days?
approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*5))

### ### ### ### ### ### ### ### ### ### ### ### 
# set INITIAL VALUES
# pick a timepoint
# initvals_sirs_model = round(matrix(ode_sol[1700,(1:(ncol(ode_sol)-1))+1])); initvals_sirs_model[100:132]=0
# init_set can be "from scratch" (all susceptible) OR last state of a previous simulation
initvals_sirs_model <- fcn_set_initconds(rsv_agegroup_sizes=rsv_age_groups$stationary_popul,
                                         init_set=c("previous","fromscratch")[1],init_cond_src=c("output","file")[1],
                                         input_from_prev_simul=ode_solution,init_seed=10,seed_vars="all",filename="")

# all parameters collected as a list (sent to fcn running ODE solver)
deaths_vector <- matrix(unlist(lapply(uk_death_rate,function(x) rep(x,n_inf*n_compartment))))
params <- list(list(birth_rates,deaths_vector),K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,delta_susc,mat_imm_inds)

# SOLVE ODE SYSTEM
# can use 'lsoda' (faster) or 'lsodes'
tm <- proc.time(); ode_sol <- lsoda(initvals_sirs_model,timesteps,
                                    func=sirs_seasonal_forc_mat_immun,parms=params); round(proc.time()-tm,2)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PROCESS OUTPUT
### ### ### ### ### ### ### ### ### ### ### ### 
# this step only needed if want to investigate other variables or want to use it as input for next simulation (initial conds)
# reshape data: get all variables in `ode_solution_tidy` 
g(final_pop,ode_solution,ode_solution_tidy) %=% fun_process_simul_output(ode_sol,varname_list,
                    incidvar="newinf",incid_only=F,init_date=simul_start_end[1],n_age,n_inf,rsv_age_groups,neg_thresh=-0.01)
# has population size converged? (percentage difference)
# c(round(100*abs(final_pop[,3]-final_pop[,2])/final_pop[,2],2))
### ### ### ### ### ### ### ### ### ### ### ### 
# Extract infection *incidence*
df_cases_infs <- fcn_process_odesol_incid(ode_sol,n_age,n_inf,n_compartment,simul_start_end) %>% mutate(value=round(value,2))

# plot *prevalence* of infections
ode_solution_tidy %>% filter(compartment %in% "I" & (agegroup %in% 1) & grepl("#1",infection) & t<=3e3) %>% #   & t<=3600
  ggplot(aes(x=t,y=value,color=infection)) + geom_line() + 
  scale_x_continuous(breaks=(0:50)*100,expand=expansion(0.01,0)) + scale_y_continuous(breaks=(0:30)*5e3) + 
  ylab("I(t)") + theme_bw() + standard_theme

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# compare to julia output
ode_sol_julia<-read_tsv("../2022_new_project/julia_SIRS/plots/model_output.csv",col_names=F)
colnames(ode_sol_julia) <- colnames(ode_sol)
ode_sol_julia <- ode_sol_julia[-(which(duplicated(ode_sol_julia$time))-1),]
g(final_pop_julia,ode_solution_julia,ode_solution_tidy_julia) %=% fun_process_simul_output(ode_sol_julia,varname_list,
            incidvar="newinf",incid_only=F,init_date=simul_start_end[1],n_age,n_inf,rsv_age_groups,neg_thresh=-0.01)
# plot together
bind_rows(
  ode_solution_tidy %>% filter(compartment %in% "I") %>% mutate(solver="R"), #  & (agegroup %in% 11) & grepl("#3",infection)
  ode_solution_tidy_julia %>% filter(compartment %in% "I") %>% #  & (agegroup %in% 11) & grepl("#3",infection)
    mutate(solver="julia")) %>% # ,t=t/1.033
  filter(t<=2400 & agegroup<=5) %>% # t>=1720 & 
  ggplot(aes(x=t,y=value,color=solver,linetype=solver,size=solver)) + geom_line() + # 
  facet_grid(infection~agegroup,scales = "free_y") + scale_size_manual(values=c(1,1.3)) +
  scale_x_continuous(expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.01,0)) + #breaks=(0:1e3)*50, #breaks=(0:30)*5e3,
  ylab("I(t)") + theme_bw() + standard_theme + theme(legend.text=element_text(size=15))
# save
ggsave(paste0("../2022_new_project/julia_SIRS/plots/R_julia_comparison_multi_age_rot90.png"),width=28,height=16,units="cm")
