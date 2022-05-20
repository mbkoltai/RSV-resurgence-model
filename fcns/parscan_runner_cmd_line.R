# CMD LINE input:
# nohup Rscript --vanilla parscan_runner_cmd_line.R 62 91 6 3 simul_output/parscan/initconds_all.csv
# load constant parameters and functions
source("load_params.R")
### ###
print("# cmd line arguments:"); print(length(commandArgs(trailingOnly=TRUE)))
k_start_end <- as.numeric(commandArgs(trailingOnly=TRUE))[1:2]; 
print(paste0(c("PARAM SETS:", k_start_end),collapse=" "))
### ###
# length of simulations
simul_length_yr <- as.numeric(commandArgs(trailingOnly=TRUE)[3]); 
post_npi_yr <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
# reduction in contact levels during NPIs
contact_reduction <- as.numeric(commandArgs(trailingOnly=TRUE)[5]) # 0.9=90%
# parameter table
partable_file_name<-commandArgs(trailingOnly=TRUE)[6]; print(partable_file_name)
partable <- read_csv(partable_file_name)[k_start_end[1]:k_start_end[2],]
# estimated attack rates
# estim_attack_rates <- read_csv(commandArgs(trailingOnly=TRUE)[6])
save_flag <- ifelse(grepl("nosave|no_save|nosavedyn|NOSAVE",commandArgs(trailingOnly=TRUE)[7]),FALSE,TRUE)
print(paste0("SAVE DYNAMICS: ",save_flag))
# start date for saving dynamics
start_date_dyn_save<-commandArgs(trailingOnly=TRUE)[8]; print(paste0("SAVE output from: ",start_date_dyn_save))
agegroup_res<-commandArgs(trailingOnly=TRUE)[9]
# SEASON LIMITS: we fix these for given RSV seasonality
# seas_start_wk <- 42; seas_stop_wk<-9
# we define season in a 26-week window so that it's comparable to Kenya study
seas_start_wk <- 40; seas_stop_wk <- 13
# save the stat sol of all param sets
stat_sol_allparsets=matrix(0,nrow=(n_compartment+1)*n_age*n_inf,ncol=nrow(partable))
# when does an epi_year start? week 40 often used
# month_day_epiyear_start <- "07-01"
week_epiyear_start <- 40
### ### ### ### ### ### ### ###
# agegroup indices for maternal immunity
# mat_imm_flag <- TRUE
# mat_imm_inds<-list(fun_sub2ind(i_inf=1,j_age=1,"R",c("S","I","R"),11,3),
#                    fun_sub2ind(i_inf=c(1,2,3),j_age=9,"R",
#               c("S","I","R"),11,3),fun_sub2ind(i_inf=c(1,2,3),j_age=9,"S",c("S","I","R"),11,3))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# LOOP
# tm<-proc.time()
serial_loop=TRUE
summ_filename <- paste0("simul_output/parscan/summ_parsets_main",
                        paste0(k_start_end[c(1,length(k_start_end))],collapse="_"),".csv")
dyn_filename <- paste0("simul_output/parscan/dyn_parsets_main",
                       paste0(k_start_end[c(1,length(k_start_end))],collapse="_"),".csv")
# print("STARTING LOOP")
all_sum_inf_epiyear_age <- data.frame(); df_cases_infs_all <- data.frame()
for (k_par in 1:nrow(partable)){ # nrow(partable)
  ### ### ###
  # first R0 samples already simulated, so skip those
  # print(paste0("R0_no CURRENT VALUE: ",partable$R0_no[k_par]))
  # assign SUSCEPT params
  exp_dep <- partable$exp_dep[k_par]; age_dep <- partable$age_dep[k_par]
  const_delta <- partable$const_delta[k_par]; delta_primary <- const_delta*exp(-exp_dep*(1:3)) # 1.24
  delta_susc <- sapply(1:n_age, function(x) {delta_primary/(exp(age_dep*x))})
  peak_week <- partable$peak_week[k_par]
  # width of season (from peak)
  seasforc_width_wks <- partable$seasforc_width_wks[k_par]
  # NPI dates
  npi_dates=as.Date(c("2020-03-26","2021-05-17"))
  # set initial conds
  # initvals_sirs_model <- fcn_set_initconds(rsv_age_groups$stationary_popul,init_set=c("previous","fromscratch")[2],
  #         init_cond_src=c("output","file")[1],NA,init_seed=10,seed_vars="all",filename="") 
  # read in a solution with stabilised age structure
  initvals_sirs_model<-readRDS("repo_data/stationary_sol.RDS")
  
  print(format(Sys.time(),"%Y/%b/%d | %X"))
  print(paste0("PARSET: ",partable$par_id[k_par], "(",k_par,")")) # print("LOAD INITVALS")
  # set length of simulation and seasonality
  l_seas <- fun_shutdown_seasforc(npi_dates,years_pre_post_npi=c(simul_length_yr-post_npi_yr,post_npi_yr),
          season_width_wks=seasforc_width_wks,init_mt_day="01-01",
          peak_week=peak_week,
          forcing_above_baseline=partable$seasforce_peak[k_par], npireduc_strength=contact_reduction)
  g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% l_seas
  # if waning is a variable parameter
  if (any(colnames(partable) %in% "omega")){
    omega=partable$omega[k_par]
    # matrix of linear terms (recovery and aging btwn compartments)
    K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,
                            agegroup_durations=rsv_age_groups$duration)
    }
  # set params 
  deaths_vector <- matrix(unlist(lapply(uk_death_rate,function(x) rep(x,n_inf*n_compartment))))
  params <- list(list(birth_rates,deaths_vector), K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,delta_susc)
  # maternal immunity? (always on)
  if (mat_imm_flag){ params[[7]] <- mat_imm_inds}
  # interpolation fcns for seas forcing & extern introds
  approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
  approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*10))
  # run simul
if (!mat_imm_flag){ 
  ode_solution <- lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params) 
    } else {
      # we always run model with maternal immunity
      ode_solution <- lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc_mat_immun,parms=params)
    }
  # process output
  df_cases_infs <- fcn_process_odesol_incid(ode_solution,n_age,n_inf,n_compartment,simul_start_end) %>% 
                          mutate(par_id=partable$par_id[k_par]) %>% # ,value=value
                          filter(date>=as.Date(start_date_dyn_save)) %>% ungroup() %>% select(!name)
  
  sel_t_point <- unique((df_cases_infs %>% filter(date == as.Date("2019-07-01")))$t)
  # print progress
  print(paste0(round(1e2*k_par/nrow(partable),1),"% "))
  stat_sol_allparsets[,k_par] <- matrix(ode_solution[sel_t_point,2:ncol(ode_solution)]) # nrow(ode_solution)-1
  # final population (it's stationary, should not change)
  final_pop <- data.frame(agegroup=1:n_age,final=round(unlist(lapply(lapply((0:(n_age-1))*(n_inf*n_compartment), 
    function(x) (1:(n_inf*n_compartment))+x ), 
    function(x_sum) sum(stat_sol_allparsets[1:(n_age*n_inf*n_compartment),k_par][x_sum])))))
  # calc attack rates
  # print(paste0("season start: ",partable$seas_start_wk[k_par]))
  
  # summary stats before broader age group aggregation
  sum_inf_epiyear_age <- left_join(df_cases_infs %>%
                          mutate(year=year(date),
                                 epi_year=ifelse(isoweek(date)>=week_epiyear_start, year(date), year(date)-1),
                                 in_out_season=ifelse(isoweek(date) >= seas_start_wk | isoweek(date) <= seas_stop_wk,"in","out")) %>% 
                          group_by(epi_year,agegroup) %>%
                                     # summing 1st, 2nd, 3rd infections
              summarise(inf_tot=sum(value,na.rm=T), # round()
                        inf_in_seas=sum(value[in_out_season=="in"]), # round()
              # for the 1st age group we only count 1st infections for the attack rate, to be comparable to Kenya study
              inf_in_seas_AR=sum(value[in_out_season=="in" & !(agegroup==1 & infection>1)]), # round()
              peak_inf=round(max(value,na.rm=T)),
              max_incid_week=mean(isoweek(date[value==max(value,na.rm=T)]),na.rm=T)) %>%
              group_by(agegroup) %>% filter(epi_year>min(epi_year)),
              final_pop,by="agegroup") %>%  # LEFT JOIN END
    mutate(par_id=partable$par_id[k_par],
           exp_dep=partable$exp_dep[k_par],age_dep=partable$age_dep[k_par],
           seasforc_width_wks=seasforc_width_wks,
           seasforce_peak=partable$seasforce_peak[k_par],
           R0=partable$R0[k_par],
           omega=partable$omega[k_par],
           attack_rate_perc=100*inf_in_seas_AR/final, # calculate in-season attack rate! # no rounding round(,1)
           seas_share=inf_in_seas/inf_tot) %>% # round(,3)
    relocate(c(inf_tot,inf_in_seas,peak_inf,max_incid_week,attack_rate_perc,seas_share),.after=omega) %>%
    relocate(par_id,.before=epi_year)
  
  # aggregate into broader age groups
  if (grepl("broad|BROAD",agegroup_res)) {
    df_cases_infs <- df_cases_infs %>% 
      mutate(agegroup_broad=findInterval(agegroup,c(2,4,7,10)+1)+1) %>% 
      filter(agegroup!=4) %>% # remove age groups 5-65y
      group_by(t,agegroup_broad,par_id,date) %>% 
      summarise(value=sum(value)) %>% rename(agegroup=agegroup_broad)
    
    final_pop <- final_pop %>% 
      mutate(agegroup_broad=findInterval(agegroup,c(2,4,7,10)+1)+1) %>% 
      filter(agegroup!=4) %>% # remove age groups 5-65y
      group_by(agegroup_broad) %>% 
      summarise(final=sum(final)) %>% rename(agegroup=agegroup_broad) 
  }
  
  
  # SAVE summary stats
  write_csv(sum_inf_epiyear_age,summ_filename,append=ifelse(k_par>1,TRUE,FALSE))
  # SAVE dynamics
  if (save_flag) { write_csv(df_cases_infs %>% select(!date),dyn_filename,append=ifelse(k_par>1,TRUE,FALSE)) }
} # end loop
