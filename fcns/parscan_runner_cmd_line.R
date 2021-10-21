# CMD LINE input:
# nohup Rscript --vanilla parscan_runner_cmd_line.R 62 91 6 3 simul_output/parscan/parallel/initconds_all.csv
# load constant parameters and functions
source("load_params.R")
# parameter table
partable_file_name<-commandArgs(trailingOnly=TRUE)[5]
print(partable_file_name)
# estimated attack rates
estim_attack_rates <- read_csv(commandArgs(trailingOnly=TRUE)[6])
# % cases within season (filtering parameter sets)
partable <- read_csv(partable_file_name)
seas_conc_lim<-unique(partable$seas_conc_lim)
# save the stat sol of all param sets
stat_sol_allparsets=matrix(0,nrow=(n_compartment+1)*n_age*n_inf,ncol=nrow(partable))
# length of simulations
simul_length_yr<-as.numeric(commandArgs(trailingOnly=TRUE)[3])
post_npi_yr <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
# agegroup indices for maternal immunity
mat_imm_flag <- TRUE
mat_imm_inds<-list(fun_sub2ind(i_inf=1,j_age=1,"R",c("S","I","R"),11,3),fun_sub2ind(i_inf=c(1,2,3),j_age=9,"R",
              c("S","I","R"),11,3),fun_sub2ind(i_inf=c(1,2,3),j_age=9,"S",c("S","I","R"),11,3))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# LOOP
# nohup Rscript --vanilla parscan_runner_cmd_line.R 62 91 6 3 simul_output/parscan/parallel/initconds_all.csv
# tm<-proc.time()
k_start_end <- as.numeric(commandArgs(trailingOnly=TRUE))[1:2]; print(paste0(c("PARAM SETS:", k_start_end),collapse=" "))
partable <- partable[k_start_end[1]:k_start_end[2],] # %>% filter(R0_no!=1)
# init conds
init_cond_file_name <- commandArgs(trailingOnly=TRUE)[5]; init_conds <- read_csv(init_cond_file_name)
serial_loop=TRUE
# print(paste0("size of partable: ",nrow(partable)))
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
  # width of season (from peak)
  seasforc_width_wks<-partable$seasforc_width_wks[k_par]
  # NPI dates
  npi_dates=as.Date(c("2020-03-26","2021-05-17"))
  # set initial conds
  # initvals_sirs_model <- t(t(((init_conds %>% filter(dep_type==partable$dep_type[k_par] & R0==partable$R0[k_par]))[,1])))
  initvals_sirs_model <- fcn_set_initconds(rsv_age_groups$stationary_popul,
                            init_set=c("previous","fromscratch")[2],init_cond_src=c("output","file")[1],
                            NA,init_seed=10,seed_vars="all",filename="") # ode_solution[1:(ncol(ode_solution)-1),]
  print(paste0("PARSET: ",partable$par_id[k_par])); print("LOAD INITVALS")
  # set length of simulation and seasonality
  l_seas<-fun_shutdown_seasforc(npi_dates,
          years_pre_post_npi=c(simul_length_yr-post_npi_yr,post_npi_yr),
          season_width_wks=seasforc_width_wks,init_mt_day="06-01",partable$peak_week[k_par],
          forcing_above_baseline=partable$seasforce_peak[k_par], npireduc_strength=0.5)
  g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% l_seas
  # set params 
  params<-list(list(birth_rates,matrix(unlist(lapply(uk_death_rate,function(x) rep(x,n_inf*n_compartment))))),
               K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,delta_susc)
  # maternal immunity?
  if (mat_imm_flag){ params[[7]] <- mat_imm_inds}
  # interpolation fcns for seas forcing & extern introds
  approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
  approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*10))
  # run simul
if (!mat_imm_flag){ ode_solution <- lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params) } else {
    ode_solution <- lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc_mat_immun,parms=params)     }
  # process output
  df_cases_infs <- fcn_process_odesol_incid(ode_solution,n_age,n_inf,n_compartment,simul_start_end) %>% 
        mutate(par_id=partable$par_id[k_par]) %>% filter(date>as.Date("2018-09-01"))
  # print progress
  print(paste0(round(1e2*k_par/nrow(partable)),"% "))
  # , dep_val=",partable$dep_val[k_par],", R0=",round(partable$R0[k_par],2),", suscept=",paste0(round(delta_primary,3),
  # collapse=","),", seas peak rel. baseline=",(partable$seasforce_peak[k_par]+1)*100,"%") )
  sel_t_point <- unique((df_cases_infs %>% filter(date == as.Date("2019-07-01")))$t)
  stat_sol_allparsets[,k_par] <- matrix(ode_solution[sel_t_point,2:ncol(ode_solution)]) # nrow(ode_solution)-1
  # final population (it's stationary, shldnt change)
  final_pop <- data.frame(agegroup=1:n_age,final=round(unlist(lapply(lapply((0:(n_age-1))*(n_inf*n_compartment), 
    function(x) (1:(n_inf*n_compartment))+x ), function(x_sum) sum(stat_sol_allparsets[1:(n_age*n_inf*n_compartment),k_par][x_sum])))))
  # calc attack rates
  print(paste0("season start: ",partable$seas_start_wk[k_par]))
  sum_inf_epiyear_age <- left_join(df_cases_infs %>% mutate(year=year(date),
   epi_year=ifelse(date>ymd(paste(year(date),"-07-01")),year(date),year(date)-1),
   in_out_season=ifelse(week(date)<=partable$seas_start_wk[k_par]|week(date)>=partable$seas_stop_wk[k_par],"in","out")) %>%
     group_by(epi_year,agegroup) %>% summarise(inf_tot=round(sum(value,na.rm=T)),
                                               inf_in_seas=round(sum(value[in_out_season=="in"])),
    max_incid_week=mean(week(date[value==max(value,na.rm=T)]),na.rm=T)) %>% group_by(agegroup) %>% 
      filter(epi_year>min(epi_year)),final_pop,by="agegroup") %>% mutate(attack_rate_perc=100*inf_tot/final,
           seas_share=inf_in_seas/inf_tot,dep_val=partable$dep_val[k_par],par_id=partable$par_id[k_par],
           seasforce_peak=partable$seasforce_peak[k_par],dep_type=partable$dep_type[k_par],R0=partable$R0[k_par])
  # store parameters # list_delta_primary[[k_par]]=delta_primary
  # store outputs
  # if (k_par==1) {all_sum_inf_epiyear_age=sum_inf_epiyear_age; print("k_par=1")} else {
  #   print("k_par>1"); all_sum_inf_epiyear_age=bind_rows(all_sum_inf_epiyear_age,sum_inf_epiyear_age)
  all_sum_inf_epiyear_age=bind_rows(all_sum_inf_epiyear_age,sum_inf_epiyear_age)
    if (k_par==nrow(partable)) { 
      all_sum_inf_epiyear_age <- left_join(all_sum_inf_epiyear_age %>%
            mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup])),estim_attack_rates,by="agegroup_name") %>% 
        mutate(attack_rate_check=ifelse(attack_rate_perc>=min_est&attack_rate_perc<=max_est,T,F),
               seas_share_check=ifelse(seas_share>seas_conc_lim,T,F)) } # }
  df_cases_infs_all=bind_rows(df_cases_infs_all,df_cases_infs)
  # save intermed output
  # write_csv(sum_inf_epiyear_age,paste0("simul_output/parscan/parallel/sum_inf_epiyear_age",k_par,".csv"))
} # end loop

# summary of simuls
write_csv(all_sum_inf_epiyear_age,paste0("simul_output/parscan/parallel/summ_parsets_main",
                                         paste0(k_start_end[c(1,length(k_start_end))],collapse="_"),".csv"))
# dynamics
par_id_vals=partable$par_id
write_csv(df_cases_infs_all %>% select(!name),paste0("simul_output/parscan/parallel/dyn_parsets_main",
                                   paste0(k_start_end[c(1,length(k_start_end))],collapse="_"),".csv") )

