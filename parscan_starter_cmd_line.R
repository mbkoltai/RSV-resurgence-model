# load constant parameters and functions
source("load_params.R")
# cmd line inputs
arg_vals_section_length_simuldur <- as.numeric(commandArgs(trailingOnly=TRUE)); print(arg_vals_section_length_simuldur)
# estimated attack rates
estim_attack_rates <- data.frame(agegroup_name=rsv_age_groups$agegroup_name, # paste0("age=",,"yr")
      median_est=c(rep(65,4),rep(40,4),10,8,5)) %>% mutate(min_est=median_est*0.25,max_est=median_est*2.5,
      median_all_inf=c(rep(70,4),rep(60,4),50,30,20),min_est_all_inf=median_all_inf*0.5,max_est_all_inf=median_all_inf*1.5)
# % cases within season (filtering parameter sets)
seas_conc_lim=0.8
# parameter sets to search through
selsets<-c(2,4:8) 
p_table <- bind_rows(expand.grid(list(dep_type="age",dep_val=seq(0.5,4,by=0.5)[selsets],R0=seq(12,14,0.5)/10,
    seasforce_peak=seq(1.1,1.5,by=0.1))),bind_rows(lapply(selsets, function(x) expand.grid(list(dep_type="exp",dep_val=x,
      seasforce_peak=list(c(0.8,1,1.2),c(1.125,1.25,1.375),c(1,1.25,1.5),c(1,1.25,1.5),
      c(1,1.125,1.25),c(1.25,1.375,1.5),c(1.25,1.375,1.5),c(1.25,1.375,1.5))[[x]],
R0=list( (12:14)/10,(12:14)/10,seq(13,16,1.5)/10,seq(15,18,1.5)/10, seq(14,20,2)/10, 
         seq(15,18,1)/10,(18:20)/10,(24:26)/10)[[x]]))))) %>% 
  arrange(dep_type,dep_val,R0,seasforce_peak) %>% rowid_to_column("par_id")
# p_table <- read_csv("simul_output/parscan/sel_parsets/sel_parsets.csv") %>% rowid_to_column("par_id")
partable <- fcn_create_partable(p_table,nstep=10, scale_age_exp=c(0.35,0.29),pop_struct=rsv_age_groups$stationary_popul,
  susc_denomin=100,susc_min=0.11,nage=11,ninf=3,rhoval=rho) %>% mutate(dep_val=ifelse(dep_type=="age",dep_val*2,dep_val)) %>%
  group_by(dep_type,R0) %>% mutate(R0_no=row_number()) %>% filter(R0_no==1)
partable_nrow <- nrow(partable)
# save the stat sol of all param sets
stat_sol_allparsets=matrix(0,nrow=(n_compartment+1)*n_age*n_inf,ncol=nrow(partable))
# NPI dates
npi_dates=as.Date(c("2020-03-26","2021-05-17"))
# width of season (from peak)
seasforc_width_wks<-8
# length of simulations
simul_length_yr=arg_vals_section_length_simuldur[3]
# agegroup indices for maternal immunity
mat_imm_flag <- TRUE
mat_imm_inds<-list(fun_sub2ind(i_inf=1,j_age=1,"R",c("S","I","R"),11,3),fun_sub2ind(i_inf=c(1,2,3),j_age=9,"R",c("S","I","R"),11,3),
                   fun_sub2ind(i_inf=c(1,2,3),j_age=9,"S",c("S","I","R"),11,3))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
serial_loop=TRUE; if (!serial_loop) {all_sum_inf_epiyear_age=list(); cl=parallel::makeCluster(6); registerDoParallel(cl)}
# as.Date <- base::as.Date
# LOOP
# tm<-proc.time()
l_sel_rows<-split(1:partable_nrow,ceiling((1:partable_nrow)/round(partable_nrow/arg_vals_section_length_simuldur[2])))
k_start_end=l_sel_rows[[arg_vals_section_length_simuldur[1]]]
print(paste0(c("this is rows",k_start_end),collapse=" "))
partable <- partable[k_start_end,]; ode_sols <- data.frame()
all_sum_inf_epiyear_age <- data.frame(); df_cases_infs_all <- data.frame()
for (k_par in 1:nrow(partable)){ # nrow(partable)
  # assign SUSCEPT params
  delta_primary<-as.numeric(partable[k_par,] %>% ungroup() %>% select(contains("delta")))
  delta_susc=sapply(1:n_age, function(x) {delta_primary/((partable$agedep_val[k_par]^(x-1)))})
  # set initial conds
  initvals_sirs_model <- fcn_set_initconds(rsv_age_groups$stationary_popul,init_set=c("previous","fromscratch")[2],
                                    init_cond_src=c("output","file")[1],NA,init_seed=10,seed_vars="all",filename="")
  # set length of simulation and seasonality
  l_seas<-fun_shutdown_seasforc(npi_dates, years_pre_post_npi=c(simul_length_yr,0),
            season_width_wks=seasforc_width_wks,init_mt_day="06-01",ifelse(grepl("exp",partable$dep_type[k_par]),45,49),
            forcing_above_baseline=partable$seasforce_peak[k_par], npireduc_strength=0.5)
  if (serial_loop){ g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% l_seas} else { 
    timesteps=l_seas[[2]]; simul_start_end=l_seas[[3]]; forcing_vector_npi=l_seas[[4]] }
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
  df_cases_infs <- fcn_process_odesol_incid(ode_solution,n_age,n_inf,n_compartment,simul_start_end) %>%
    mutate(partable$par_id[k_par])
  # print progress
  print(paste0(round(1e2*k_par/nrow(partable)),"% , dep_val=",partable$dep_val[k_par],", R0=",round(partable$R0[k_par],2),
      ", suscept=",paste0(round(delta_primary,3),collapse=","),", seas peak rel. baseline=",(partable$seasforce_peak[k_par]+1)*100,"%") )
  sel_t_point <- unique((df_cases_infs %>% filter(date == as.Date("2019-07-01")))$t)
  stat_sol_allparsets[,k_par] <- matrix(ode_solution[sel_t_point,2:ncol(ode_solution)]) # nrow(ode_solution)-1
  # final population (it's stationary, shldnt change)
  final_pop <- data.frame(agegroup=1:n_age,final=round(unlist(lapply(lapply((0:(n_age-1))*(n_inf*n_compartment), 
      function(x) (1:(n_inf*n_compartment))+x), 
      function(x_sum) sum(stat_sol_allparsets[1:(n_age*n_inf*n_compartment),k_par][x_sum])))))
  ode_sols <- bind_rows(ode_sols,data.frame(ode_solution[nrow(ode_solution)-1,2:ncol(ode_solution)]) %>% 
      mutate(var_id=row_number(),par_id=partable$par_id[k_par],dep_type=partable$dep_type[k_par],R0=partable$R0[k_par]))
  # full simul output
  # calc attack rates
  sum_inf_epiyear_age <- left_join(df_cases_infs %>% mutate(year=year(date),epi_year=ifelse(date>ymd(paste(year(date),"-07-01")),
             year(date),year(date)-1),in_out_season=ifelse(week(date)<=9 | week(date)>=41,"in","out")) %>%
  group_by(epi_year,agegroup) %>% summarise(inf_tot=round(sum(value)),inf_in_seas=round(sum(value[in_out_season=="in"])),
  max_incid_week=mean(week(date[value==max(value,na.rm=T)]))) %>% group_by(agegroup) %>% 
    filter(epi_year>min(epi_year)),final_pop,by="agegroup") %>% mutate(attack_rate_perc=100*inf_tot/final,
           seas_share=inf_in_seas/inf_tot,dep_val=partable$dep_val[k_par],par_id=partable$par_id[k_par],
           seasforce_peak=partable$seasforce_peak[k_par],dep_type=partable$dep_type[k_par],R0=partable$R0[k_par])
  # store outputs
  all_sum_inf_epiyear_age=bind_rows(all_sum_inf_epiyear_age,sum_inf_epiyear_age)
    if (k_par==nrow(partable)) { all_sum_inf_epiyear_age <- left_join(all_sum_inf_epiyear_age %>%
            mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup])),estim_attack_rates,by="agegroup_name") %>% 
        mutate(attack_rate_check=ifelse(attack_rate_perc>=min_est&attack_rate_perc<=max_est,T,F),
               seas_share_check=ifelse(seas_share>seas_conc_lim,T,F)) }
  df_cases_infs_all=bind_rows(df_cases_infs_all,df_cases_infs)
} # end loop
# rename column name
colnames(ode_sols)[1] <- "value"
# 
par_id_vals=partable$par_id

# summary simul output
write_csv(all_sum_inf_epiyear_age,paste0("simul_output/parscan/parallel/summ_parsets_start",
                                         paste0(par_id_vals,collapse="_"),".csv"))
# dynamics
write_csv(df_cases_infs_all %>% select(!name),
          paste0("simul_output/parscan/parallel/dyn_parsets_start",paste0(par_id_vals,collapse="_"),".csv") )
# init conds
write_csv(ode_sols,paste0("simul_output/parscan/parallel/init_cond_start",paste0(par_id_vals,collapse="_"),".csv"))
