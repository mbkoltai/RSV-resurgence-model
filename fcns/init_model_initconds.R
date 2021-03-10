if (exists("df_ode_solution")){initsrc="output"} else {initsrc="file"; df_ode_solution=c()}
initvals_sirs_model=fcn_init_susc_vals(stationary_init=TRUE,from_file_or_output=initsrc,simul_output=df_ode_solution,
                                       susc_vars_inds,agegr_sizes=rsv_age_groups$value,sim_filepath="simul_output/statsol_10years.RDS")
# INITIAL INFECTION (taking stationary sol normally contains some Inf.s so no need to re-seed it)
initvals_sirs_model[inf_vars_inds[[1]][1]]=10 # all first infection groups: sapply(inf_vars_inds, '[[',1)
### integrate ODE --------------------------------------------------------
# deSolve input

# solve with interpolation
npi_year=4; preseas_npi_on=2; postseas_npi_off=2 # shutdown_scale=0.2
basal_rate=0.15; peak_week=49; season_width=3.8 # (in weeks)
shutdown_list=fun_shutdown_seasforc(timesteps,elem_time_step,basal_rate,npi_year,peak_week,season_width,
                                    preseas_npi_on,postseas_npi_off,n_prec=0.01,st_devs=2:3,n_sd=2)
forcing_vector_npi=shutdown_list[[1]]
# RUN
params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,forcing_vector_npi,elem_time_step,delta_susc)
approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
ptm<-proc.time(); ode_solution<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc_interpol,parms=params); proc.time()-ptm
df_ode_solution=fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups,neg_thresh=-0.001)[[1]]
# 
if (exists("df_ode_solution") & ncol(df_ode_solution)>0){initsrc="output"} else {initsrc="file"; df_ode_solution=c()}
initvals_sirs_model=fcn_init_susc_vals(stationary_init=TRUE,from_file_or_output=initsrc,simul_output=df_ode_solution,
                                       susc_vars_inds,agegr_sizes=rsv_age_groups$value,sim_filepath="simul_output/statsol_10years.RDS")
# INITIAL INFECTION (taking stationary sol normally contains some Inf.s so no need to re-seed it)
# initvals_sirs_model[inf_vars_inds[[1]][1]]=10 # all first infection groups: sapply(inf_vars_inds, '[[',1)