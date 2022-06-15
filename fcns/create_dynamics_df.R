# create dataframes for dynamics 

all_dynamics_accepted <- bind_rows(lapply(
  list.files("simul_output/2e4_parsets/all_dyn",pattern="dyn_parsets.*csv",full.names=T), 
  function(x) read_csv(x) %>% 
    filter(par_id %in% partable_regular_dyn$par_id) # %>%
  # filter(par_id %in% sample(unique(par_id),size=sample(max(length(unique(par_id)),94:96),size=1))) 
) ) 
all_dynamics_accepted <- all_dynamics_accepted %>% 
  filter(par_id %in% sample(x=unique(all_dynamics_accepted$par_id),size=2e3,replace=F) )
write_csv(all_dynamics_accepted,"simul_output/2e4_parsets/all_dynamics_accepted_2e3.csv")
# we need to also take rejected parsets to make comparison (likelihood vs filtering)
all_dynamics_rejected <- bind_rows(lapply(
  list.files("simul_output/2e4_parsets/all_dyn",pattern="dyn_parsets",full.names=T),
  function(x) read_csv(x) %>% filter(!par_id %in% unique(results_summ_all_reg_dyn$par_id)) %>%
    filter(par_id %in% sample(unique(par_id),sample(c(83,84),1),replace=F) ) ))
# save
write_csv(all_dynamics_rejected,"simul_output/2e4_parsets/all_dynamics_rejected_2e3.csv") 