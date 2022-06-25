
age_exp_dep_uniqvals <- list(exp_dep=seq(3/10,1.25,0.19),
                             age_dep=seq(1/15,1/3,1/15),age=1:11,exp=1:3) %>% 
  expand.grid %>% bind_rows %>% 
  mutate(suscept_unscaled=exp(-(exp_dep*exp+age_dep*age)))

age_exp_dep_uniqvals <- age_exp_dep_uniqvals %>% 
  mutate(const_delta=1/unlist(lapply(lapply(1:nrow(age_exp_dep_uniqvals), 
                                            function(n_p) { sapply(1:n_age,function(x) {
                                              (1*exp(-age_exp_dep_uniqvals$exp_dep[n_p]*(1:3)))/
                                                (exp(age_exp_dep_uniqvals$age_dep[n_p]*x))})}),
                                     function(x) R0_calc_SIRS(C_m,x,rho,n_inf))),
         susc_scaled=suscept_unscaled*const_delta)

# PLOT
ggplot(age_exp_dep_uniqvals %>% 
         filter(exp_dep %in% seq(3/10,1.25,0.19)[c(1,3,5)] 
                & age_dep %in% seq(1/15,1/3,1/15)[c(1,3,5)]) %>%
         mutate(age_dep=round(age_dep,2)) %>%
         rename(`exposure-dependence`=exp_dep,`age-dependence`=age_dep) %>% 
         mutate(age=factor(rsv_age_groups$agegroup_name[age],
                           levels=unique(rsv_age_groups$agegroup_name))) )  + 
  geom_line(aes(x=age,color=factor(exp),group=exp,y=susc_scaled),size=1.06) + 
  facet_grid(`exposure-dependence`~`age-dependence`,
             labeller=labeller(`exposure-dependence`=label_both,
                               `age-dependence`=label_both)) + 
  scale_y_log10() + labs(color="exposure") + 
  xlab("age group") + ylab(expression(delta[exp]^(age))) +
  theme_bw() + standard_theme
# ggsave
ggsave(here(figs_folder,"SI_FIG2.png"),width=22,height=18,units="cm")