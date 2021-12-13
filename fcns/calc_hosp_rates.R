# hospitalisations data
uk_rsv_hospitalisation_estimate <- read_csv("repo_data/uk_rsv_hospitalisation_estimate.csv")
# scatter
n_source<-length(unique(uk_rsv_hospitalisation_estimate$source))
hosp_coef <- data.frame(matrix(0,nrow=n_source,ncol=2)) %>% rename(a=X1,r=X2)
for (k_source in 1:n_source) {
  source_name <- unique(uk_rsv_hospitalisation_estimate$source)[k_source]
  x_midpoints <- (uk_rsv_hospitalisation_estimate %>% filter(source %in% source_name))$midpoint
  y_vals <- (uk_rsv_hospitalisation_estimate %>% filter(source %in% source_name))$rate_per_100e3_person_year
  # if (k_source!=3|k_source!=5){
  #   coeffs <- coef(nls(y_vals ~ a*exp(r*x_midpoints), start=list(a=5e3,r=-1))) } else {
  lin_coef <- coef(lm(log_y_vals~age, data = data.frame(log_y_vals=log(y_vals),age=x_midpoints)))
  coeffs <- as.numeric(c(exp(lin_coef[1]),lin_coef[2]))
  # }
  hosp_coef[k_source,] <- coeffs  }

fit_x<-list((0:71)/4,(72:360)/4)
fit_vals_under18 <- bind_rows(lapply(1:nrow(hosp_coef), 
                                     function(k) data.frame(source=unique(uk_rsv_hospitalisation_estimate$source)[k],
                    age=fit_x[[ifelse(k<4,1,2)]],fitvals=hosp_coef[k,"a"]*exp(hosp_coef[k,"r"]*fit_x[[ifelse(k<4,1,2)]]) ))) %>%
  group_by(age) %>% mutate(average_fit=mean(fitvals)) %>% ungroup() %>% mutate(upper_limit=age+1/4) 

# ggplot(uk_rsv_hospitalisation_estimate,aes(x=midpoint,y=rate_per_100e3_person_year,color=source)) + geom_point() + 
#   geom_line(data=fit_vals_under18,aes(x=age,y=fitvals,color=source)) + ylab("hospitalisations/100e3 population") +
#   scale_x_continuous(expand=expansion(0.01,0)) + theme_bw() + scale_y_log10(expand=expansion(0.01,0))
# # save
# ggsave(paste0("data//uk_rsv_hospitalisation_fit_5sources_log10.png"),width=32,height=20,units="cm")

# join popul estimates from ONS
ons_2020_midyear_estimates_uk <- read_csv("repo_data/ons_2020_midyear_estimates_uk.csv") %>% 
  mutate(age_num=as.numeric(gsub("\\+","",age)))
fit_vals_under18_aver <- fit_vals_under18 %>% pivot_wider(names_from=source,values_from=fitvals)
fit_vals_under18_aver <- fit_vals_under18_aver %>% 
  mutate(size=ons_2020_midyear_estimates_uk$value[findInterval(
    fit_vals_under18_aver$age,ons_2020_midyear_estimates_uk$age_num)]/4,size=ifelse(age==max(age),size*4,size),
    perc_pop=size/sum(size))
# plot rates
# ggplot(fit_vals_under18_aver %>% pivot_longer(!c(age,upper_limit,size,perc_pop,average_fit)),
#        aes(x=age,y=average_fit))+geom_line()+geom_point(fill=NA,shape=21)+ylab("hospitalisations/100K population/year")+theme_bw() 
# + scale_y_log10()

# calculate averages for model age brackets
hosp_rates <- data.frame()
for (k_hosp in 1:nrow(rsv_age_groups)) {
  floor_ind <- which(fit_vals_under18_aver$age %in% rsv_age_groups$age_low[k_hosp])
  ceil_ind <- which(fit_vals_under18_aver$age %in% 
                      (rsv_age_groups$age_low[k_hosp]+rsv_age_groups$duration[k_hosp]))-1
  ceil_ind <- ifelse(rsv_age_groups$age_low[k_hosp]+rsv_age_groups$duration[k_hosp]>max(fit_vals_under18_aver$age),
                     nrow(fit_vals_under18_aver),ceil_ind)
  age_bracket_inds<-floor_ind:ceil_ind
  hosp_rates[k_hosp,1] <- rsv_age_groups$agegroup_name[k_hosp]
  hosp_rates[k_hosp,2] <- round(sum(fit_vals_under18_aver$average_fit[age_bracket_inds]*(
    fit_vals_under18_aver$perc_pop[age_bracket_inds]/sum(fit_vals_under18_aver$perc_pop[age_bracket_inds]))),1) }
# shape dataframe
hosp_rates <- hosp_rates %>% rename(agegroup_name=V1,hosp_rate_per_100K=V2) %>%
  mutate(hosp_number_from_pop_estim=round(hosp_rate_per_100K*rsv_age_groups$value/1e5))
# if (sum(hosp_rates$hosp_number)>50e3){ corr_r <- 50e3/sum(hosp_rates$hosp_number)
#   hosp_rates[,2:3] <- hosp_rates[,2:3]*corr_r }

hosp_probabilities <- read_csv("repo_data/hosp_probabilities_hodgson_adjusted.csv") 
# %>% mutate(size=rsv_age_groups$value,
#   hosp_num=size*(estim_attack_rates$median_all_inf/100)*prob_hosp_per_infection,agegroup_name=gsub("_","-",agegroup_name)) %>%
#  rename(hosp_num_from_per_inf_prob=hosp_num)

# left_join(hosp_probabilities,hosp_rates,by="agegroup_name") %>% 
#   select(c(agegroup_name,hosp_num_from_per_inf_prob,hosp_number_from_pop_estim)) %>% pivot_longer(!agegroup_name)