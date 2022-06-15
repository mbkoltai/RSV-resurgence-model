partable <- data.frame(randomLHS(n=32e3,k=7))
colnames(partable) <- c("exp_dep","age_dep","seasforc_width_wks","R0","seasforce_peak","omega","peak_week")
# convert to relevant ranges, columns: 
# 1) expdep 2) agedep 3) peak_width 4) R0 5) peak_height 6) waning 7) peak week
partable[,1] <- qunif(partable[,1],min=3/10,max=1.25) # 1) expdep
partable[,2] <- qunif(partable[,2],min=1/15,max=1/3) # 2) agedep
partable[,3] <- qgamma(partable[,3],shape=10,rate=2) # 3) peak_width
partable[,4] <- qgamma(partable[,4],shape=14,rate=8) # 4) R0
partable[,5] <- qunif(partable[,5],min=0.2,max=1.5) # 5) peak height
partable[,6] <- qnorm(partable[,6],mean=350,sd=50); partable[,6] <- 1/partable[,6] # 6) waning
partable[,7] <- qunif(partable[,7],min=43,max=48) # 7) peak week
if (any(partable<0)) { partable=abs(partable); message("negative values") }

# initial sampling showed that accepted parsets satisfy the condition: exp_dep < 1.65 - 4.5*age_dep
partable <- partable %>% filter(exp_dep+4.5*age_dep<1.85)
# calculating the susceptibility parameter (delta_i_j)
l_delta_susc <- 1:nrow(partable) %>%
  lapply( function(n_p) {
    sapply(1:n_age, function(x) { (1 * exp(-partable$exp_dep[n_p] * (1:3))) / (exp(partable$age_dep[n_p]*x)) }) } )
# create partable with scaling parameters (this parameter scales 'expdep' and 'agedep' to have the desired R0 value)
partable <- partable %>% 
  mutate(par_id=row_number(),
         const_delta=R0/unlist(lapply(l_delta_susc, function(x) R0_calc_SIRS(C_m,x,rho,n_inf)))) %>%
  relocate(par_id,.before=exp_dep)
# clear list
rm(l_delta_susc)
# SAVE filtered parameter table
# write_csv(partable,"repo_data/partable_full_lhs.csv")
