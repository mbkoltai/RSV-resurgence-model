C_m_HU=Reduce('+',cm_parameters_SEI3R("Hungary")$pop[[1]]$matrices); 
C_m_HU_symm=C_m_HU
HU_pop_agegroups=cm_parameters_SEI3R("Hungary")$pop[[1]]$size
# (popF[popF$name %in% "Hungary","2020"] + popM[popM$name %in% "Hungary","2020"])*1e3
# HU_pop_agegroups=c(HU_pop_agegroups[1:15],sum(HU_pop_agegroups[16:length(HU_pop_agegroups)]))
# make matrix symmetric
for (k in 1:nrow(all_perms)) { i=all_perms[k,1]; j=all_perms[k,2]; pair_popul=HU_pop_agegroups[j]+HU_pop_agegroups[i]
  # sym[i,j] <- ( m[i,j] + m[j,i]*( N[j]/N[i] ) )/2
  C_m_HU_symm[i,j]=(C_m_HU[i,j]+C_m_HU[j,i]*HU_pop_agegroups[j]/HU_pop_agegroups[i])/2
}
# 
HU_age_groups=list(1,2:3,4:6,7:12,13:14,15,16); HU_pop_coarse_agegroups=sapply(HU_age_groups, function(x){sum(HU_pop_agegroups[x])})
C_M_summed=matrix(0,length(HU_age_groups),length(HU_age_groups))
for (k_row in 1:length(HU_age_groups)) {
  for (k_col in 1:length(HU_age_groups)) {
    # popul_factor=(HU_pop_coarse_agegroups[k_row]/sum(HU_pop_coarse_agegroups))
  C_M_summed[k_row,k_col]=sum(C_m_HU_symm[HU_age_groups[[k_row]],HU_age_groups[[k_col]]]) # /popul_factor
  }
}