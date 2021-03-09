source('RSV_model_functions.R')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SET PARAMETERS --------------------------------------------------------
# selected country
country_sel="Germany"
# time resolution (in days)
elem_time_step=0.5
# population data
standard_age_groups <- fun_cntr_agestr("Germany","2015",seq(0,75,5),c(seq(4,74,5),99))
# RSV age groups (population data from wpp2019)
rsv_age_groups_low=c(0,0.5,1,1.5, 2,3,4, 5,10,15, 20); rsv_age_group_sizes=c(rep(0.4,4),rep(0.9,3),rep(4,3),79)
rsv_age_groups=fun_rsv_agegroups(standard_age_groups,rsv_age_groups_low,rsv_age_group_sizes)
# population by age group
N_tot=sum(rsv_age_groups$value) # S_by_age=rsv_age_groups$value
# number of age groups, reinfections and variables (S,I,R)
n_age=nrow(rsv_age_groups); varname_list=c('S','I','R'); n_compartment=length(varname_list); n_inf=3
dim_sys=n_age*n_compartment*n_inf; n_days_year=365
# query variables: fun_sub2ind(1:3,11,"R",varname_list,n_age,n_inf)
# force of infection terms
# linear indices of the I & S variables
l_inf_susc=fun_inf_susc_index_lists(n_age,n_inf,varname_list); inf_vars_inds=l_inf_susc[[1]]; susc_vars_inds=l_inf_susc[[2]]
# CONTACT MATRIX
# contact matrix from covidm ("home","work","school","other")
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
list_contmatrs=fun_covidm_contactmatrix(country_sel="Germany",currentdir_path=currentdir_path,cm_path=cm_path)
C_m_full=Reduce('+',list_contmatrs)
# make matrix symmetric
C_m_fullrecipr=fun_recipr_contmatr(C_m_full,age_group_sizes=standard_age_groups$values)
# create for our age groups
C_m=fun_create_red_C_m(C_m_fullrecipr,rsv_age_groups,
                       orig_age_groups_duration=standard_age_groups$duration,orig_age_groups_sizes=standard_age_groups$values)
# make it reciprocal for the larger group
C_m=fun_recipr_contmatr(C_m,age_group_sizes=rsv_age_groups$value)
# bc of reinfections we need to input contact matrix repeatedly
contmatr_rowvector=t(do.call(cbind, lapply(1:nrow(C_m), function(x){diag(C_m[x,]) %*% matrix(1,n_age,n_inf)})))

# build kinetic matrix --------------------------------------------------------
# WANING (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j
omega=1/1.25e2 # 1/runif(1,60,200)
# RECOVERY
rho=1/7; # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,rsv_age_groups)