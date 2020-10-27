# RSV_model_functions
# generate linear index from 2-dim index (infection-age) ----------------------------------------------------------
fun_sub2ind=function(i_inf,j_age,varname,varname_list,n_var,n_age,n_inf){
  varnum=which(varname_list %in% varname); k=(j_age-1)*n_var*n_inf + (varnum-1)*n_inf + i_inf; k }

# generate all model names ----------------------------------------------------------
fun_sirs_varnames=function(varname_list,n_age,n_inf){
  array(sapply(1:n_age, function(x_age) {sapply(varname_list, function(x) {paste0(paste0(x,'_',1:n_inf),'_',x_age)})} ))
}

# set up kinetic matrix ----------------------------------------------------------
fun_K_m_sirs_multiage=function(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list){
K_m=matrix(0,nrow=dim_sys,ncol=dim_sys)
# S_i_j -> S_1_1 is S, subscript=1, superscript=1. subscript: # infection, superscript= # age group
# conversion between i,j and X_k, when variables are stacked as S_i_1,I_i_1,R_i_1, S_i_2,I_i_2,R_i_2 ...
# varname_list=c('S','I','R')
# waning terms: R_i_j -> S_min(i+1,n_inf)_j
# omega=1/1e2; # 1/runif(1,60,200)
for (j_age in 1:n_age) {
  for (i_inf in 1:n_inf) { if (j_age==1 & i_inf==1) {waning_terms_source_target=data.frame()}
    wanevals=c(fun_sub2ind(i_inf,j_age,'R',varname_list,n_compartment,n_age,n_inf),
               fun_sub2ind(min(i_inf+1,n_inf),j_age,'S',varname_list,n_compartment,n_age,n_inf))
    # waning_terms_source_target=rbind(waning_terms_source_target,wanevals)
    K_m[wanevals[2],wanevals[1]]=omega } }
# aging terms between compartments: S_i_j -> S_i_(j+1), R_i_j -> S_(i+1)_(j+1)
duration_age_groups=rep(1,n_age); # eta_a=1/(365*d_a); 
for (j_age in 1:(n_age-1)) {
  for (i_inf in 1:n_inf) { if (j_age==1 & i_inf==1) {aging_terms_source_target=data.frame()}
    agevals=rbind(c(fun_sub2ind(i_inf,j_age,'S',varname_list,n_compartment,n_age,n_inf),
                    fun_sub2ind(i_inf,j_age+1,'S',varname_list,n_compartment,n_age,n_inf)),
                  c(fun_sub2ind(i_inf,j_age,'R',varname_list,n_compartment,n_age,n_inf),
                    fun_sub2ind(min(i_inf+1,n_inf),j_age+1,'S',varname_list,n_compartment,n_age,n_inf))) ### end of rbind
    aging_terms_source_target=rbind(aging_terms_source_target,agevals); d_a=duration_age_groups[j_age] } }
for (k in 1:nrow(aging_terms_source_target)) {K_m[aging_terms_source_target[k,2],aging_terms_source_target[k,1]]=1/(365*d_a)}

# recovery terms
# rho=1/6; # 1/rho=rweibull(1, shape=4.1,scale=8.3)
for (j_age in 1:n_age) {
  for (i_inf in 1:n_inf) { if (j_age==1 & i_inf==1) {recov_terms_source_target=data.frame()}
    recov_vals=c(fun_sub2ind(i_inf,j_age,'I',varname_list,n_compartment,n_age,n_inf),
                 fun_sub2ind(i_inf,j_age,'R',varname_list,n_compartment,n_age,n_inf))
    recov_terms_source_target=rbind(recov_terms_source_target,recov_vals)
    K_m[recov_vals[2],recov_vals[1]]=rho } }

# diagonal terms
# outflow terms that represent aging 'out of the model' from the highest age groups
n_days_year=365
for (j_age in n_age) {
  for (i_inf in 1:n_inf) { if (i_inf==1) {ageout_terms=data.frame()}
    ageout_terms=rbind(ageout_terms, rbind(fun_sub2ind(i_inf,j_age,'S',varname_list,n_compartment,n_age,n_inf),
                                           fun_sub2ind(i_inf,j_age,'R',varname_list,n_compartment,n_age,n_inf))) } }
for (k in 1:nrow(ageout_terms)) {K_m[ageout_terms[k,1],ageout_terms[k,1]]=-1/(n_days_year*d_a) }
# diagonal terms balancing the outgoing terms, these are the (sums of the off diagonal terms)x(-1) 
if (any(diag(K_m)==0)){ diag(K_m)=diag(K_m)-colSums(K_m-diag(diag(K_m))) }
#### return output
K_m
}