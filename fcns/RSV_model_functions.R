# RSV_model_functions
# Model description ----------------------
# without reinfections the model would be: dx/dt = [-id_matr;id_matr;0]*diag(I_vect)*C_m*S_vect + K_m*x(t)
# with reinfections more complicated, abstractly denoted: 
# dx/dt = F(I_vect,C_m,S_vect) + K_m*x(t)
# I_vect: vector of infectious compartments, S_vect: vector of susceptible compartments, C_m: contact matrix
# K_m: kinetic matrix for linear terms (aging, recovery, waning of immunity)
#
# steps of constructing infection terms (lambda):
# linear list of state variables
# X_vars=matrix(abs(rnorm(n_inf*n_age*n_compartment)),n_inf*n_age*n_compartment,1) # matrix(0,n_inf*n_age*n_compartment,1)
# inf_vars_stacked=do.call(cbind,lapply(inf_vars_inds, function(x){X_vars[x]})) # do.call(cbind,inf_vars_inds)
# inf_vars_stacked_fullsize=t(matrix(1,1,n_inf)%*%inf_vars_stacked) 
# # full lambda column vector
# lambda_vect=diag(array(delta_susc))%*%contmatr_rowvector%*%inf_vars_stacked_fullsize
# # infection vector
# infection_vect=diag(X_vars[unlist(susc_vars_inds)])%*%lambda_vect
# 
# # put together RHS of ODEs (this is to test, we need to do it within ODE function)
# # dX/dt = F(delta,S,lambda) + K_m*X
# F_vect=matrix(0,dim_sys,1)
# F_vect[c(unlist(susc_vars_inds),unlist(inf_vars_inds))]=rbind(-infection_vect,infection_vect)
# rhs_odes=birth_term + F_vect + K_m%*%X_vars
# sirs_seasonal_forc <- function(t,X,parms){ 
#   birth_term=parms[[1]];K_m=parms[[2]];contmatr_rowvector=parms[[3]];inf_vars_inds=parms[[4]];susc_vars_inds=parms[[5]]
#   forcing_vector_npi=parms[[6]]; elem_time_step=parms[[7]]; # event_vector=parms[[8]]
#   # stack I vars
#   inf_vars_stacked=do.call(cbind,lapply(inf_vars_inds, function(x){X[x]}))
#   inf_vars_stacked_fullsize=t(matrix(1,1,n_inf)%*%inf_vars_stacked)
#   lambda_vect=diag(forcing_vector_npi[(t/elem_time_step)+1]*array(delta_susc))%*%contmatr_rowvector%*%inf_vars_stacked_fullsize
#   infection_vect=diag(X[unlist(susc_vars_inds)])%*%lambda_vect
#   F_vect=matrix(0,dim_sys,1); F_vect[c(unlist(susc_vars_inds),unlist(inf_vars_inds))]=rbind(-infection_vect,infection_vect)
#   dXdt=birth_term + F_vect + K_m%*%X; list(dXdt) }
#
### clinical fraction as a fcn of age and exposure -----------------------------------
# list_symptom_agegroups=list(1:2,3:7,8:9,10:11); prop_symptom=1-c(mean(rbeta(1e3,shape1=3,shape2=30)),
#   mean(rbeta(1e3,shape1=9,shape2=43)),mean(rbeta(1e3,shape1=38,shape2=35)),mean(rbeta(1e3,shape1=36,shape2=11))) 
#
# SUSCEPT
# rbeta(35.583,11.417)~0.75; B(22.829,3.171)~0.9; B(6.117,12.882)~0.32

# set plotting theme
standard_theme=theme(# panel.grid=element_line(linetype="solid",colour="black",size=0.1),
                     plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),
                     axis.text.y=element_text(size=9),
                     axis.title=element_text(size=14), text=element_text(family="Calibri"))

# standard cols
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Set up SIRS ODE model --------------------------------------------------------
# see description in "RSV_model_functions.R"
## with seasonal forcing by subsetting t (SLOW!!!)
# sirs_seasonal_forc <- function(t,X,parms){
#   birth_term=parms[[1]];K_m=parms[[2]];contmatr_rowvector=parms[[3]];inf_vars_inds=parms[[4]];susc_vars_inds=parms[[5]]
#   forcing_vector_npi=parms[[6]]; elem_time_step=parms[[7]]; delta_susc=parms[[8]]
#   # stack I vars
#   inf_vars_stacked=do.call(cbind,lapply(inf_vars_inds, function(x){X[x]}))
#   inf_vars_stacked_fullsize=t(matrix(1,1,n_inf)%*%inf_vars_stacked)
#   lambda_vect=diag(forcing_vector_npi[(t/elem_time_step)+1]*array(delta_susc))%*%contmatr_rowvector%*%inf_vars_stacked_fullsize
#   infection_vect=diag(X[unlist(susc_vars_inds)])%*%lambda_vect
#   F_vect=matrix(0,dim_sys,1); F_vect[c(unlist(susc_vars_inds),unlist(inf_vars_inds))]=rbind(-infection_vect,infection_vect)
#   dXdt=birth_term + F_vect + K_m %*% X; list(dXdt) }

## with seasonal forcing using interpolation
sirs_seasonal_forc <- function(t,X,parms){
  birth_term=parms[[1]];K_m=parms[[2]];contmatr_rowvector=parms[[3]];inf_vars_inds=parms[[4]];susc_vars_inds=parms[[5]]
  elem_time_step=parms[[6]]; delta_susc=parms[[7]] # forcing_vector_npi=parms[[6]]; 
  # stack I vars
  inf_vars_stacked=do.call(cbind,lapply(inf_vars_inds, function(x){X[x]}))
  inf_vars_stacked_fullsize=t(matrix(1,1,n_inf)%*%inf_vars_stacked)
  lambda_vect=diag(approx_seas_forc(t)*array(delta_susc))%*%contmatr_rowvector%*%inf_vars_stacked_fullsize 
  infection_vect=diag(X[unlist(susc_vars_inds)])%*%lambda_vect
  F_vect=matrix(0,dim_sys,1); F_vect[c(unlist(susc_vars_inds),unlist(inf_vars_inds))]=rbind(-infection_vect,infection_vect+approx_introd(t))
  dXdt=birth_term + F_vect + K_m %*% X; list(dXdt) }

# linear index from 2-dim index (infection-age) ----------------------------------------------------------
fun_sub2ind=function(i_inf,j_age,varname,varname_list,n_age,n_inf){
  varnum=which(varname_list %in% varname); k=(j_age-1)*length(varname_list)*n_inf + (varnum-1)*n_inf + i_inf; k }

# generate all model names ----------------------------------------------------------
fun_sirs_varnames=function(varname_list,n_age,n_inf){
  array(sapply(1:n_age, function(x_age) {sapply(varname_list, function(x) {paste0(paste0(x,'_',1:n_inf),'_',x_age)})} ))
}

# set up kinetic matrix ----------------------------------------------------------
fun_K_m_sirs_multiage=function(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,rsv_age_groups){
K_m=matrix(0,nrow=dim_sys,ncol=dim_sys)
# S_i_j -> S_1_1 is S, subscript=1, superscript=1. subscript: # infection, superscript= # age group
# varname_list=c('S','I','R') # conversion between i,j and X_k, when variables are stacked as S_i_1,I_i_1,R_i_1, S_i_2,I_i_2,R_i_2 ...
# waning terms: omega=1/1e2. flow: R_i_j -> S_min(i+1,n_inf)_j;  eg. R_1_1 -> S_2_1
for (j_age in 1:n_age) {
  for (i_inf in 1:n_inf) { if (j_age==1 & i_inf==1) {waning_terms_source_target=data.frame()}
    wanevals=c(fun_sub2ind(i_inf,j_age,'R',varname_list,n_age,n_inf),
               fun_sub2ind(min(i_inf+1,n_inf),j_age,'S',varname_list,n_age,n_inf))
    # waning_terms_source_target=rbind(waning_terms_source_target,wanevals)
    K_m[wanevals[2],wanevals[1]]=omega } }
# aging terms between AGE GROUPS: S_i_j -> S_i_(j+1), R_i_j -> S_(i+1)_(j+1), I_i_j->I_(i+1)_j
for (j_age in 1:(n_age-1)) {
  for (i_inf in 1:n_inf) { if (j_age==1 & i_inf==1) {aging_terms_source_target=data.frame()}
  agevals=rbind( 
  # S_i_j -> S_i_(j+1)
  c(fun_sub2ind(i_inf,j_age,'S',varname_list,n_age,n_inf),fun_sub2ind(i_inf,j_age+1,'S',varname_list,n_age,n_inf),j_age ),
  # R_i_j -> S_(i+1)_(j+1) OR R_i_j -> R_i_(j+1)
  # c(fun_sub2ind(i_inf,j_age,'R',varname_list,n_age,n_inf),fun_sub2ind(min(i_inf+1,n_inf),j_age+1,'S',varname_list,n_age,n_inf),j_age)
  c(fun_sub2ind(i_inf,j_age,'R',varname_list,n_age,n_inf),fun_sub2ind(i_inf,j_age+1,'R',varname_list,n_age,n_inf),j_age),
  # I_i_j -> I_i_(j+1)
  c(fun_sub2ind(i_inf,j_age,'I',varname_list,n_age,n_inf),fun_sub2ind(i_inf,j_age+1,'I',varname_list,n_age,n_inf),j_age)
  )
  # bundle
  aging_terms_source_target=rbind(aging_terms_source_target,agevals) } }
for (k in 1:nrow(aging_terms_source_target)) {
  duration_scale=rsv_age_groups$duration[aging_terms_source_target[k,3]]
  # print(c(aging_terms_source_target[k,2],aging_terms_source_target[k,1],duration_scale))
  K_m[aging_terms_source_target[k,2],aging_terms_source_target[k,1]]=1/(365*duration_scale)}
# recovery terms
# rho=1/6; # 1/rho=rweibull(1, shape=4.1,scale=8.3)
for (j_age in 1:n_age) {
  for (i_inf in 1:n_inf) { if (j_age==1 & i_inf==1) {recov_terms_source_target=data.frame()}
    recov_vals=c(fun_sub2ind(i_inf,j_age,'I',varname_list,n_age,n_inf),
                 fun_sub2ind(i_inf,j_age,'R',varname_list,n_age,n_inf))
    recov_terms_source_target=rbind(recov_terms_source_target,recov_vals)
    K_m[recov_vals[2],recov_vals[1]]=rho } }

# diagonal terms
# outflow terms that represent 'aging out' of the model from the highest age groups
for (j_age in n_age) {
  for (i_inf in 1:n_inf) { if (i_inf==1) {ageout_terms=data.frame()}
    ageout_terms=rbind(ageout_terms, rbind(fun_sub2ind(i_inf,j_age,'S',varname_list,n_age,n_inf),
                                           # if infections have the aging term on them too
                                           fun_sub2ind(i_inf,j_age,'I',varname_list,n_age,n_inf),
                                           fun_sub2ind(i_inf,j_age,'R',varname_list,n_age,n_inf))) } }
for (k in 1:nrow(ageout_terms)) {K_m[ageout_terms[k,1],ageout_terms[k,1]]=-1/(365*rsv_age_groups$duration[nrow(rsv_age_groups)]) }
# diagonal terms balancing the outgoing terms, these are the (sums of the off diagonal terms)x(-1) 
diag(K_m)=diag(K_m)-colSums(K_m-diag(diag(K_m)))

#### return output
# print(diag(K_m))
K_m
}

### R0 calculation --------------
R0_calc_SIRS <- function(C_m,delta_susc_prop,rho,n_inf) {
  # po = parameters$pop[[population]];
  # dIp = sum(po$dIp * seq(0, by = parameters$time_step, length.out = length(po$dIp)));
  # dIs = sum(po$dIs * seq(0, by = parameters$time_step, length.out = length(po$dIs)));
  # dIa = sum(po$dIa * seq(0, by = parameters$time_step, length.out = length(po$dIa)));
  # cm = Reduce('+', mapply(function(c, m) c * m, po$contact, po$matrices, SIMPLIFY = F));
  susc_matrs=lapply(1:n_inf, function(x) {matrix(rep(delta_susc_prop[x,],ncol(delta_susc_prop)),ncol = ncol(delta_susc_prop))})
  # ngm = po$u * t(t(cm) * (    po$y * (po$fIp * dIp + po$fIs * dIs) + (1 - po$y) * po$fIa * dIa) )
  ngm=Reduce('+',lapply(susc_matrs, function(x) {x* C_m}))*(1/rho)*(1/n_inf)
  abs(eigen(ngm)$values[1])
}

### ode w/o seasonal forcing --------------
sirs_template <- function(t,X,parms){
  birth_term=parms[[1]]; K_m=parms[[2]]; contmatr_rowvector=parms[[3]]; inf_vars_inds=parms[[4]]; susc_vars_inds=parms[[5]]
  # stack I vars
  inf_vars_stacked=do.call(cbind,lapply(inf_vars_inds, function(x){X[x]}))
  inf_vars_stacked_fullsize=t(matrix(1,1,n_inf)%*%inf_vars_stacked)
  lambda_vect=diag(array(delta_susc))%*%contmatr_rowvector%*%inf_vars_stacked_fullsize
  infection_vect=diag(X[unlist(susc_vars_inds)])%*%lambda_vect
  F_vect=matrix(0,dim_sys,1); F_vect[c(unlist(susc_vars_inds),unlist(inf_vars_inds))]=rbind(-infection_vect,infection_vect)
  dXdt=birth_term + F_vect + K_m%*%X; list(dXdt) }

# age structure of country --------------
fun_cntr_agestr=function(i_cntr,i_year,age_low_vals,age_high_vals){
  age_groups=data.frame(age_low=seq(0,75,5), age_high=c(seq(4,74,5),100))
  if (!any((.packages()) %in% "wpp2019")) {library(wpp2019)}; if (!exists("popF")) {data("pop")}
  cntr_agestr=data.frame(agegroups=popF[popF$name %in% i_cntr,"age"],values=popF[popF$name %in% i_cntr,i_year] +
                           popM[popM$name %in% i_cntr,i_year])
  agegr_truthvals=sapply(strsplit(as.character(cntr_agestr$agegroups),"-"),"[[",1) %in% age_groups$age_low
  N_tot=cntr_agestr$values[agegr_truthvals]
  N_tot[length(N_tot)]=N_tot[length(N_tot)]+sum(cntr_agestr$values[!agegr_truthvals]); N_tot=N_tot*1e3; # N_tot
  data.frame(age_low=age_low_vals, age_high=age_high_vals,values=N_tot, duration=(age_high_vals-age_low_vals)+1) %>%
    mutate(proportion=values/sum(values))
}

# create inf and susc vars indices --------------
fun_inf_susc_index_lists<-function(n_age,n_inf,varname_list){   n_compartment=length(varname_list)
  inf_vars_inds=lapply(1:n_age, function(x_age){ sapply(1:n_inf, function(x_inf){
  fun_sub2ind(x_inf,x_age,varname='I',varname_list,n_age,n_inf) }) })
# linear indices of S variables
susc_vars_inds=lapply(1:n_age, function(x_age){ sapply(1:n_inf, function(x_inf){
  fun_sub2ind(x_inf,x_age,varname='S',varname_list,n_age,n_inf) }) })
list(inf_vars_inds,susc_vars_inds)
}

### get objects larger than 1M (memory) --------------
fcn_objs_mem_use <- function(min_size){
mem_use_df=round(data.frame(unlist(sapply(ls(envir=.GlobalEnv), function(n) object.size(get(n)), simplify = FALSE)))/1e6,1)
colnames(mem_use_df)[1]<-"size (Mb)"; mem_use_df[,"objs"]=rownames(mem_use_df)
mem_use_df<-mem_use_df[order(mem_use_df$size,decreasing=T),]; rownames(mem_use_df)<-c()
mem_use_df[mem_use_df$size>min_size,c(2,1)]
}

### functn call covidm contact matrix --------------
fun_covidm_contactmatrix <- function(country_sel,currentdir_path,cm_path){if (!exists("covid_params")){
  setwd(cm_path); cm_force_rebuild<<-F; cm_build_verbose<<-T; cm_version<<-2; source(file.path(cm_path,"R","covidm.R"))
  covid_params=cm_parameters_SEI3R(country_sel); setwd(currentdir_path)} 
  # covidm_contactm=Reduce('+',covid_params$pop[[1]]$matrices)
#  "home"   "work"   "school" "other" 
covid_params$pop[[1]]$matrices}

# funcn create reciprocal matrix  --------------
fun_recipr_contmatr<-function(C_m_full,age_group_sizes){
  all_perms=permutations(n=nrow(C_m_full),r=2,repeats.allowed=T); N_tot=sum(age_group_sizes)
  C_m_full_symm=matrix(0,nrow=nrow(C_m_full),ncol=nrow(C_m_full))
  for (k in 1:nrow(all_perms)) { 
    i=all_perms[k,1]; j=all_perms[k,2]
    C_m_full_symm[i,j]=(C_m_full[i,j] + C_m_full[j,i]*(age_group_sizes[j]/age_group_sizes[i]))/2
  }
  colnames(C_m_full_symm)=colnames(C_m_full); rownames(C_m_full_symm)=rownames(C_m_full) 
  C_m_full_symm
}

# create reduced contact matrix --------------
fun_create_red_C_m=function(C_m_full,rsv_age_groups,orig_age_groups_duration,orig_age_groups_sizes){
  C_m=matrix(0,nrow=nrow(rsv_age_groups),ncol=nrow(rsv_age_groups)); rownames(C_m)=rsv_age_groups$agegroup_name
  colnames(C_m)=rsv_age_groups$agegroup_name
  for (i_row in 1:n_age){
    for (j_col in 1:n_age){
      # we are merging or splitting age groups, there are 3 possibilities for a new age group:
      # same OR smaller than (ST) OR larger than (LT) the original
      # (it is an *average* of contacts per person, and we have no resolution within age bands)
      #
      # if the 'i' group (C[i,j]) is the *same* as original or *smaller*, this (in itself) does not change the contact rate
      if (rsv_age_groups$wpp_agegroup_low[i_row]==rsv_age_groups$wpp_agegroup_high[i_row]) {
        # 'j' group same or smaller as original
        if (rsv_age_groups$wpp_agegroup_low[j_col]==rsv_age_groups$wpp_agegroup_high[j_col]) {
            f_dur=rsv_age_groups$duration[j_col]/orig_age_groups_duration[rsv_age_groups$wpp_agegroup_high[j_col]]
            C_m[i_row,j_col]=(C_m_full[rsv_age_groups$wpp_agegroup_low[i_row],rsv_age_groups$wpp_agegroup_low[j_col]])*f_dur
        } else { # if 'j' is larger than original group
          group_span=rsv_age_groups$wpp_agegroup_low[j_col]:rsv_age_groups$wpp_agegroup_high[j_col]
          agegroup_weights=orig_age_groups_sizes[group_span]/sum(orig_age_groups_sizes[group_span])
          C_m[i_row,j_col]=sum(agegroup_weights*C_m_full[i_row,group_span])
        } # end of 'i' smaller or same as original
      } else { # if 'i' in C[i,j] is a bigger age band -> weighted average of the contact rates of constituent groups
        group_span=rsv_age_groups$wpp_agegroup_low[i_row]:rsv_age_groups$wpp_agegroup_high[i_row]
        agegroup_weights=orig_age_groups_sizes[group_span]/sum(orig_age_groups_sizes[group_span])
        # if 'j' is same/smaller -> contact rate with original group proportionally divided
        if (rsv_age_groups$wpp_agegroup_low[j_col]==rsv_age_groups$wpp_agegroup_high[j_col]) {
          f_dur=rsv_age_groups$duration[j_col]/orig_age_groups_duration[rsv_age_groups$wpp_agegroup_high[j_col]]
  C_m[i_row,j_col]=sum((orig_age_groups_sizes[group_span]/sum(orig_age_groups_sizes[group_span]))*C_m_full[group_span,j_col])*f_dur
        } else {# if 'j' larger -> weighted average of the contact rates of the constituent groups
          C_m[i_row,j_col]=sum(rep(agegroup_weights,length(agegroup_weights))*unlist(
            lapply(group_span,function(x) {agegroup_weights*C_m_full[x,group_span]})))
        }
      }
      # C_m[i_row,j_col]=mean(C_m_full[rsv_age_groups$wpp_agegroup_low[i_row]:rsv_age_groups$wpp_agegroup_high[i_row],
      #                                rsv_age_groups$wpp_agegroup_low[j_col]:rsv_age_groups$wpp_agegroup_high[j_col]])   
    } }
  C_m }

### fcn RSV age groups --------------
fun_rsv_agegroups<-function(standard_age_groups,rsv_age_groups_low,rsv_age_group_sizes){
  rsv_age_groups=data.frame(age_low=rsv_age_groups_low,age_high=rsv_age_groups_low+rsv_age_group_sizes)
  truthvals=which(match(rsv_age_groups$age_low,standard_age_groups$age_low)==match(rsv_age_groups$age_high,standard_age_groups$age_high))
  rsv_age_groups[,c("wpp_agegroup_low","wpp_agegroup_high")]=NA
  rsv_age_groups[,c("wpp_agegroup_low","wpp_agegroup_high")]=data.frame(t(sapply(1:length(rsv_age_groups_low), function(x) 
  {c(max(which(rsv_age_groups_low[x]>=standard_age_groups$age_low)),
   max(which(rsv_age_groups_low[x]+rsv_age_group_sizes[x]>=standard_age_groups$age_low)))})))
  agelim_diffs=rsv_age_groups$age_high-rsv_age_groups$age_low; agelim_diffs_increm=rep(NA,length(agelim_diffs))
  agelim_diffs_increm[agelim_diffs %% 1==0]=1; agelim_diffs_increm[agelim_diffs %% 1>0]=0.1 
  rsv_age_groups[,"duration"]=(rsv_age_groups$age_high-rsv_age_groups$age_low) + agelim_diffs_increm
  scaling_fact=rsv_age_groups$duration/sapply(1:nrow(rsv_age_groups),function(x) {sum(standard_age_groups$duration[
    rsv_age_groups$wpp_agegroup_low[x]:rsv_age_groups$wpp_agegroup_high[x]])})
  popul_custom_agegroups=sapply(1:nrow(rsv_age_groups),function(x) {sum(standard_age_groups$values[
    rsv_age_groups$wpp_agegroup_low[x]:rsv_age_groups$wpp_agegroup_high[x]])}); rsv_age_groups[,"value"]=NA
# rsv_age_groups$value[truthvals]=standard_age_groups$value[match(rsv_age_groups$age_low,standard_age_groups$age_low)[truthvals]]
  rsv_age_groups$value=popul_custom_agegroups*scaling_fact; rsv_age_groups[,"fraction"]=round(rsv_age_groups$value/sum(rsv_age_groups$value),4)
  rsv_age_groups[,"agegroup_name"]=paste(rsv_age_groups$age_low,rsv_age_groups$age_low+rsv_age_groups$duration,sep='-'); rsv_age_groups
}

### process simul output -----------------
fun_process_simul_output=function(ode_solution,varname_list,n_age,n_inf,rsv_age_groups,neg_thresh){
df_ode_solution=ode_solution %>% as.data.frame() %>% setNames(c("t",fun_sirs_varnames(varname_list,n_age,n_inf)))
df_ode_solution = df_ode_solution[1:(nrow(df_ode_solution)-1),] %>% filter(t %% 1 ==0) 
# neg_thresh=-1e-3
if (any(rowSums(df_ode_solution<neg_thresh)>0)){
  print(paste0("negative values in ", sum(rowSums(df_ode_solution<neg_thresh)>0), " rows!"))
  neg_rows=rowSums(df_ode_solution<neg_thresh)>0
  df_ode_solution[neg_rows,!(colnames(df_ode_solution) %in% "t")]=NA
}
# df_ode_solution_nonzero=df_ode_solution[,colSums(df_ode_solution)>0]
df_ode_solution_tidy=df_ode_solution[,colSums(df_ode_solution,na.rm=T)>0] %>% pivot_longer(!t) # ,id.vars='t')
df_ode_solution_tidy[c('compartment','infection','agegroup')]=
  sapply(1:3, function(x) {sapply(strsplit(as.character(df_ode_solution_tidy$name),'_'),'[[',x)}) 
df_ode_solution_tidy = df_ode_solution_tidy # !grepl("R",compartment) & 
df_ode_solution_tidy$compartment=factor(df_ode_solution_tidy$compartment,levels=varname_list)
df_ode_solution_tidy$agegroup=as.numeric(df_ode_solution_tidy$agegroup)
finalvals=df_ode_solution_tidy %>% group_by(agegroup) %>% filter(t==max(t)) %>% summarise(agegroup_sum_popul=sum(value)) #,na.rm = T
if (any(is.na(finalvals$agegroup_sum_popul))) {
  finalvals=df_ode_solution_tidy %>% group_by(agegroup) %>% filter(t==0) %>% summarise(agegroup_sum_popul=sum(value))
} # print(finalvals)
df_ode_solution_tidy[,"value_fract"]=df_ode_solution_tidy$value/finalvals$agegroup_sum_popul[df_ode_solution_tidy$agegroup]
# df_ode_solution_tidy[,"t_years"]=df_ode_solution_tidy$t/365
df_ode_solution_tidy[,"agegroup_name"]=rsv_age_groups$agegroup_name[df_ode_solution_tidy$agegroup]
df_ode_solution_tidy$agegroup_name=factor(df_ode_solution_tidy$agegroup_name,
                                          levels=unique(df_ode_solution_tidy$agegroup_name))
df_ode_solution_tidy$agegroup_name=factor(paste0("age=",df_ode_solution_tidy$agegroup_name,"yr"),
                                          levels=unique(paste0("age=",df_ode_solution_tidy$agegroup_name,"yr")))
df_ode_solution_tidy$infection=paste0("infection #",df_ode_solution_tidy$infection)

list(df_ode_solution,df_ode_solution_tidy) }

### seasonal forcing term ------------------
fun_seas_forc=function(timesteps,peak_day,st_dev_season,forcing_above_baseline){
# peak_day=60; st_dev_season=27; forcing_above_baseline=0.1; 
dist_from_peak=apply(data.frame( abs(timesteps %% 365-peak_day),365-peak_day+(timesteps %% 365) ),1,min)
forcing_vector= 1 + forcing_above_baseline*exp(-0.5*(dist_from_peak/st_dev_season)^2); forcing_vector }

# shutdown term ----------------------
fun_shutdown_seasforc=function(timesteps,elem_time_step,forcing_above_baseline,npi_strength,npi_year,peak_week,season_width, # shutdown_scale,
                               npi_on,npi_off,n_prec,n_sd){
  # seasonal forcing without shutdown
  seas_forcing=fun_seas_forc(timesteps,peak_day=peak_week*7,st_dev_season=season_width*7,forcing_above_baseline); max_time=max(timesteps)
  # shutdown
  seas_lims=t(data.frame(lapply(n_sd, 
                          function(n) {lapply(0:(max_time/365),function(x){x+(peak_week+c(-n*season_width,n*season_width))/52})})))
  rownames(seas_lims)=c(); colnames(seas_lims)=c("on","off"); seas_lims=subset(data.frame(seas_lims),on<=max(timesteps)/365)
  seas_lims[,"season"]=ceiling(seas_lims$on)
  # forcing with NPI vector
  forcing_vector_npi=seas_forcing
  shutdown_start_day=seas_lims$on[seas_lims$season==(npi_year+1)]*365 - npi_on*7
  shutdown_stop_day=seas_lims$off[seas_lims$season==(npi_year+1)]*365 + npi_off*7
  # shutdown: flat
  forcing_vector_npi[round(shutdown_start_day/elem_time_step):round(shutdown_stop_day/elem_time_step)]=1-npi_strength
  if (sum(forcing_vector_npi>seas_forcing)>0){
     forcing_vector_npi[forcing_vector_npi>seas_forcing] <- seas_forcing[forcing_vector_npi>seas_forcing]}
  # first period should not be on-season
  if (forcing_vector_npi[1]>forcing_above_baseline){
    onsetpoint=max(which(timesteps/365 < seas_lims$on[1]-0.1)); forcing_vector_npi[1:onsetpoint]=forcing_vector_npi[onsetpoint]
    seas_forcing[1:onsetpoint]=seas_forcing[onsetpoint]
  }
  list(forcing_vector_npi,timesteps[c(round(shutdown_start_day/elem_time_step),round(shutdown_stop_day/elem_time_step))],seas_forcing,seas_lims)
}

### plot seasonal forcing ----------------------
fcn_plot_seas_forc <- function(timesteps,seas_force,forcing_vector_npi,shutdwn_lims,seas_lims){
  date_forcing=data.frame(date=as.Date(timesteps+as.numeric(as.Date("01-01-01"))),
            year=as.numeric(sapply(strsplit(as.character(yearweek(as.Date(timesteps)))," "),"[[",1))-1969,
            week=yearweek(as.Date(timesteps)) %>% str_replace(".*W",""),forcing_vector_npi) %>% mutate(year_week=paste0(year,"-",week))
  seas_lims_date=data.frame(on=as.Date(seas_lims$on*365 + as.numeric(date_forcing$date[1])),
                            off=as.Date(seas_lims$off*365 + as.numeric(date_forcing$date[1])))
  shutdwn_lims_date=as.Date(shutdwn_lims + as.numeric(date_forcing$date[1]) )
ggplot(date_forcing) + geom_line(aes(x=date,y=forcing_vector_npi,group=1)) +
    geom_vline(data=seas_lims_date,aes(xintercept=on),color="blue",linetype="dashed",size=0.5) + 
    geom_vline(data=seas_lims_date,aes(xintercept=off),color="blue",linetype="dashed",size=0.5) + 
    geom_rect(aes(xmin=shutdwn_lims_date[1],xmax=shutdwn_lims_date[2],ymin=-Inf,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
    scale_x_date(date_breaks="2 months",expand=expansion(0,0)) + scale_y_continuous(breaks=(0:40)/10,expand=expansion(0,0.2)) +
    theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5,size=6),axis.text.y=element_text(size=7)) +
    ylab("seasonal forcing") + ggtitle(paste0("NPI: [",paste0(shutdwn_lims_date,collapse=", "),"]")) 
# "base line:",round(min(seas_force),1),
}

### set up initial condition ------------------
fcn_set_initconds=function(init_set,init_cond_src,input_from_prev_simul,init_seed,seed_vars,filename){
  if (grepl("previous",init_set)){ 
    # print("using init cond from previous simulation")
    initvals_sirs_model=fcn_init_susc_vals(stationary_init=TRUE,from_file_or_output=init_cond_src,simul_output=input_from_prev_simul,
                                           susc_vars_inds,agegr_sizes=rsv_age_groups$value,sim_filepath=filename)} else {
    # INITIAL INFECTION (taking stationary sol should contain [I]s so no need to re-seed it)
    # set up susceptibles from scratch
    initvals_sirs_model=matrix(0,nrow=dim_sys); initvals_sirs_model[sapply(susc_vars_inds,'[[',1)]=rsv_age_groups$value
    # all first infection groups: sapply(inf_vars_inds, '[[',1) | first infection in first age group: inf_vars_inds[[1]][1]
    if (seed_vars=="all") {seed_vars=sapply(inf_vars_inds,'[[',1)} else {seed_vars=inf_vars_inds[[1]][1]}
    initvals_sirs_model[seed_vars]=init_seed
  }
  round(initvals_sirs_model)
}

### initial susceptible populs -----------------
fcn_init_susc_vals<-function(stationary_init,from_file_or_output,simul_output,susc_vars_inds,agegr_sizes,sim_filepath){
initvals_sirs_model=matrix(0,dim_sys,1); # stationary_init=FALSE
if (stationary_init){
  if (grepl("file",from_file_or_output)) {  x=readRDS(sim_filepath); initvals_sirs_model=as.numeric(x[nrow(x),2:ncol(x)]) } else {
    initvals_sirs_model[,1]=as.numeric(simul_output[nrow(simul_output),2:ncol(simul_output)]) } } else {
    # at t=0 entire popul into susceptibles
initvals_sirs_model[sapply(susc_vars_inds,"[[",1)]=agegr_sizes } # rsv_age_groups$value
round(matrix(initvals_sirs_model))
}

# initial-final totals by age group ----------------------
fun_agegroup_init_final_pop<-function(df_ode_solution_tidy){
df_initial_final_totals=cbind(
    df_ode_solution_tidy %>% group_by(agegroup) %>% filter(t==min(t)) %>% summarise(agegroup_sum_popul=sum(value)),
    df_ode_solution_tidy %>% group_by(agegroup) %>% filter(t==max(t)) %>% summarise(agegroup_sum_popul=sum(value)) )
df_initial_final_totals=df_initial_final_totals[,c(1,2,4)]; colnames(df_initial_final_totals)=c("agegroup","initial_pop","final_pop")
df_initial_final_totals[,"fract_final_init"]=df_initial_final_totals$final_pop/df_initial_final_totals$initial_pop
df_initial_final_totals }

### age distrib for seasons ----------------------
fcn_seas_agedistrib <- function(df_ode_sol_cases_sum,max_time,timesteps,seas_case_threshold){
# AUC: AUC <- sum(diff(x[1:n])*rollmean(y[1:n],2))
season_peaks_AUC=df_ode_sol_cases_sum %>% group_by(agegroup,agegroup_name,season) %>% 
  summarise(max_case=max(symptom_cases),auc_case=sum(diff(t*(1/unique(diff(timesteps))))*rollmean(symptom_cases,2))) %>% 
    group_by(season) %>% mutate(max_case_share_season=max_case/sum(max_case),auc_case_share_season=auc_case/sum(auc_case),
    pre_post_shtdn=ifelse(season>=round(max(shutdwn_lims/365))+1, "post-NPI", "pre-NPI")) %>% 
  group_by(season) %>% filter(any(max_case>seas_case_threshold))
# add pre/npi/post
  season_peaks_AUC$pre_post_shtdn[season_peaks_AUC$season==ceiling(max(shutdwn_lims/365))]="NPI"
  season_peaks_AUC$pre_post_shtdn=factor(season_peaks_AUC$pre_post_shtdn,levels=c("pre-NPI","NPI","post-NPI"))
  season_peaks_AUC$agegr_size=unlist(lapply(lapply(strsplit(gsub("yr","",gsub("age=","",season_peaks_AUC$agegroup_name)),"-"),as.numeric),diff))
  season_peaks_AUC$mean_age=unlist(lapply(lapply(strsplit(gsub("yr","",gsub("age=","",season_peaks_AUC$agegroup_name)),"-"),as.numeric),mean))
season_peaks_AUC}

### create symptom fract table ----------------------
fcn_clin_fract_age_exp <- function(agedep_fact,clinfract_expos,rsv_age_groups,n_age,n_inf){
  clin_fract_age_exp=data.frame(agegroup=rep(1:n_age,n_inf),n_inf=sort(rep(1:n_inf,n_age)),
  infection=unlist(lapply(1:n_inf,function(x) {paste0("infection #",rep(x,n_age))})),
  agegroup_name=factor(rep(rsv_age_groups$agegroup_name,n_inf),levels=unique(rsv_age_groups$agegroup_name)),
  fraction=rep(rsv_age_groups$fraction,n_inf),
  sympt_value=matrix(t(sapply(1:n_age, function(x) {clinfract_expos/((agedep_fact^(x-1)))})),ncol=1))
  clin_fract_age_exp
}

### plot clin fract
fcn_plot_clinfract_table=function(clin_fract_age_exp){
  ggplot(clin_fract_age_exp,aes(x=agegroup_name,y=sympt_value,group=infection,color=infection,linetype=infection)) + geom_line(size=2) +
    theme_bw() + standard_theme + xlab("age group (years)") + ylab("clinical fraction") + scale_y_continuous(breaks=(0:10)/10) +
    theme(legend.position="top",axis.text.x=element_text(size=14),axis.text.y=element_text(size=12),legend.text=element_text(size=12))
}

fcn_plot_suscept_table=function(suscept_table){
ggplot(suscept_table,aes(x=agegroup,y=value,group=name,color=name,linetype=name)) + 
  geom_line(size=2) + theme_bw() + standard_theme + xlab("age group (years)") + ylab("susceptibility") + 
  theme(legend.position="top",axis.text.x=element_text(size=14),axis.text.y=element_text(size=12),legend.text=element_text(size=12))
}

### fcn plot agegroup total populs ----------------------
fcn_suscept_agedeptable <- function(rsv_age_groups,delta_susc,n_inf){
suscept_agedep=data.frame(agegroup=factor(rsv_age_groups$agegroup_name,levels=unique(rsv_age_groups$agegroup_name)),
                          t(delta_susc*matrix(rep(rsv_age_groups$value,n_inf),nrow=n_inf,byrow=T))) %>% pivot_longer(!agegroup)
suscept_agedep$name=gsub("X","infection #",suscept_agedep$name); suscept_agedep
}


### fcn plot agegroup total populs ----------------------
fcn_plotagegroup_totals <- function(df_ode_solution_tidy,scale_val){
  ggplot(df_ode_solution_tidy %>% group_by(t,agegroup_name) %>% summarise(agegroup_total=sum(value)),
       aes(x=t/365,y=round(agegroup_total/1e6,4),group=agegroup_name)) + geom_line() + facet_wrap(~agegroup_name,scales=scale_val) +
  scale_y_log10() + theme_bw() + standard_theme + xlab("year") + ylab("million persons")
}

### age distribution by season ----------------------
fcn_agedistrib_calc_season <- function(df_calc,seaslims,selvar,agegr_merge_min,rsv_age_groups){
subset(df_calc,season>seaslims[1] & season<=seaslims[2] & grepl(selvar,name)) %>% #  max_case
  mutate(agegroup_merged=as.character(agegroup_name)) %>% 
  mutate(agegroup_merged=factor(replace(agegroup_merged,as.numeric(agegroup_name)>agegr_merge_min,
                                        paste0("age>",rsv_age_groups$age_low[agegr_merge_min+1],"yr") ))) %>%
  group_by(season,name,agegroup_merged,pre_post_shtdn) %>% summarise(value=sum(value)) %>% 
  group_by(season,name) %>% mutate(value=replace(value,grepl("share",name),value/sum(value))) %>%
  mutate(value=replace(value,value>1,value/1e6)) %>% 
  mutate(name=replace(name,name %in% "max_case_share_season","share of cases at peak")) %>%
  mutate(name=replace(name,name %in% "auc_case_share_season","share of cases (AUC)")) %>%
  mutate(name=replace(name,name %in% "max_case","number of cases at peak (million)")) %>%
  mutate(name=replace(name,name %in% "auc_case","number of cases (AUC, million)"))
}

### fcn plot age distribution ----------------------
fcn_plot_agedistrib_perseas <- function(df_plot,nrowval,yexpval,xlimvals,subtitlestr,captiontxt,textsize){
ggplot(df_plot,aes(x=season,y=value,group=rev(agegroup_merged))) +
  geom_bar(aes(fill=factor(pre_post_shtdn),alpha=agegroup_merged),color="black",stat="identity") + 
  scale_fill_manual(values=gg_color_hue(3)[colorvals]) + facet_wrap(~name,scales="free",nrow=nrowval) + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(size=12,angle=0),legend.position="bottom",legend.text=element_text(size=12),
        legend.title=element_text(size=14),axis.text.y=element_text(size=12),strip.text=element_text(size=13)) + 
  scale_x_continuous(breaks=seq(round(xlimvals)[1],round(xlimvals)[2]+1,1)) + scale_y_continuous(expand=expansion(yexpval)) + 
  geom_text(aes(x=season,y=value,label=paste0(gsub("age=","",agegroup_merged),"\n",round(value,2))),size=textsize,
            position=position_stack(vjust=0.5),check_overlap=T) +labs(fill="",alpha="",subtitle=subtitlestr,caption=captiontxt) + 
  guides(linetype=guide_legend(override.aes=list(fill=c(NA,NA,NA)))) + xlab("epi-year") + ylab("") + coord_flip()
}

### fcn calc mean age season ----------------------
fcn_calc_mean_age <- function(df_calc,season_min){
subset(df_calc,grepl("share",name)) %>% group_by(season,name) %>%
  summarise(under2yr=sum(value[agegroup<5]*mean_age[agegroup<5]/sum(value[agegroup<5])),
            under3yr=sum(value[agegroup<6]*mean_age[agegroup<6]/sum(value[agegroup<6])),
            under5yr=sum(value[agegroup<7]*mean_age[agegroup<7]/sum(value[agegroup<7])),all_agegroups=sum(value*mean_age/sum(value)),
            pre_post_shtdn=unique(pre_post_shtdn)) %>% pivot_longer(cols=!c(season,name,pre_post_shtdn),names_to="mean_age_type") %>%
  mutate(mean_age_type=factor(mean_age_type,unique(mean_age_type))) %>% group_by(name,mean_age_type) %>% 
  mutate(mean_val=mean(value[pre_post_shtdn=="pre-NPI" & season>season_min-1])) %>% 
  mutate(norm_val=value/mean_val,dep=gsub("exp","previous exposure",gsub("_dep","",gsub("suscept_","",foldername))))
}

### fcn plot mean age season ----------------------
fcn_plot_mean_age <- function(df_plot,seaslims,npiyear,highl_lims,yexpand,ylabtag){
ggplot(subset(left_join(df_plot,seaslims %>% mutate(season=season+1)),
              season>npiyear-1 & season<npiyear+6 & grepl("auc",name) & as.numeric(mean_age_type)>2) %>% mutate(year=season-1)) + 
  geom_segment(aes(x=on,xend=on+1-0.05,y=get(var_sel),yend=get(var_sel),color=dep),size=2) + facet_wrap(~mean_age_type,scales="free") + 
  geom_rect(aes(xmin=highl_lims[1]/365,xmax=highl_lims[2]/365,ymin=-Inf,ymax=Inf),fill="pink",color=NA,alpha=0.1,show.legend=TRUE) +
  geom_vline(data=subset(seaslims %>% pivot_longer(cols=!season),value>npiyear-2 & value<npiyear+4 & name %in% "on"),
             aes(xintercept=value),linetype="dashed",size=0.3) +
  scale_x_continuous(breaks=0:ceiling(xval_lims[2]),expand=expansion(0,0)) + scale_y_continuous(expand=expansion(yexpand[1],yexpand[2])) + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(size=12,vjust=0.75),axis.text.y=element_text(size=12), 
                                      legend.title=element_text(size=12),legend.text=element_text(size=12)) + 
  xlab("epi-year") +ylab(paste0("<age symptomatic cases>",ylabtag)) + labs(color="dependence",fill="NPI",caption=caption_txt)
}

### fun timecourse plot tags ----------------------
fun_tcourse_plottags <- function(k,nval,rval,n_inf,n_age,colvar,agegr_lim,delta_susc_prop,delta_primary,
                                 npi_on,npi_off,shutdown_scale,forcing_above_baseline){

all_perms=permutations(n=nval,r=rval,repeats.allowed=T)
g(k1,k2,k3) %=% all_perms[k,] # assign multiple variables
scale_val=c('fixed','free_y')[k1]; facet2tag=c('','infection')[k2]; value_type=c("value","value_fract")[k3]
# set tags for filename
if (grepl("fract",value_type)){y_axis_tag='fraction'} else {y_axis_tag="# cases"}
# if (nchar(facet2tag)){nrow_val=6} else{nrow_val=3}; 
if (nchar(facet2tag)>0) {height_div=1} else{height_div=2}
# no. of cols
if (grepl("age",colvar)){
ncol_val=agegr_lim # round(agegr_lim/as.numeric(height_div))
facet_formula=paste('~',gsub("^\\+","",paste(facet2tag,'+agegroup_name',sep='')),sep='')} else {
  ncol_val=n_inf # round(agegr_lim/as.numeric(height_div))
  facet_formula=paste('~',paste("agegroup_name+",gsub("^\\+","",facet2tag),sep=''),sep='') }
# name of folder: susceptibility
if (length(unique(round(delta_susc_prop[,1],4)))==1) { 
  if (length(unique(round(delta_susc_prop[1,],4)))>1) { foldername="suscept_age_dep"; subtitle_str="suscept~f(age)"} else {
    foldername="suscept_const"; subtitle_str="suscept~const"}
} else {
  if (length(unique(round(delta_susc_prop[1,],4)))>1) { foldername="suscept_ageexp_dep"; subtitle_str="suscept~f(age,expos)"} else {
    foldername="suscept_exp_dep"; subtitle_str="suscept~f(expos)"}
}
# if (grepl("_noagedep",foldername)) {subtitle_str="suscept~f(expos)"} else {subtitle_str="suscept~f(age,expos)"}
# caption 
caption_txt<-paste0('NPI (',npi_on,',',npi_off,") weeks from CI95 season, -",round(1e2*shutdown_scale),
    "% baseline, ",1e2*round(1/(1+forcing_above_baseline),2),"% off-season activity")
# filename
timecourse_filename=fun_create_filename(paste0("simul_output/",foldername),facet2tag,value_type,
                                        n_age,gsub("_y","",scale_val),"png")
# add npi timing to filename
timecourse_filename=gsub('.png',paste0('_npi',npi_year,'y_on',npi_on,"w_off",npi_off,
                                       'w_maxforc',forcing_above_baseline*1e2,'pct.png'),timecourse_filename)

c(scale_val,facet2tag,value_type,y_axis_tag,ncol_val,facet_formula,foldername,caption_txt,subtitle_str,timecourse_filename)
}

## create table of symptomatic cases ----------------------
fun_symptomcases_table <- function(df_ode_solution_tidy,df_clin_fract,bindcolnames,seas_lims){
  df_ode_solution_tidy_cases=left_join(df_ode_solution_tidy[grepl('I_',df_ode_solution_tidy$name),],
          df_clin_fract[,!grepl("agegroup_name",colnames(df_clin_fract))],by=bindcolnames) %>%
      mutate(symptom_cases=value*sympt_value,symptom_cases_fract=value_fract*sympt_value)
# plot sum of 1,2,3rd infections
df_ode_sol_cases_sum=df_ode_solution_tidy_cases %>% group_by(t,compartment,agegroup,agegroup_name) %>% 
  summarise(symptom_cases=round(sum(symptom_cases)),symptom_cases_fract=sum(symptom_cases_fract)) %>%
  mutate(season=findInterval(t/365,c(0,seas_lims$on)))
df_ode_sol_cases_sum
}

### fun sumcase plot tags
fun_sumcase_plot_tags <- function(n_val,r_val,k_plot,df_symptom_prop,delta_susc_prop,npi_on,
                                  npi_off,shutdown_scale,forcing_above_baseline){
  plot_perms=permutations(n=n_val,r=r_val,repeats.allowed=T)
scale_val=c("free","fixed")[plot_perms[k_plot,1]]; plotvar=c("symptom_cases_fract","symptom_cases")[plot_perms[k_plot,2]]
if (grepl("fract",c("symptom_cases_fract","symptom_cases")[plot_perms[k_plot,2]])){
  y_axis_tag='fraction'} else {y_axis_tag="# cases"}
symptvals_by_age=length(unique(round(subset(df_symptom_prop,n_inf==1)$sympt_value,4)))
symptvals_by_exp=length(unique(round(subset(df_symptom_prop,agegroup==1)$sympt_value,4)))

if (symptvals_by_age>1 & symptvals_by_exp>1){dep_tag="_age_exp_dep"; subtitle_str_exp="clin.fract~f(age,expos)"} else {
  if (symptvals_by_age>1 & symptvals_by_exp==1) {
    dep_tag="_age_dep"; subtitle_str_exp="clin.fract~f(age)"} else if (symptvals_by_age==1 & symptvals_by_exp>1) {
    dep_tag="_exp_dep"; subtitle_str_exp="clin.fract~f(expos)"} else {dep_tag="_const"; subtitle_str_exp="clin.fract~const"}}

if (length(unique(round(delta_susc_prop[,1],4)))==1) { # ncol(matrix(apply(delta_susc_prop,2,unique))
  if (length(unique(round(delta_susc_prop[1,],4)))>1) { foldername="suscept_age_dep"; subtitle_str_susc="suscept~f(age)"} else {
    foldername="suscept_const"; subtitle_str_susc="suscept~const"}
} else {
  if (length(unique(round(delta_susc_prop[1,],4)))>1) { foldername="suscept_ageexp_dep"; subtitle_str_susc="suscept~f(age,expos)"} else {
    foldername="suscept_exp_dep"; subtitle_str_susc="suscept~f(expos)"}
}

full_filename=paste0("simul_output/",foldername,paste0("/sever",dep_tag),"/symptomcases_y",
                     scale_val,'_',gsub("# cases","absval",y_axis_tag),".png")
full_filename=gsub('.png',paste0('_npi',npi_year,'y_on',npi_on,"w_off",npi_off,'w_basrate',
                                 forcing_above_baseline*1e2,'pct.png'),full_filename)
# gsub('.png',paste0('_npi',npi_year,'y_week',npi_on,'.png'),full_filename)

caption_txt<-paste0('NPI (',npi_on,',',npi_off,") weeks from CI95 season, -",round(1e2*shutdown_scale),
                    "% baseline, ",1e2*round(1/(1+forcing_above_baseline),2),"% off-season activity")
if (grepl("_expdep_only",strsplit(full_filename,"/")[[1]][3])) {sever_str="clin.fract.~f(expos)"} else {
  if (grepl("_agedep_only",strsplit(full_filename,"/")[[1]][3])) {sever_str="clin.fract.~f(age)"} else {
    sever_str="clin.fract.~f(age,expos)"}
}

c(scale_val,y_axis_tag,plotvar,foldername,full_filename,caption_txt,paste0(subtitle_str_susc,", ", subtitle_str_exp))
}

# create file name ----------------------
fun_create_filename=function(foldername,facet2tag,value_type,n_age,scale_val,filetype){
if (nchar(facet2tag)==0) {overlay_tag='_overlaid'} else {overlay_tag=''}
if (grepl('fract',value_type)) {fract_abs="fractional"} else {fract_abs="absval"}
timecourse_filename=paste0(foldername,"/",fract_abs,'_y',scale_val,overlay_tag,".",filetype)
timecourse_filename
}

### fcn age distrib plot tags -----------------
fcn_agedistrib_plot_tags <- function(delta_susc_prop,delta_primary,plot_season_peaks,df_symptom_prop,
                                     npi_on,npi_off,npi_strength,seas_case_threshold){
if (length(unique(round(delta_susc_prop[,1],4)))==1) { # ncol(matrix(apply(delta_susc_prop,2,unique))
    if (length(unique(round(delta_susc_prop[1,],4)))>1) { foldername="suscept_age_dep"; subtitle_str_susc="suscept~f(age)"} else {
      foldername="suscept_const"; subtitle_str_susc="suscept~const"}
} else {
    if (length(unique(round(delta_susc_prop[1,],4)))>1) { foldername="suscept_ageexp_dep"; subtitle_str_susc="suscept~f(age,expos)"} else {
      foldername="suscept_exp_dep"; subtitle_str_susc="suscept~f(expos)"}
}
  
symptvals_by_age=length(unique(subset(df_symptom_prop,n_inf==1)$sympt_value))
symptvals_by_exp=length(unique(subset(df_symptom_prop,agegroup==1)$sympt_value))

if (symptvals_by_age>1 & symptvals_by_exp>1){dep_tag="_age_exp_dep"; subtitle_str_exp="clin.fract~f(age)"} else {
  if (symptvals_by_age>1 & symptvals_by_exp==1) {
      dep_tag="_age_dep"; subtitle_str_exp="clin.fract~f(age)"} else if (symptvals_by_age==1 & symptvals_by_exp>1) {
        dep_tag="_exp_dep"; subtitle_str_exp="clin.fract~f(expos)"} else {dep_tag="_const"; subtitle_str_exp="clin.fract~const"} }

agedistr_filename=paste0("simul_output/",foldername,"/sever",dep_tag,"/age_distrib_npi_on",npi_on,"w_off",
                         npi_off,"_NPIred",1e2*npi_strength,"pct.png") # _seasmincase,seas_case_threshold/1e3,"e3
#,"_seas", paste0(unique(plot_season_peaks$season),collapse="_")
agedistr_filename
}

# PLOT time course normalised by max per age group  --------------------------------------------------------
fcn_plot_allcases_absval_stackedinfs <- function(df_plot,valuetype,x_lims,t_subset,agegrlim,ncolval,yaxistag,
                                                 scaleval,vertline_x,highl_lims,xvalbreaks,subtitle_str,caption_txt){
ggplot(subset(df_plot,grepl('I',name) & agegroup<=agegrlim & t/365>x_lims[1] & t/365<x_lims[2]& t%%t_subset==0 ), # 
    aes(x=t/365,y=get(valuetype),group=name)) + geom_area(aes(fill=infection),position=position_stack(reverse=T),color="black",size=0.25) +
  facet_wrap(~agegroup_name,ncol=ncolval,scales=scaleval) + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=11,vjust=0.5),
        axis.text.y=element_text(size=12),legend.position="top",legend.title=element_blank(),strip.text=element_text(size=12)) +
  scale_x_continuous(breaks=xvalbreaks,expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.01,0)) +
  geom_rect(aes(xmin=highl_lims[1]/365,xmax=highl_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=vertline_x %>% pivot_longer(cols=!season),aes(xintercept=value),color="blue",linetype="dashed",size=0.3) +
  xlab('years') + ylab(yaxistag) + labs(subtitle=subtitle_str,caption=caption_txt)
}

# PLOT time course normalised by max per age group  --------------------------------------------------------
fcn_plot_norm_max_byage <- function(df_plot,x_lims,vertline_x,highl_lims,xvalbreaks,subtitle_str,caption_txt){
ggplot(subset(df_plot,grepl('I',name) & agegroup<=agegr_lim & t/365>x_lims[1] & t/365<x_lims[2]) %>%
         group_by(t,agegroup) %>% mutate(sum_val=sum(value)) %>% group_by(agegroup) %>% mutate(value_max_norm=value/max(sum_val)),
    aes(x=t/365,y=value_max_norm,group=name)) + geom_area(aes(fill=infection),color="black",size=0.25,position=position_stack(reverse=T)) +
  facet_wrap(~agegroup_name,ncol=2,scales=scale_val) + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=11,vjust=0.5),
  axis.text.y=element_text(size=12),legend.title=element_blank(),strip.text=element_text(size=12)) +
  scale_x_continuous(breaks=xvalbreaks,expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.01,0)) +
  geom_rect(aes(xmin=highl_lims[1]/365,xmax=highl_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=vertline_x %>% pivot_longer(cols=!season),aes(xintercept=value),color="blue",linetype="dashed",size=0.3) + 
  xlab('years') + ylab("% of maximum incidence") + labs(subtitle=subtitle_str,caption=caption_txt)
}

### PLOT time course faceted by infection --------------------------------------------------------
fcn_plot_norm_max_byage <- function(df,x_lims,vertline_x,highl_lims,xvalbreaks,subtitle_str,caption_txt){
# xval_lims=c(npi_year-1.27,npi_year+3.25); seas_lims_plot=subset(seas_lims,on>xval_lims[1] & off<xval_lims[2])
ggplot(subset(df,grepl('I',name) & t/365>x_lims[1] & t/365<x_lims[2]) %>% 
         group_by(t,infection) %>% mutate(sum_val=sum(value)) %>% group_by(infection) %>% mutate(value_max_norm=value/max(sum_val)),
       aes(x=t/365,y=value_max_norm,group=name)) + geom_area(aes(fill=agegroup_name),color="black",size=0.25,position=position_stack(reverse=T)) + 
  facet_wrap(~infection,ncol=1,scales=scale_val) + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=11,vjust=0.5),
        axis.text.y=element_text(size=12),legend.position="top",legend.title=element_blank(),strip.text=element_text(size=12)) +
  scale_x_continuous(breaks=xvalbreaks,expand=expansion(0,0)) + scale_y_continuous(expand=expansion(0.01,0)) +
  geom_rect(aes(xmin=highl_lims[1]/365,xmax=highl_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=vertline_x %>% pivot_longer(cols=!season),aes(xintercept=value),color="blue",linetype="dashed",size=0.3) + 
  xlab('years') + ylab("% of maximum incidence") + labs(subtitle=subtitle_str,caption=caption_txt)
} 

### fcn plot share of infections --------------------------------------------------------
# xval_lims=c(npi_year-2.25,npi_year+3.25)
fcn_plot_share_infs_agefaceted <- function(df,x_lims,vertline_x,highl_lims,xvalbreaks,n_aver,subtitle_str,caption_txt){
ggplot(subset(df,grepl('I',name) & t/365>x_lims[1] & t/365<x_lims[2] & agegroup<9) %>% group_by(t,agegroup) %>% 
  mutate(value_fract_inf=value/sum(value)) %>% group_by(name) %>% 
  mutate(value_fract_inf_smooth=rollmean(value_fract_inf,k=n_aver,align="center",fill=NA)), 
  aes(x=t/365,y=value_fract_inf_smooth,group=name)) + geom_area(aes(fill=infection),color="black",size=0.25,position=position_stack(reverse=T)) + 
  facet_wrap(~agegroup_name,ncol=4,scales=scale_val) + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=11,vjust=0.5),
      axis.text.y=element_text(size=12),legend.position="top",legend.title=element_blank(),strip.text=element_text(size=12)) +
  scale_x_continuous(breaks=xvalbreaks,expand=expansion(0,0)) + scale_y_continuous(expand=expansion(0,0)) +
  geom_rect(aes(xmin=highl_lims[1]/365,xmax=highl_lims[2]/365,ymin=0,ymax=Inf),fill="grey",color=NA,alpha=0.01) +
  geom_vline(data=subset(vertline_x %>% pivot_longer(cols=!season),name %in% "on"),aes(xintercept=value),color="black",linetype="dashed",size=0.3) +
  geom_vline(data=subset(vertline_x %>% pivot_longer(cols=!season),name %in% "off"),aes(xintercept=value),color="grey30",linetype="dashed",size=0.3) +
  xlab('years') + ylab("% of cases at t") + labs(subtitle=subtitle_str,caption=caption_txt)
}

### fcn plot symptom cases, age faceted --------------------------------------------------------
fcn_plot_symptomcases_agefacet <- function(df_plot,x_lims,y_plotvar,yaxistag,scale_val,vertline_x,highl_lims,xvalbreaks,n_col,n_agelim,
                                           subtitle_str,caption_txt){
ggplot(subset(df_plot,t/365>x_lims[1] & t/365<x_lims[2] & agegroup<=n_agelim),aes(x=t/365,y=get(y_plotvar))) + geom_line() +
  facet_wrap(~agegroup_name,scales=scale_val,ncol=n_col) + theme_bw() + standard_theme +
  theme(axis.text.x=element_text(size=12,vjust=0.5),axis.text.y=element_text(size=13),strip.text=element_text(size=15),legend.position='none') +
  scale_x_continuous(breaks=xvalbreaks,minor_breaks=seq(0,max_time/365,by=1/12),expand=expansion(0,0)) + 
  scale_y_continuous(expand=expansion(0.01,0)) + # 
  geom_rect(aes(xmin=highl_lims[1]/365,xmax=highl_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=vertline_x,aes(xintercept=value),color="blue",linetype="dashed",size=0.3) +
  xlab('years') + ylab(yaxistag) + labs(caption=caption_txt,subtitle=subtitle_str)
}

### plot symptom cases stacked as area plot -----------
fcn_plot_symptomcases_agestacked <- function(df_plot,x_lims,y_plotvar,yaxistag,vertline_x,highl_lims,xvalbreaks,subtitle_str,caption_txt){
ggplot(subset(df_plot,t/365>x_lims[1] & t/365<x_lims[2]),aes(x=t/365,y=get(y_plotvar))) + #  & agegroup<=9
  geom_area(aes(fill=agegroup_name),position=position_stack(reverse=T),color="black",size=0.25) + theme_bw() + standard_theme +
  theme(axis.text.x=element_text(size=12,vjust=0.5),axis.text.y=element_text(size=13),strip.text=element_text(size=15)) +
  scale_x_continuous(breaks=xvalbreaks,minor_breaks=seq(0,max_time/365,by=1/12),expand=expansion(0,0)) + 
  scale_y_continuous(expand=expansion(0,0)) +
  geom_rect(aes(xmin=highl_lims[1]/365,xmax=highl_lims[2]/365,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_vline(data=vertline_x,aes(xintercept=value),color="blue",linetype="dashed",size=0.3) +
  xlab('years') + ylab(yaxistag) + labs(fill="",caption=caption_txt,subtitle=subtitle_str)
}

### bar plot age distrib -----------------
# ggplot(subset(mean_age_perseason,season>3 & grepl("max",name)), #  & mean_age_type %in% "mean_age_under5"
#   aes(x=season,y=value,fill=factor(pre_post_shtdn),group=mean_age_type,alpha=factor(mean_age_type))) +
#   geom_bar(stat="identity",color="black",position="dodge") + scale_alpha_manual(values=1-(0:(length(unique(mean_age_perseason$mean_age_type))-1))*0.2) +
#   scale_x_continuous(breaks=0:max(mean_age_perseason$season)) + scale_y_continuous(breaks=0:max(mean_age_perseason$value),expand=expansion(c(0,0.1))) +
#   scale_fill_manual(values=gg_color_hue(3)[colorvals],guide='none') + theme_bw() + standard_theme +
#   theme(axis.text.x=element_text(size=9,vjust=0.5),axis.text.y=element_text(size=8),legend.position="bottom",legend.text=element_text(size=12)) +
#         labs(alpha="",subtitle=subtitle_str,caption=caption_txt) + ylab("mean age of symptomatic cases") +
#   geom_text(aes(label=round(value,2)),size=4,color="black",alpha=1,vjust=-0.5,position=position_dodge(width=1)) # + coord_flip()
# # SAVE
# ggsave(gsub("age_distrib","mean_age",agedistr_filename),width=22,height=14,units="cm")

### assign multiple variables -----------------
# use as: g(a,b,c) %=% c(1,2,3)
'%=%' = function(l, r, ...) UseMethod('%=%')
# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}
# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

# scale_y_log10(limits=c(0.1,ymaxval)) + scale_size_manual(values=rev(as.numeric(unique(df_ode_solution_tidy$infection)))*0.3)
# facet_grid(compartment~infection,scales='free',labeller=label_both) + # 1 panel: 1 vartype, 1 # infection, 2 age groups
# facet_grid(compartment~agegroup+infection,scales='free',labeller=label_both) + 
# facet_wrap(~compartment+infection,ncol=3,scales='free',labeller=label_both) + # 1 panel: 1 vartype, 1 # inf, 2 age groups  
# facet_wrap(~compartment+agegroup,ncol=2,scales='free') + # 1 panel: 1 vartype, 1 age groups, 3 #s infection, 
# scale_x_log10(limits=c(1,t_maxval),breaks=scales::trans_breaks("log10", function(x) 10^x),
# labels=scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(sides='b') + 