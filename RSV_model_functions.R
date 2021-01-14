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
#   forcing_vector=parms[[6]]; elem_time_step=parms[[7]]; # event_vector=parms[[8]]
#   # stack I vars
#   inf_vars_stacked=do.call(cbind,lapply(inf_vars_inds, function(x){X[x]}))
#   inf_vars_stacked_fullsize=t(matrix(1,1,n_inf)%*%inf_vars_stacked)
#   lambda_vect=diag(forcing_vector[(t/elem_time_step)+1]*array(delta_susc))%*%contmatr_rowvector%*%inf_vars_stacked_fullsize
#   infection_vect=diag(X[unlist(susc_vars_inds)])%*%lambda_vect
#   F_vect=matrix(0,dim_sys,1); F_vect[c(unlist(susc_vars_inds),unlist(inf_vars_inds))]=rbind(-infection_vect,infection_vect)
#   dXdt=birth_term + F_vect + K_m%*%X; list(dXdt) }

# set plotting theme
standard_theme=theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
                     plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),
                     axis.text.y=element_text(size=9),
                     axis.title=element_text(size=14), text=element_text(family="Calibri"))

# generate linear index from 2-dim index (infection-age) ----------------------------------------------------------
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

# ode w/o seasonal forcing --------------
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
  data.frame(age_low=age_low_vals, age_high=age_high_vals,values=N_tot, duration=(age_high_vals-age_low_vals)+1)
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

# functn call covidm contact matrix --------------
fun_covidm_contactmatrix <- function(country_sel,currentdir_path,cm_path){if (!exists("covid_params")){
  cm_force_rebuild=F; cm_build_verbose=T; cm_version=2; setwd(cm_path); source(file.path(cm_path,"R","covidm.R"))
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
  rsv_age_groups$value=popul_custom_agegroups*scaling_fact
  rsv_age_groups[,"agegroup_name"]=paste(rsv_age_groups$age_low,rsv_age_groups$age_low+rsv_age_groups$duration,sep='-'); rsv_age_groups
}

# process output
fun_process_simul_output=function(ode_solution,varname_list,n_age,n_inf,rsv_age_groups){
df_ode_solution=ode_solution %>% as.data.frame() %>% setNames(c("t",fun_sirs_varnames(varname_list,n_age,n_inf)))
# df_ode_solution_nonzero=df_ode_solution[,colSums(df_ode_solution)>0]
df_ode_solution_tidy=df_ode_solution[,colSums(df_ode_solution)>0] %>% pivot_longer(!t) # ,id.vars='t')
df_ode_solution_tidy[c('compartment','infection','agegroup')]=
  sapply(1:3, function(x) {sapply(strsplit(as.character(df_ode_solution_tidy$name),'_'),'[[',x)})
df_ode_solution_tidy$compartment=factor(df_ode_solution_tidy$compartment,levels=varname_list)
df_ode_solution_tidy$agegroup=as.numeric(df_ode_solution_tidy$agegroup)
finalvals=df_ode_solution_tidy %>% group_by(agegroup) %>% filter(t==max(t)) %>% summarise(agegroup_sum_popul=sum(value))
df_ode_solution_tidy[,"value_fract"]=df_ode_solution_tidy$value/finalvals$agegroup_sum_popul[df_ode_solution_tidy$agegroup]
df_ode_solution_tidy[,"t_years"]=df_ode_solution_tidy$t/365
df_ode_solution_tidy[,"agegroup_name"]=rsv_age_groups$agegroup_name[df_ode_solution_tidy$agegroup]
df_ode_solution_tidy$agegroup_name=factor(df_ode_solution_tidy$agegroup_name,
                                          levels=unique(df_ode_solution_tidy$agegroup_name))
df_ode_solution_tidy$agegroup_name=factor(paste0("age=",df_ode_solution_tidy$agegroup_name,"yr"),
                                          levels=unique(paste0("age=",df_ode_solution_tidy$agegroup_name,"yr")))
df_ode_solution_tidy$infection=paste0("infection #",df_ode_solution_tidy$infection)

list(df_ode_solution,df_ode_solution_tidy) }

# seasonal forcing term ------------------
fun_seas_forc=function(timesteps,peak_day,st_dev_season,basal_rate){
# peak_day=60; st_dev_season=27; basal_rate=0.1; 
dist_from_peak=apply(data.frame( abs(timesteps %% 365-peak_day),365-peak_day+(timesteps %% 365) ),1,min)
forcing_vector=basal_rate + exp(-0.5*(dist_from_peak/st_dev_season)^2); forcing_vector }

# shutdown term ----------------------
fun_shutdown_seasforc=function(shutdown_start_week,shutdown_stop_week,shutdown_scale,forcing_vector,elem_time_step,basal_rate,n_prec){
shutdown_timewindow=c(shutdown_start_week*7/elem_time_step,shutdown_stop_week*7/elem_time_step)
start_repl=max(which(forcing_vector[1:shutdown_timewindow[1]]-basal_rate<n_prec))
stop_repl=min(which(forcing_vector[shutdown_timewindow[2]:length(forcing_vector)]-basal_rate<n_prec))+shutdown_timewindow[2]
forcing_vector[start_repl:stop_repl]=forcing_vector[start_repl:stop_repl]*shutdown_scale+basal_rate
list(forcing_vector,timesteps[shutdown_timewindow])}

# initial-final totals by age group ----------------------
fun_agegroup_init_final_pop<-function(df_ode_solution_tidy){
df_initial_final_totals=cbind(
    df_ode_solution_tidy %>% group_by(agegroup) %>% filter(t==min(t)) %>% summarise(agegroup_sum_popul=sum(value)),
    df_ode_solution_tidy %>% group_by(agegroup) %>% filter(t==max(t)) %>% summarise(agegroup_sum_popul=sum(value)) )
df_initial_final_totals=df_initial_final_totals[,c(1,2,4)]; colnames(df_initial_final_totals)=c("agegroup","initial_pop","final_pop")
df_initial_final_totals[,"fract_final_init"]=df_initial_final_totals$final_pop/df_initial_final_totals$initial_pop
df_initial_final_totals }

### create symptom fract table ----------------------
fun_propsymptom_table <- function(list_symptom_agegroups,expos_dep_val,agegroupname,n_inf){
for (k in 1:n_inf){
  a=data.frame(t(rbind(unlist(list_symptom_agegroups),
          unlist(sapply(1:length(list_symptom_agegroups), function(x){rep(x,length(list_symptom_agegroups[[x]]))})),
          prop_symptom[unlist(sapply(1:length(list_symptom_agegroups), 
          function(x){rep(x,length(list_symptom_agegroups[[x]]))}))]*(1/(k^expos_agedep)) )),
          infection=paste0("infection #",as.character(k)), n_inf=k )
if (k==1){df_symptom_prop=a} else{df_symptom_prop=rbind(df_symptom_prop,a)} 
}
colnames(df_symptom_prop)[1:n_inf]=c("agegroup","sympt_group","sympt_value"); 
df_symptom_prop[,"agegroup_name"]=factor(agegroupname[df_symptom_prop$agegroup],levels=agegroupname)
df_symptom_prop
}

### fcn plot agegroup total populs ----------------------
fcn_suscept_agedeptable <- function(rsv_age_groups,delta_susc,n_inf){
suscept_agedep=data.frame(agegroup=factor(rsv_age_groups$agegroup_name,levels=unique(rsv_age_groups$agegroup_name)),
                          t(delta_susc*matrix(rep(rsv_age_groups$value,n_inf),nrow=n_inf,byrow=T))) %>% pivot_longer(!agegroup)
suscept_agedep$name=gsub("X","infection #",suscept_agedep$name); suscept_agedep
}

### fcn plot agegroup total populs ----------------------
fcn_plotagegroup_totals <- function(df_ode_solution_tidy,scale_val){
  ggplot(df_ode_solution_tidy %>% group_by(t_years,agegroup_name) %>% summarise(agegroup_total=sum(value)),
       aes(x=t_years,y=agegroup_total,group=agegroup_name)) + geom_line() + facet_wrap(~agegroup_name,scales=scale_val) +
  scale_y_log10() + theme_bw() + standard_theme + xlab("year") + ylab("million popul")
}

### fun timecourse plot tags ----------------------
fun_tcourse_plottags <- function(k,nval,rval){
all_perms=permutations(n=nval,r=rval,repeats.allowed=T)
g(k1,k2,k3) %=% all_perms[k,] # assign multiple variables
scale_val=c('fixed','free_y')[k1]; facet2tag=c('','infection')[k2]; value_type=c("value","value_fract")[k3]
# set tags for filename
if (grepl("fract",value_type)){y_axis_tag='fraction'} else {y_axis_tag="# cases"}
if (nchar(facet2tag)){nrow_val=6} else{nrow_val=3}; if (nchar(facet2tag)>0) {height_div=1} else{height_div=2}
facet_formula=paste('~',gsub("^\\+","",paste(facet2tag,'+agegroup_name',sep='')),sep='')
c(scale_val,facet2tag,value_type,y_axis_tag,nrow_val,height_div,facet_formula)
}

## create table of symptomatic cases ----------------------
fun_symptomcases_table <- function(df_ode_solution_tidy,df_symptom_prop,bindcolnames){
  # bindcolnames=c("infection","agegroup")
  df_ode_solution_tidy_cases=left_join(df_ode_solution_tidy[grepl('I_',df_ode_solution_tidy$name),],
                                     df_symptom_prop,by=bindcolnames)
df_ode_solution_tidy_cases[,"symptom_cases"]=df_ode_solution_tidy_cases$value*df_ode_solution_tidy_cases$sympt_value
df_ode_solution_tidy_cases[,"symptom_cases_fract"]=df_ode_solution_tidy_cases$value_fract*df_ode_solution_tidy_cases$sympt_value
# plot sum of 1,2,3rd infections
# df_ode_solution_tidy_cases_sum=df_ode_solution_tidy_cases %>% group_by(t_years,compartment,agegroup,agegroup_name) %>% 
#   summarise(symptom_cases=sum(symptom_cases),symptom_cases_fract=sum(symptom_cases_fract))
df_ode_solution_tidy_cases
}

# create file name ----------------------
fun_create_filename=function(foldername,facet2tag,value_type,n_age,scale_val,filetype){
if (nchar(facet2tag)==0) {overlay_tag='_overlaid'} else {overlay_tag=''}
if (grepl('fract',value_type)) {fract_abs="fractional"} else {fract_abs="absval"}
timecourse_filename=paste0(foldername,"/RSV_DE_",paste0(n_age,'agegroups_'),fract_abs,'_y',scale_val,overlay_tag,".",filetype)
timecourse_filename
}

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