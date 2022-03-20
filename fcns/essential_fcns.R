# essential functions needed to run individual simulations

### assign multiple variables ----------------------------------------------------------
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


# linear index from 2-dim index (infection-age) ----------------------------------------------------------
fun_sub2ind=function(i_inf,j_age,varname,varname_list,n_age,n_inf){
  varnum=which(varname_list %in% varname); k=(j_age-1)*length(varname_list)*n_inf + (varnum-1)*n_inf + i_inf; k }

# generate all model names ----------------------------------------------------------
fun_sirs_varnames=function(varname_list,n_age,n_inf){
  array(sapply(1:n_age, function(x_age) {sapply(varname_list, function(x) {paste0(paste0(x,'_',1:n_inf),'_',x_age)})} ))
}

# set up kinetic matrix ----------------------------------------------------------
fun_K_m_sirs_multiage=function(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,agegroup_durations){
  K_m=matrix(0,nrow=dim_sys,ncol=dim_sys)
  # S_i_j -> S_1_1 is S, subscript=1, superscript=1. subscript: # infection, superscript= # age group
  # varname_list=c('S','I','R') 
  # conversion between i,j and X_k, when variables are stacked as S_i_1,I_i_1,R_i_1, S_i_2,I_i_2,R_i_2 ...
  # waning terms: omega=1/1e2. flow: R_i_j -> S_min(i+1,n_inf)_j;  eg. R_1_1 -> S_2_1
  for (j_age in 1:n_age) {
    for (i_inf in 1:n_inf) { if (j_age==1 & i_inf==1) {waning_terms_source_target=data.frame()}
      wanevals=c(fun_sub2ind(i_inf,j_age,'R',varname_list,n_age,n_inf),
                 fun_sub2ind(min(i_inf+1,n_inf),j_age,'S',varname_list,n_age,n_inf))
      # waning_terms_source_target=rbind(waning_terms_source_target,wanevals)
      K_m[wanevals[2],wanevals[1]]=omega } }
  # aging terms between AGE GROUPS: S_i_j -> S_i_(j+1), R_i_j -> R_i_(j+1), I_i_j->I_(i+1)_j
  for (j_age in 1:(n_age-1)) {
    for (i_inf in 1:n_inf) { if (j_age==1 & i_inf==1) {aging_terms_source_target=data.frame()}
      agevals=rbind( 
        # S_i_j -> S_i_(j+1)
        c(fun_sub2ind(i_inf,j_age,'S',varname_list,n_age,n_inf),
          fun_sub2ind(i_inf,j_age+1,'S',varname_list,n_age,n_inf),j_age ),
        # R_i_j -> R_i_(j+1)
        c(fun_sub2ind(i_inf,j_age,'R',varname_list,n_age,n_inf),
          fun_sub2ind(i_inf,j_age+1,'R',varname_list,n_age,n_inf),j_age),
        # I_i_j -> I_i_(j+1)
        c(fun_sub2ind(i_inf,j_age,'I',varname_list,n_age,n_inf),
          fun_sub2ind(i_inf,j_age+1,'I',varname_list,n_age,n_inf),j_age)
      )
      # bundle
      aging_terms_source_target=rbind(aging_terms_source_target,agevals) } }
  for (k in 1:nrow(aging_terms_source_target)) {
    duration_scale=agegroup_durations[aging_terms_source_target[k,3]]
    # print(c(aging_terms_source_target[k,2],aging_terms_source_target[k,1])) # ,duration_scale
    K_m[aging_terms_source_target[k,2],aging_terms_source_target[k,1]]=1/(365*duration_scale)}
  # print(aging_terms_source_target) # [,1:2]
  # recovery terms
  # rho=1/6; # 1/rho=rweibull(1, shape=4.1,scale=8.3)
  for (j_age in 1:n_age) {
    for (i_inf in 1:n_inf) { if (j_age==1 & i_inf==1) {recov_terms_source_target=data.frame()}
      recov_vals=c(fun_sub2ind(i_inf,j_age,'I',varname_list,n_age,n_inf),
                   fun_sub2ind(i_inf,j_age,'R',varname_list,n_age,n_inf))
      recov_terms_source_target=rbind(recov_terms_source_target,recov_vals)
      K_m[recov_vals[2],recov_vals[1]]=rho } }
  
  # diagonal terms balancing the outgoing terms, these are the (sums of the off diagonal terms)x(-1) 
  diag(K_m)=diag(K_m)-colSums(K_m-diag(diag(K_m)))
  
  #### return output
  # print(diag(K_m))
  K_m
}

###

# age structure of country ----------------------------------------------------------
fun_cntr_agestr <- function(i_cntr,i_year,age_low_vals,age_high_vals){
  age_groups=data.frame(age_low=seq(0,75,5), age_high=c(seq(4,74,5),100))
  if (!any((.packages()) %in% "wpp2019")) {library(wpp2019)}; if (!exists("popF")) {data("pop")}
  cntr_agestr=data.frame(agegroups=popF[popF$name %in% i_cntr,"age"],values=popF[popF$name %in% i_cntr,i_year] +
                           popM[popM$name %in% i_cntr,i_year])
  agegr_truthvals=sapply(strsplit(as.character(cntr_agestr$agegroups),"-"),"[[",1) %in% age_groups$age_low
  N_tot=cntr_agestr$values[agegr_truthvals]
  N_tot[length(N_tot)]=N_tot[length(N_tot)]+sum(cntr_agestr$values[!agegr_truthvals])
  N_tot=N_tot*1e3; # N_tot
  data.frame(age_low=age_low_vals, age_high=age_high_vals,values=N_tot, duration=(age_high_vals-age_low_vals)+1) %>%
    mutate(proportion=values/sum(values))
}


### country full popul struct ----------------------------------------------------------
fcn_cntr_fullpop <- function(n_year,country_sel){
  uk_popul=left_join(subset(popF,name %in% country_sel)[,c("age",n_year)],
                     subset(popM,name %in% country_sel)[,c("age",n_year)],by="age",suffix=c("F","M"))
  uk_popul[,"totalpop"]=uk_popul[,2]+uk_popul[,3]
  uk_popul=uk_popul %>% mutate(lower=as.numeric(gsub("-\\d+","",age)),upper=as.numeric(gsub("\\d+-","",age))+0.9)
  if (any(is.na(uk_popul$lower))){
    uk_popul[grepl("95",uk_popul$age),c(paste0(n_year,"F"),paste0(n_year,"M"),"totalpop")]=
      (uk_popul[grepl("95",uk_popul$age),c(paste0(n_year,"F"),paste0(n_year,"M"),"totalpop")]+
         uk_popul[uk_popul$age=="100+",c(paste0(n_year,"F"),paste0(n_year,"M"),"totalpop")])
    uk_popul=uk_popul[-which(uk_popul$age=="100+"),] %>% 
      mutate(mean_age=(lower+upper+0.1)/2,fraction_pop=totalpop/sum(totalpop)) }
  uk_popul
}

### fcn RSV age groups ----------------------------------------------------------
fun_rsv_agegroups<-function(standard_age_groups,popul_struct,rsv_age_groups_low,rsv_age_group_sizes){
  rsv_age_groups=data.frame(age_low=rsv_age_groups_low,age_high=rsv_age_groups_low+rsv_age_group_sizes)
  truthvals=which(match(rsv_age_groups$age_low,standard_age_groups$age_low)==
                    match(rsv_age_groups$age_high,standard_age_groups$age_high))
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
  rsv_age_groups$value=popul_custom_agegroups*scaling_fact; rsv_age_groups[,"fraction"]=
    round(rsv_age_groups$value/sum(rsv_age_groups$value),4)
  rsv_age_groups[,"agegroup_name"]=paste(rsv_age_groups$age_low,rsv_age_groups$age_low+rsv_age_groups$duration,sep='-'); 
  rsv_age_groups
  
  agegroup_match=data.frame(model_agegroup=1:nrow(rsv_age_groups),
                            age_low=rsv_age_groups$age_low,age_high=rsv_age_groups$age_high,
                            wpp_agegroup_low=unlist(lapply(rsv_age_groups$age_low,
                                          function(x){which(x>=popul_struct$lower & x<=popul_struct$upper)})),
                            wpp_agegroup_high=unlist(lapply(rsv_age_groups$age_high,
                                          function(x){which(x>=popul_struct$lower & x<=popul_struct$upper)}))) %>%
    mutate(age_high=ifelse(model_agegroup<max(model_agegroup),age_low[model_agegroup+1],age_high),
           mean_age_arithm=(age_low+age_high)/2, mean_age_weighted=sapply(1:max(model_agegroup),function(x) {
             sum((popul_struct$totalpop[wpp_agegroup_low[x]:wpp_agegroup_high[x]]*
                    popul_struct$mean_age[wpp_agegroup_low[x]:wpp_agegroup_high[x]])/
                   sum(popul_struct$totalpop[wpp_agegroup_low[x]:wpp_agegroup_high[x]]))})) %>%
    mutate(mean_age_weighted=ifelse(wpp_agegroup_low==wpp_agegroup_high,mean_age_arithm,mean_age_weighted)) %>% 
    select(-mean_age_arithm)
  rsv_age_groups %>% mutate(mean_age_weighted=agegroup_match$mean_age_weighted)
}

### calculate stationary population age structure -----------------
fcn_calc_stat_popul <- function(rsv_agegroups,agegroup_time_size,dailybirth,dailydeath_per_agegrp,output_type){
  # dailybirth=2314
  # all rates PER DAY
  D_m = -diag(dailydeath_per_agegrp) # # uk_death_rate$deaths_per_1000person_peryear/(365*1e3)
  B_m = matrix(0,nrow = nrow(rsv_agegroups),ncol=nrow(rsv_agegroups))
  diag(B_m) = -c(1/(agegroup_time_size[1:(nrow(rsv_agegroups)-1)]*365),0)
  # subdiag terms
  subdiag_inds=(nrow(rsv_agegroups)+1)*(0:(nrow(rsv_agegroups)-2))+2
  B_m[subdiag_inds]=1/(agegroup_time_size[1:(nrow(rsv_agegroups)-1)]*365)
  K_matr=D_m+B_m
  # has almost zero determinant
  if (output_type=="matrix") {-round(inv(K_matr) %*% matrix(c(dailybirth,rep(0,nrow(rsv_agegroups)-1))) )} else {
    c(-round(inv(K_matr) %*% matrix(c(dailybirth,rep(0,nrow(rsv_agegroups)-1))) ))}
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

# create reduced contact matrix --------------
fun_create_red_C_m=function(C_m_full,rsv_agegroups,orig_age_groups_duration,orig_age_groups_sizes){
  C_m=matrix(0,nrow=nrow(rsv_agegroups),ncol=nrow(rsv_agegroups)); rownames(C_m)=rsv_agegroups$agegroup_name
  colnames(C_m)=rsv_agegroups$agegroup_name
  for (i_row in 1:n_age){
    for (j_col in 1:n_age){
      # we are merging or splitting age groups, there are 3 possibilities for a new age group:
      # same OR smaller than (ST) OR larger than (LT) the original
      # (it is an *average* of contacts per person, and we have no resolution within age bands)
      #
      # if the 'i' group (C[i,j]) is the *same* as original or *smaller*, this (in itself) does not change the contact rate
      if (rsv_agegroups$wpp_agegroup_low[i_row]==rsv_agegroups$wpp_agegroup_high[i_row]) {
        # 'j' group same or smaller as original
        if (rsv_agegroups$wpp_agegroup_low[j_col]==rsv_agegroups$wpp_agegroup_high[j_col]) {
          f_dur=rsv_agegroups$duration[j_col]/orig_age_groups_duration[rsv_agegroups$wpp_agegroup_high[j_col]]
          C_m[i_row,j_col]=(C_m_full[rsv_agegroups$wpp_agegroup_low[i_row],rsv_agegroups$wpp_agegroup_low[j_col]])*f_dur
        } else { # if 'j' is larger than original group
          group_span=rsv_agegroups$wpp_agegroup_low[j_col]:rsv_agegroups$wpp_agegroup_high[j_col]
          agegroup_weights=orig_age_groups_sizes[group_span]/sum(orig_age_groups_sizes[group_span])
          C_m[i_row,j_col]=sum(agegroup_weights*C_m_full[i_row,group_span])
        } # end of 'i' smaller or same as original
      } else { # if 'i' in C[i,j] is a bigger age band -> weighted average of the contact rates of constituent groups
        group_span=rsv_agegroups$wpp_agegroup_low[i_row]:rsv_agegroups$wpp_agegroup_high[i_row]
        agegroup_weights=orig_age_groups_sizes[group_span]/sum(orig_age_groups_sizes[group_span])
        # if 'j' is same/smaller -> contact rate with original group proportionally divided
        if (rsv_agegroups$wpp_agegroup_low[j_col]==rsv_agegroups$wpp_agegroup_high[j_col]) {
          f_dur=rsv_agegroups$duration[j_col]/orig_age_groups_duration[rsv_agegroups$wpp_agegroup_high[j_col]]
          C_m[i_row,j_col]=
            sum((orig_age_groups_sizes[group_span]/sum(orig_age_groups_sizes[group_span]))*C_m_full[group_span,j_col])*f_dur
        } else {# if 'j' larger -> weighted average of the contact rates of the constituent groups
          # print(c(i_row,j_col)); print(agegroup_weights); print(group_span)
          C_m[i_row,j_col]=sum(rep(agegroup_weights,length(agegroup_weights))*unlist(
            lapply(group_span,function(x) {agegroup_weights*C_m_full[x,group_span]})))
        }
      }
      # C_m[i_row,j_col]=mean(C_m_full[rsv_agegroups$wpp_agegroup_low[i_row]:rsv_agegroups$wpp_agegroup_high[i_row],
      #                                rsv_agegroups$wpp_agegroup_low[j_col]:rsv_agegroups$wpp_agegroup_high[j_col]])   
    } }
  C_m }


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

### R0 calculation ----------------------
R0_calc_SIRS <- function(Cm,deltasusc,rho,n_inf) {
  susc_matrs=lapply(1:n_inf, function(x) {matrix(rep(deltasusc[x,],ncol(deltasusc)),ncol = ncol(deltasusc))})
  ngm=Reduce('+',lapply(susc_matrs, function(x) {x*Cm}))*(1/rho)*(1/n_inf)
  abs(eigen(ngm)$values[1])
}

# shutdown term ----------------------
fun_shutdown_seasforc <- function(npidates,years_pre_post_npi,season_width_wks,init_mt_day,
                                  peak_week,forcing_above_baseline,npireduc_strength){
  start_end=as.Date(c(paste0(year(npidates[1]-years_pre_post_npi[1]*365),"-",init_mt_day),
                      paste0(year(npidates[2]+years_pre_post_npi[2]*365),"-",init_mt_day)))
  # print(start_end)
  simul_dates=seq(start_end[1],start_end[2],by=1)
  forcingvector_npi=fun_seas_forc(yday(simul_dates),peak_day=peak_week*7,
                                  st_dev_season=season_width_wks*7,forcing_above_baseline,day_num=T)
  npi_ind=simul_dates>npidates[1] & simul_dates<npidates[2]
  forcingvector_npi[npi_ind]=forcingvector_npi[npi_ind]*(1-npireduc_strength)
  # n_years=length(forcing_vector_npi)/365; timesteps <- seq(1,length(forcing_vector_npi),by=1)
  list(length(forcingvector_npi)/365,seq(1,length(forcingvector_npi),by=1),start_end,forcingvector_npi) 
}


### plot seasonal forcing ----------------------
fcn_plot_seas_forc <- function(simul_startend,forcingvector_npi,seas_lims_wks,npidates,date_resol){
  df_seas_forc_npi=data.frame(date=seq(simul_startend[1],simul_startend[2],by=1),
                              forcing=forcingvector_npi) %>% mutate(year=year(date),week=week(date)) 
  ggplot(df_seas_forc_npi,aes(x=date,y=forcing)) + geom_line() + xlab("") +
    geom_vline(data=df_seas_forc_npi %>% filter(week %in% seas_lims_wks) %>% group_by(week,year) %>% filter(date==min(date)),
               aes(xintercept=date),linetype="dashed",color="blue",size=1/4) + 
    scale_x_date(date_breaks=date_resol,expand=expansion(0.01,0)) +
    geom_rect(xmin=npidates[1],xmax=npidates[2],ymin=-Inf,ymax=Inf,fill="pink",alpha=0.01) + theme_bw() + standard_theme 
}

### set up initial condition ----------------------
fcn_set_initconds<-function(rsv_agegroup_sizes,init_set,init_cond_src,input_from_prev_simul,init_seed,seed_vars,filename){
  if (grepl("previous",init_set)){     # print("using init cond from previous simulation")
    initvals_sirs_model=fcn_init_susc_vals(stationary_init=TRUE,
                                           from_file_or_output=init_cond_src,simul_output=input_from_prev_simul,
                                        susc_vars_inds,agegr_sizes=rsv_agegroup_sizes,sim_filepath=filename) } else {
  # INITIAL INFECTION (taking stationary sol should contain [I]s so no need to re-seed it)
 # set up susceptibles from scratch
              initvals_sirs_model=matrix(0,nrow=dim_sys*4/3); 
              initvals_sirs_model[sapply(susc_vars_inds,'[[',1)]=rsv_agegroup_sizes
 # all first infection groups: sapply(inf_vars_inds, '[[',1) | first infection in first age group: inf_vars_inds[[1]][1]
          if (seed_vars=="all") {seed_vars=sapply(inf_vars_inds,'[[',1)} else {seed_vars=inf_vars_inds[[1]][1]}
                                             initvals_sirs_model[seed_vars]=init_seed
                                           }
  round(initvals_sirs_model)
}

## model with MATERNAL IMMUNITY (I am only using this!) ---------------------- 
sirs_seasonal_forc_mat_immun <- function(t,X,parms){
  birthrates=parms[[1]][[1]]; deathrates=parms[[1]][[2]]; Km=parms[[2]]; contmatr_row=parms[[3]]; infvars_inds=parms[[4]]; 
  suscvars_inds=parms[[5]]; deltasusc=parms[[6]]; 
  prot_inf_ind=parms[[7]][[1]]; prot_adults_childb=parms[[7]][[2]]; susc_adults_childb=parms[[7]][[3]]
  dimsys=nrow(Km); proport_adult_susc=sum(X[susc_adults_childb])/(sum(X[susc_adults_childb])+sum(X[prot_adults_childb]))
  birthrates[prot_inf_ind,]=(1-proport_adult_susc)*birthrates[1,]
  birthrates[1,]=proport_adult_susc*birthrates[1,]
  # stack I vars
  inf_vars_stacked=do.call(cbind,lapply(infvars_inds, function(x){X[x]}))
  inf_vars_stacked_fullsize=t(matrix(1,1,n_inf)%*%inf_vars_stacked)
  lambda_vect=diag(approx_seas_forc(t)*array(deltasusc)) %*% contmatr_row %*% inf_vars_stacked_fullsize 
  infection_vect=diag(X[unlist(suscvars_inds)])%*%lambda_vect
  F_vect=matrix(0,dimsys,1)
  F_vect[c(unlist(suscvars_inds),unlist(infvars_inds))]=rbind(-infection_vect,infection_vect+approx_introd(t))
  # append infection vector for incidence
  dXdt=birthrates + F_vect + Km %*% X[1:dimsys] - deathrates*X[1:dimsys]; list(rbind(dXdt,infection_vect)) 
}

### process simul output----------------------
fun_process_simul_output <- function(ode_solution,varnamelist,incidvar,incid_only,
                                     init_date,n_age,n_inf,rsvagegroups,neg_thresh){
  # init_date=as.Date(paste0(as.numeric(format(Sys.Date(),"%Y"))-(npiyear+1),"-01-01"))
  df_ode_solution=ode_solution %>% as.data.frame() %>% 
    setNames(c("t",c(fun_sirs_varnames(varnamelist,n_age,n_inf),fun_sirs_varnames(incidvar,n_age,n_inf)) ))
  df_ode_solution = round(df_ode_solution[1:(nrow(df_ode_solution)-1),] %>% filter(t %% 1 ==0))
  # neg_thresh=-1e-3
  if (any(rowSums(df_ode_solution<neg_thresh)>0)){
    print(paste0("negative values in ", sum(rowSums(df_ode_solution<neg_thresh)>0), " rows!"))
    neg_rows=rowSums(df_ode_solution<neg_thresh)>0
    df_ode_solution[neg_rows,!(colnames(df_ode_solution) %in% "t")]=NA
  }
  df_ode_solution_tidy=df_ode_solution[,colSums(df_ode_solution,na.rm=T)>0] %>% pivot_longer(!t)
  df_ode_solution_tidy[c('compartment','infection','agegroup')]=
    sapply(1:3, function(x) {sapply(strsplit(as.character(df_ode_solution_tidy$name),'_'),'[[',x)})
  df_ode_solution_tidy = df_ode_solution_tidy %>% mutate(compartment=factor(compartment,levels=c(varnamelist,incidvar)),
                                                         agegroup=as.numeric(agegroup))
  finalvals=df_ode_solution_tidy %>% filter(t==max(t) & compartment!=incidvar) %>% group_by(agegroup) %>% 
    summarise(agegroup_sum_popul=sum(value))
  if (any(is.na(finalvals$agegroup_sum_popul))) {
    finalvals=df_ode_solution_tidy %>% group_by(agegroup) %>% filter(t==0) %>% summarise(agegroup_sum_popul=sum(value)) 
    }
  df_ode_solution_tidy = df_ode_solution_tidy %>% mutate(value_fract=value/finalvals$agegroup_sum_popul[agegroup]) %>%
    mutate(agegroup_name=paste0(rsvagegroups$agegroup_name[agegroup]),date=init_date+t,
           agegroup_name=factor(agegroup_name,levels=unique(agegroup_name))) %>% 
    group_by(name) %>% mutate(infection=paste0("infection #",infection)) %>% 
    mutate(value=ifelse(compartment==incidvar,c(0,diff(value)),value),
                     value_fract=ifelse(compartment==incidvar,c(0,diff(value_fract)),value_fract))
  if (incid_only) {df_ode_solution_tidy = df_ode_solution_tidy %>% filter(compartment==incidvar)}
  finpop=fun_agegroup_init_final_pop(df_ode_solution,n_age*n_inf*length(varnamelist))
  
  list(finpop,df_ode_solution,df_ode_solution_tidy) 
}

### process incidence only, fast ----------------------
fcn_process_odesol_incid <- function(odesol,n_agegr,n_infect,n_comp,date_start_end){
  odesol[,c(1,(n_agegr*n_infect*n_comp+2):ncol(odesol))] %>% as.data.frame() %>% 
    setNames(c("t",paste0(rep(1:3,n_agegr),"_", unlist(lapply(1:n_agegr, function(x) rep(x,n_infect)))))) %>%
    pivot_longer(!t) %>% group_by(name) %>% 
    mutate(value=value-lag(value,n=1,order_by=t),date=t+date_start_end[1],
                agegroup=as.numeric(sapply(strsplit(name,"_"),"[[",2)),
                infection=as.numeric(sapply(strsplit(name,"_"),"[[",1)))
}

### sub-function to create seasonal forcing term ----------------------
fun_seas_forc <- function(time_input,peak_day,st_dev_season,forcing_above_baseline,day_num){
  # peak_day=60; st_dev_season=27; forcing_above_baseline=0.1; 
  if (!day_num) { time_input = time_input %% 365}
  dist_from_peak=apply(data.frame( abs(time_input - peak_day),365-peak_day+time_input ),1,min)
  forcing_vector= 1 + forcing_above_baseline*exp(-0.5*(dist_from_peak/st_dev_season)^2); forcing_vector 
}

### initial susceptible populs -----------------
fcn_init_susc_vals <- function(stationary_init,from_file_or_output,simul_output,susc_vars_inds,agegr_sizes,sim_filepath){
  if (stationary_init){
    if (grepl("file",from_file_or_output)) { 
      x=readRDS(sim_filepath); initvals_sirs_model=as.numeric(x[nrow(x),2:ncol(x)]) 
    } else { initvals_sirs_model[,1]=as.numeric(simul_output[nrow(simul_output),2:ncol(simul_output)]) } } else {
      # at t=0 entire popul into susceptibles
      initvals_sirs_model=matrix(0,ncol(simul_output)-1,1); # stationary_init=FALSE
      initvals_sirs_model[sapply(susc_vars_inds,"[[",1)]=agegr_sizes 
      }
  round(matrix(initvals_sirs_model)) 
}

# initial-final totals by age group ----------------------
fun_agegroup_init_final_pop <- function(odesol,dimsys){
  finalpop=data.frame(t(odesol[c(1,nrow(odesol)),2:(dimsys+1)])) %>% rownames_to_column(var="name")
  colnames(finalpop)=c("name","init","final"); 
  finalpop <- finalpop %>% separate(name,sep="_",into=c("compartment","inf","agegroup")) %>% group_by(agegroup) %>% 
    summarise(init=sum(init),final=sum(final)) %>% mutate(agegroup=as.numeric(agegroup)) %>% arrange(agegroup)
  finalpop 
}
