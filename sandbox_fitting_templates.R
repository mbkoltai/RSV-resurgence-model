### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### single simulation with unchanging susceptibility --------------------------------------------------------
n_years=15; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=elem_time_step)
# seasonal forcing
basal_rate=0.3; R0_val=round(basal_rate*R0_calc_SIRS(C_m,delta_susc_prop,rho,n_inf),2)
forcing_vector_npi=rep(basal_rate,length(timesteps))

params=list(birth_term,K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,forcing_vector_npi,elem_time_step,delta_susc)
# SIMULATE
ptm<-proc.time(); ode_solution<-lsoda(initvals_sirs_model,timesteps,func=sirs_seasonal_forc,parms=params); proc.time()-ptm
# reshape data
list_simul_output=fun_process_simul_output(ode_solution,varname_list,n_age,n_inf,rsv_age_groups)
df_ode_solution=list_simul_output[[1]]; df_ode_solution_tidy=list_simul_output[[2]]; rm(list_simul_output)
# plot
ggplot(subset(df_ode_solution_tidy,grepl('I',name)&agegroup<=agegr_lim ),
       aes_string(x="t_years",y=value_type,group="name",color="infection")) + geom_line(size=1.025) + theme_bw() + standard_theme +
  theme(axis.text.x=element_text(size=9,vjust=0.5),axis.text.y=element_text(size=9),legend.position="top",legend.title=element_blank()) +
  facet_wrap(as.formula(facet_formula),ncol=as.numeric(ncol_val),scales=scale_val) + # scale_x_continuous(breaks=xval_breaks) +
  ggtitle('RSV infections by age group') + xlab('years') + ylab(y_axis_tag) + labs(caption=paste0("R0=",R0_val)) + xlim(c(0,5))

# save
ggsave(paste0("simul_output/suscept_noagedep/RSV_DE_R0_",R0_val,"_longterm.png"), 
       width=plot_x_y[1],height=plot_x_y[2],units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# simple ODE fcn examples --------------------------------------------------------
#
# linear chain model: -->x1-->x2-->x3 flow, x2 accumulates bc its outflow parameter is small, its value exceeds x1
parvect=c(1e1,1e1,1,1e1); 
linear_chain_ode <- function(t,x,parms) { p0=parms[1];p1=parms[2];p2=parms[3];p3=parms[4]; birth_vect=c(p0,0,0)
kin_matr=matrix(c(-p1,0,0, p1,-p2,0, 0,p2,-p3),3,3,byrow=T); dxdt=birth_vect+kin_matr%*%x; list(dxdt) } #####
ode_solution_toymodel <- lsoda(c(0,0,0),seq(0,5e1,0.1),func=linear_chain_ode,parms=parvect)
df_ode_solution_toymodel=ode_solution_toymodel %>% as.data.frame() %>% setNames(c("t",paste0('x',1:3)))
ggplot(melt(df_ode_solution_toymodel,id.vars='t'),aes(x=t,y=value,color=variable)) + geom_line() + standard_theme + theme_bw()

### ageing process model --------------------------------------------------------
births=10; agegr_sizes=c(1,2,3); g(d1,d2,d3) %=% c(rep(1/n_days_year,3)/agegr_sizes)
K_small=-diag(c(d1,d2,d3)); K_small[2,1]=abs(diag(K_small)[1]); K_small[3,2]=abs(diag(K_small)[2]); inits=c(365*births*agegr_sizes)
ode_aging=function(t,X,parms){births=parms[[1]]; K_small=parms[[2]]; dxdt=matrix(c(births,0,0)) + K_small%*%X; list(dxdt)}
ptm<-proc.time();ode_aging_sol<-lsoda(inits,times=(0:4e4)/10,func=ode_aging,parms=list(births,K_small));proc.time()-ptm
varnames=c('01','12','36'); df_ode_aging_sol=ode_aging_sol %>% as.data.frame() %>% setNames(c("t",varnames))
matplot(df_ode_aging_sol$t,df_ode_aging_sol[,2:ncol(df_ode_aging_sol)],type="l",col=1:3,xlab="days",ylab="# ppl")
legend("topright",legend=varnames,col=1:length(varnames),lty=1:3,cex=0.5,pt.cex=5)

# check if kinet matrix OK (aging/waning terms)
# source_vars=fun_sub2ind(1:3,10,varname = "R",varname_list,n_age,n_inf); target_vars=fun_sub2ind(c(2,3,3),11,"S",varname_list,n_age,n_inf)
source_vars=fun_sub2ind(1:3,10,"S",varname_list,n_age,n_inf); target_vars=fun_sub2ind(1:3,11,"S",varname_list,n_age,n_inf)
for (x in 1:3) {print(K_m[target_vars[x],source_vars[x]])}

# sir model
sir <- function(t,y,parms) {
  beta <- parms[1]; gamma <- parms[2]; S <- y[1]; I <- y[2]
  return(list(c(S = -beta*S*I, I = beta*S*I - gamma*I))) }
# sir with birth/death
sir_birth_death <- function(t,y,parms){
  mu=parms[1]; beta=parms[2]; gamma=parms[3]; S=y[1]; I=y[2]
  dxdt = c(-beta*S*I-mu*S+mu,beta*S*I - (gamma+mu)*I)
  list(dxdt) }
# t - time; y - current state vector of the ODE at time t; parms - Parameter vector used by the ODE system
# Returns list with 1 component being a vector of length two containing: dS(t)/dt and dI(t)/dt
####

# simulate simple ODE examples (fractional) -------------------------------
gamma=1/3; mu=1/60; beta=1.05; sigma=beta/(gamma+mu)
initvals_S_I=c(0.99,0.01); max_time<-150; timesteps <- seq(0, max_time, by=0.1)
# ode_solution <- rk4(y=initvals_S_I,times=timesteps,func=sir_birth_death,parms=c(mu,beta,gamma))
params=c(mu,beta,gamma)
ode_solution <- lsoda(initvals_S_I,timesteps,func=sir_birth_death,parms=params)
df_ode_solution=ode_solution %>% as.data.frame() %>% setNames(c("t", "S", "I"))

# time vs state vars or statevars-statevars
ggplot(df_ode_solution, aes(x=S,y=I,color=t)) + # df_ode_solution
  # geom_line(aes(y=S),color="darkred") + geom_line(aes(y=I), color="steelblue", linetype="twodash") +
  geom_path() + scale_colour_gradient() + # scale_colour_gradientn(colours=terrain.colors(length(timesteps))) + 
  ggtitle('SIR+birth model') + theme_bw() + theme(axis.title=element_text(size=9),
                                                  panel.grid=element_line(linetype="dashed",colour="black",size=0.1)) + 
  # xlab('time') + ylab('S,I') + xlim(0,max(timesteps)) + ylim(0,1)
  xlab('S') + ylab(('I'))

####
# calculate sum of square error ----------------------------------------------------
data=df_ode_solution # 
data[,c('S','I')] = data[,c('S','I')]*(1+matrix(rnorm(2*nrow(df_ode_solution),mean=0,sd=0.03),ncol=2))
# cbind(df_ode_solution$t,matrix(rnorm(2*nrow(df_ode_solution),mean=0,sd=0.01),ncol=2))
# sum of squared error
sir_sse = function(logparams,y){
  df_simul <- lsoda(initvals_S_I,timesteps,func=sir_birth_death,parms=exp(logparams)) %>% 
    as.data.frame() %>% setNames(c("t", "S", "I"))
  I_sqerr = (data$I - df_simul$I)^2; S_sqerr = (data$S - df_simul$S)^2
  sse=sum(I_sqerr) + sum(S_sqerr); sse }

sir_sse(log(params),data)

# Negative log likelihood ----------------------------------------------------
# assuming poisson distrib
N=1e6 # popul size
# calculate negative loglikelihood
sir_nll = function(logparams,y,popul_size){
  df_simul <- lsoda(initvals_S_I,timesteps,func=sir_birth_death,parms=exp(logparams)) %>% 
    as.data.frame() %>% setNames(c("t", "S", "I"))
  logdensities_I=dpois(x=round(popul_size*data$I),lambda=round(popul_size*df_simul$I),log=TRUE)
  logdensities_S=dpois(x=round(popul_size*data$S),lambda=round(popul_size*df_simul$S),log=TRUE)
  nll=-sum(logdensities_I) + (-sum(logdensities_S))
  nll
}

# R squared (least squares)
SS_tot=sum((df_ode_solution[,c('S','I')]-colMeans(df_ode_solution[,c('S','I')]))^2)
SS_res=sir_sse(log(params),data) # sum((data[,c('S','I')]-df_ode_solution[,c('S','I')])^2)
R_sq = 1-SS_res/SS_tot

# NegLogLikelihood
par_initguess=c(0.01,1,0.1)
sir_nll(logparams=log(par_initguess),data,N)/N
# with true value
sir_nll(logparams=log(params),data,N)/N

# plot
ggplot(sq_err_llh) + # geom_point(aes(x=pois_err,y=sq_err)) + 
  geom_point(aes(x=negbinom_err_size1,y=sq_err),color="blue") + geom_point(aes(x=negbinom_err_size10,y=sq_err),color="red") +  
  geom_point(aes(x=negbinom_err_size100,y=sq_err),color="green") +
  theme_bw() + standard_theme + scale_y_log10()
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# least-square fitting with Nelder-Mead ----------------------------------------------------
# true values: params=c(1/60, 1.05,1/3)
par_initguess=c(0.01,1,0.1); parnames=c('mu','beta','gamma')
# Initialize parameters object
optim_proc=capture.output(optim(log(par_initguess),fn=sir_sse,method='Nelder-Mead',y=data,control=c(trace=2)))
####
# extract convergence process
fcn_extract_optimresults <- function(optim_proc,parnames){
  optim_outputs=as.data.frame(
    t(trimws(gsub("\\s+", " ",str_replace_all(optim_proc[which(grepl('\\$',optim_proc))+1],'\\[1\\]','')))),stringsAsFactors=FALSE)
  colnames(optim_outputs)=str_replace_all(optim_proc[grepl('\\$',optim_proc)],'\\$','')
  optim_outputs=optim_outputs %>% separate(par,into=paste0('log(',parnames,')'),sep=' ')
  optim_outputs[,1:length(parnames)]=as.numeric(optim_outputs[,1:length(parnames)])
  optim_outputs
}
####
optim_outputs=fcn_extract_optimresults(optim_proc,parnames)
# convergence process
conv_proc=gsub("\\s+", " ",optim_proc[(which(grepl('Stepsize|Exit',optim_proc))[1]+1):(which(grepl('Stepsize|Exit',optim_proc))[2]-1)])
conv_proc=as.data.frame(conv_proc) %>% separate('conv_proc',c('step_type','n_step','sse1','sse2'),sep=' ')
# post-optim error
sir_sse(as.numeric(optim_outputs[,1:length(parnames)]),data)
# initguess error
sir_sse(log(par_initguess),data)

# with loglikelihood
optim_proc_max_llh=optim(log(par_initguess),fn=sir_nll,method='Nelder-Mead',y=data,popul_size=N,control=c(trace=2))
# yes!
# fractional error
abs(exp(optim_proc_max_llh$par)-params)/params

####
# compare poisson with squared error
singlefit=left_join(subset(df_ode_solution_scan,parset_ID==3000 & !is.na(sum_weekly)), # 3379
                    newburials_weekly[,!grepl("epi_year|week",colnames(newburials_weekly))],by="date",suffix=c("_model","_satel_data"))
adj_factor=mogad_weekly_death_rate/basal_burial_rate
sq_err_llh=data.frame(data=round(singlefit$new_graves_best_ipol-basal_burial_rate)*adj_factor, pred=round(singlefit$sum_weekly),
                      pois_err=-dpois(x=round((singlefit$new_graves_best_ipol-basal_burial_rate)*adj_factor),lambda=round(singlefit$sum_weekly),log=TRUE),
                      sq_err=(round((singlefit$sum_weekly-(singlefit$new_graves_best_ipol-basal_burial_rate))*adj_factor ))^2,
                      negbinom_err_size1=-dnbinom(x=round((singlefit$new_graves_best_ipol-basal_burial_rate)*adj_factor),mu=round(singlefit$sum_weekly),size=1,log=TRUE),
                      negbinom_err_size10=-dnbinom(x=round((singlefit$new_graves_best_ipol-basal_burial_rate)*adj_factor),mu=round(singlefit$sum_weekly),size=10,log=TRUE),
                      negbinom_err_size100=-dnbinom(x=round((singlefit$new_graves_best_ipol-basal_burial_rate)*adj_factor),mu=round(singlefit$sum_weekly),size=100,log=TRUE))

###
# Amadillo sample code ----------------------------------------------------
# (code from Stefan: https://github.com/StefanFlasche/simpleCorePneumoModel)
library("RcppArmadillo"); library("coda"); library("rootSolve") # install.packages("rootSolve")
sourceCpp("toymodel_rcpp/functions.cpp")

df_kilifi=data.frame(Setting="kilifi",
                     Age_groups=c("<1y","1-5y","6-14y","15-20y","21-49y","50+"),
                     Age_group_upper = c(1,6,15,21,50,62),
                     Population=c(9617,45170,68547,33289,72143,24214),
                     VT.prev=c(0.41,0.338,.146,.142,.072,.038),#PCV10
                     NVT.prev=c(.46,.44,.39,.25,.22,.20),
                     N.carr=c(39,127,82,56,97,104),
                     VT.clear=c(.062,.12,.34,.34,.34,.34)/7, #weekly conversion to daily
                     NVT.clear=c(.086,.15,.34,.34,.34,.34)/7)

# set parameters
df=df_kilifi
competition=0.1
no.agegps=dim(df)[1]
agegp.l=df$Age_group_upper - c(0,df$Age_group_upper[-no.agegps])
state_prop=c(agegp.l*100-2,rep(1,no.agegps),rep(1,no.agegps),rep(0,no.agegps))
res = runsteady(y=state_prop, fun=SIS_ode_cpp, parms=param.list, times=c(0,1e5))