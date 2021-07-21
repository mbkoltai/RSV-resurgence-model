# simulations with Hodgson model
# may need to install these
library(Rcpp)       # For c++ integration
library(RcppEigen)  # Ditto
library(tidyverse)  # For dataframe manipulation
library(gridExtra)  # For plotting purposes library(reshape2)
file.path = "rsv_model_R-master"
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
load(file=paste0(file.path, '/Rdata/posteriors.Rda')) # The posteriors from the true model
# Load in model data, including observation data (which is made up i.e. not the same as I used in the paper).
load(file=paste0(file.path, '/Rdata/rsv_data.Rda'))   
# Load in the epidemic model from c++ and update important parameters
# Need to have boost installed on computer. Can install via homebrew.
# Might throw loads of [-Wunknown-pragmas] warnings, just ignore, usually an issue with the coompiler.
sourceCpp(paste0(file.path, "/src/logLikelihoodModule.cpp")) # ensure c++14 is enabled
# Calls cpp class
classEvaluateLogLikelihood <- new(EvaluateLogLikelihood, numberDailyLiveBirths, population, ageGroupBoundary) 
classEvaluateLogLikelihood$contactMatrixPhy <- contactMatrixPhy # Physical contact matrix 
classEvaluateLogLikelihood$contactMatrixCon <- contactMatrixCon # Conversational contact matrix 
classEvaluateLogLikelihood$observedData <- as.matrix(observationalData) 
classEvaluateLogLikelihood$lowerParamSupport <- fit.par$lowerParSupport # lower support of parameters 
classEvaluateLogLikelihood$upperParamSupport <- fit.par$upperParSupport # upper support of parameters 
classEvaluateLogLikelihood$run_full <- nrow(observationalData)*7       # number of days to fit the data and model to
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
source(paste0(file.path,"/R/plot_helpers.R"))

output <- as.data.frame(classEvaluateLogLikelihood$getWeeklySampleCpp(as.numeric(post[1,]),FALSE)) %>% 
  mutate(t=1:nrow(output)) %>% pivot_longer(!t)
# plot
ggplot(output, aes(x=t,y=value)) + geom_line() + facet_wrap(~name) + theme_bw() + standard_theme

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## run_sample is in plot_helpers.R file
est_ep <- run_sample(post, TRUE)
# Predicted number of positive samples (useful for comparing the model-predicted incidence with the data)
est <- run_sample(post, FALSE)  # Predicted incidence

# Output the figs for each age group
for (a in 1:ageGroupNumber){
  dataplot <- data.frame(model_m=est_ep[[a]][2,], model_l=est_ep[[a]][1,], model_u=est_ep[[a]][3,], 
                         data=observationalData[,a+1], time=c(1:(52*7)))
  ggplot() + geom_point(data = dataplot, aes(x = time, y = data, color = 'black'), size=0.5, alpha=0.5) +
    geom_line(data=dataplot, aes(x = time, y = model_m, color = 'red'), size=0.4, linetype="dotted") +
    geom_line(data=dataplot, aes(x = time, y = model_l), color = 'red', size=0.4, alpha=0.5) +
    geom_line(data=dataplot, aes(x = time, y = model_u), color = 'red', size=0.4, alpha=0.5) +
    geom_ribbon(data=dataplot, aes(x = time, ymin=model_l, ymax=model_u), fill="red", alpha=0.5) +
    scale_colour_manual(name="",values=c('black'='black','red'='red'),labels=c('Observational data','Model-predicted')) +
    xlab('Week') + ylab('Estimated number of weekly new positive samples') + ggtitle(paste("Age group", a))
  # save
  ggsave(paste0(file.path, "/figs/compar/inci_", a, ".pdf"))
}

# all agegroup outputs
data.frame(t(bind_rows(data.frame(est)))) %>% 
  mutate(time=rep(1:length(est_ep[[1]][1,]),length(est_ep)),
         agegroup=unlist(lapply(1:length(est_ep), function(x) rep(x,length(est_ep[[1]][1,]))))) %>%
ggplot(aes(x=time)) + geom_line(aes(y=X2)) + geom_ribbon(aes(ymin=X1,ymax=X3),fill="red",alpha=0.5) +
  facet_wrap(~agegroup,scales="free_y") + theme_bw() + standard_theme

# Output the figs for each age group
for (a in 1:ageGroupNumber){
  dataplot <- data.frame(model_m=est[[a]][2,], model_l=est[[a]][1,], model_u=est[[a]][3,], time=c(1:(52*7)))
  ggplot() + geom_line(data = dataplot, aes(x = time, y = model_m, color = 'red'), size=0.4, linetype="dotted") +
    geom_line(data = dataplot, aes(x = time, y = model_l), color = 'red', size=0.4, alpha=0.5) +
    geom_line(data = dataplot, aes(x = time, y = model_u), color = 'red', size=0.4, alpha=0.5) +
    geom_ribbon(data = dataplot, aes(x = time, ymin=model_l, ymax=model_u), fill="red", alpha=0.5) + 
    xlab('Week') + ylab('Estimated number of weekly new infections') + ggtitle(paste("Age group", a))
  # save
  ggsave(paste0(file.path, "/figs/pred/inci_", a, ".pdf"))
}

### 3.4 Plot the annual incidences using the posterior samples.
ann_est <- run_sample_annual(post)
dataplot <- data.frame(model_m=ann_est[2,], model_l=ann_est[1,], model_u=ann_est[3,], age=1:ageGroupNumber)

# Output the annual incidence fig
ggplot() + geom_line(data = dataplot, aes(x = age, y = model_m, color = 'red'), size=0.4, linetype="dotted") +
  geom_line(data = dataplot, aes(x = age, y = model_l), color = 'red', size=0.4, alpha=0.5) +
  geom_line(data = dataplot, aes(x = age, y = model_u), color = 'red', size=0.4, alpha=0.5) +
  geom_ribbon(data = dataplot, aes(x = age, ymin=model_l, ymax=model_u), fill="red", alpha=0.5) + 
  xlab('Age group') + ylab('Proportion of age group who acquire infection  annually') + ggtitle("Annual incidence")
# save
ggsave(paste0(file.path, "/figs/inci_all.pdf"))

### 3.5 Plot the proportion of infants born with protection
pRoutput <- run_sample_maternal_protected(post)
pR <- data.frame(model_m=pRoutput[2,], model_l=pRoutput[1,], model_u=pRoutput[3,], time=1:364)

ggplot() + geom_line(data=pR, aes(x = time, y = model_m), color = 'red', size=0.4, alpha=0.5) +
  geom_ribbon(data = pR, aes(x = time, ymin=model_l, ymax=model_u), fill="red", alpha=0.5) + theme_bw() + standard_theme

write.csv(pR, file=paste0(file.path, "/figs/pRdata.csv"))