################################################################
# Vector-borne transmission of plague
# Model function (wrapper)

# Author: Fabienne Krauer, University of Oslo
# Created: 05.10.2017
# Last updated: 19.04.2020
################################################################

meta <- function(theta, inits, time, solver, diseasemodel, M) {
  

  traj <- data.frame(ode(y=inits,
                         times=time, 
                         func=match.fun(diseasemodel), 
                         parms=theta,
                         method=solver,
                         M=M))
  
  # calculate total cumulative deaths
  traj$Dhtot <- rowSums(exp(traj[, grep("Dh", colnames(traj)), drop=FALSE]))
  # calculate total incident deaths
  traj$Dhinc <- c(theta[["Dh0"]], diff(exp(traj$Dhtot), lag=1))

  return(traj)
  
}

