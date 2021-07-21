# Function to run a sample
# classEvaluateLogLikelihood$getWeeklySampleCpp(pars, bool),
# run this function with a vector of parameters (pars) to get the incidence output

# run_sample: runs the model with 1000 random samples from the posterior, outputs the median, and 95% confidence intervals for the weekly incidence per age group
# post: the posterior samples as imported (post)
# bool: flag to output the model multipled by the ascertainment probability (epsilon). 
# If FALSE then outputs the true predicted values.
run_sample <- function(post, bool){
  a <- rep(NA, 25*(52*7)*1000); arr <- array(a, c(52*7,25,1000))
  for (i in 1:1000) {
    j <- rdunif(1, 1, nrow(post))
    arr[,, i] <- classEvaluateLogLikelihood$getWeeklySampleCpp(as.numeric(post[j,]), bool)   }
  out <- lapply(1:25, function(y) sapply(1:(52*7), function(x) sort(arr[x,y,])[c(25,500,975)]))  
  out }

# run_sample_annual: runs the model with 1000 random samples from the posterior, outputs the median, and 95% confidence intervals for the annual incidence per age group
run_sample_annual <- function(post){
  a <- rep(NA, 25*1000); arr<-array(a, c(25,  1000))
  for (i in 1:1000) { j <- rdunif(1, 1, nrow(post))
    arr[, i] <- classEvaluateLogLikelihood$getAnnualIncidenceCpp(as.numeric(post[j,])) }
  out <- sapply(1:25, function(s ) sort(arr[s,])[c(25, 500 , 975)])
  out
}

run_sample_maternal_protected <- function(post){
  a <- rep(NA, 364*1000); arr <- array(a,c(364,1000))
  for (i in 1:1000){
    j <- rdunif(1, 1, nrow(post))
    arr[,i] <- classEvaluateLogLikelihood$getProportionBornProtectedCpp(as.numeric(post[j,]))   }
  out <- sapply(1:364, function(s) sort(arr[s,])[c(25,500,975)])
  out }
