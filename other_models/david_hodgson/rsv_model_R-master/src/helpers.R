logit <- function(x, a, b)
{
  log((x-a)/(b-x))
} 

inv.logit <- function(x, a, b)
{
  a + (b-a)*exp(x)/(exp(x)+1)
} 

D.inv.logit <- function(x, a, b)
{
  (b-a)*exp(x)/(exp(x)+1)^2
} 

logit.prior <- function (x, a, b, FUNC, bool, ...)
{
    if (bool == TRUE)
        out = log(FUNC(inv.logit(x, a, b), ...)*D.inv.logit(x, a, b))
    else
        out = FUNC(inv.logit(x, a, b), ...)*D.inv.logit(x, a, b)
    
    out
}

logit.prior.sample <- function(a, b, FUNC, ...)
{
  x <- FUNC(1, ...)
  while ((x < a) | (x > b)){
    x <- FUNC(1, ...)
  }
  
  logit(x, a, b)
}

prior <- function (x, a, b, FUNC, bool, ...)
{
    if (bool == TRUE)
        out = log(FUNC(x, ...))
    else
        out = FUNC(x, ...)
    
    out
}

prior.sample <- function(a, b, FUNC, ...)
{
    x <- FUNC(1, ...)
    while ((x < a) | (x > b)){
        x <- FUNC(1, ...)
    }
    
    x
}
