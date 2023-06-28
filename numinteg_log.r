# numerical integation
# log routines
# trapez, romberg
# base: cmna


################################################################################
# trapez numerical integration
trap.nl <- function(f, lower, upper, Nsteps=100, type="normal")
{
  x <- seq(lower, upper, length.out=Nsteps + 1)
  if(type == "normal")
  {
    fx <- f(x)
    p.area <- sum((fx[2:(Nsteps + 1)] + fx[1:Nsteps]))
    res <- p.area * abs(upper - lower)/(2 * Nsteps)
  } else if(type == "log")
  {
    fx.log <- f.log(x)
	p.area.log <- .llog.2add.short( c( fx.log[2:(Nsteps + 1)], fx.log[1:Nsteps] ) )
    res <- p.area.log + log( abs(upper - lower) / (2 * Nsteps) )
  }
return(res)
}
# calls:
#
# f <- function(x) { sin(x)^2 + log(x) }
# f.log <- function(x) { log( sin(x)^2 + log(x) ) }
# a <- lower <- 1
# b <- upper <- 10
# m <- Nsteps <- 3
#
# trap.nl(f, lower, upper, Nsteps, type="normal")
# exp(trap.nl(f.log, lower, upper, Nsteps, type="log"))
################################################################################


################################################################################
# romberg numerical integration
romberg.nl <- function(f, lower, upper, Nsteps, tab=FALSE, type="normal")
{ 
  if(type == "normal")
  {
    R <- matrix(NA, nrow=Nsteps, ncol=Nsteps)
    R[1, 1] <- trap.nl(f, lower, upper, Nsteps=1, type="normal")
    for(j in 2:Nsteps)
    {
      R[j, 1] <- trap.nl(f, lower, upper, Nsteps=2^(j - 1), type="normal")
      for(k in 2:j)
      {
        k4 <- 4^(k - 1)
        R[j, k] <- k4 * R[j, k - 1] - R[j - 1, k - 1]
        R[j, k] <- R[j, k] / (k4 - 1)
        
      }
    }
    ifelse(tab == TRUE, res <- R, res <- R[Nsteps, Nsteps])
  } else if(type == "log")
  {
    R.log <- matrix(NA, nrow=Nsteps, ncol=Nsteps)
    R.log[1, 1] <- trap.nl(f.log, lower, upper, Nsteps=1, type="log")
    for(j in 2:Nsteps)
    {
      R.log[j, 1] <- trap.nl(f.log, lower, upper, Nsteps=2^(j - 1), type="log")
      for(k in 2:j)
      {
        k4.log <- (k - 1) * log(4)
        R.log[j, k] <- .llog2sub.short( (k4.log + R.log[j, k - 1]), R.log[j - 1, k - 1] )
        R.log[j, k] <- R.log[j, k] - .llog2sub.short(k4.log, log(1))
      }
    }
    ifelse(tab == TRUE, res <- R.log, res <- R.log[Nsteps, Nsteps])
  }
return(res)
}
# calls:
#
# f <- function(x) { sin(x)^2 + log(x) }
# f.log <- function(x) { log( sin(x)^2 + log(x) ) }
# a <- lower <- 1
# b <- upper <- 10
# m <- Nsteps <- 3
#
# log(romberg.nl(f, lower, upper, Nsteps, tab=TRUE, type="normal"))
# romberg.nl(f.log, lower, upper, Nsteps, tab=TRUE, type="log")
################################################################################

