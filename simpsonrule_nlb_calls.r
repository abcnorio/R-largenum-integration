# example calls for Simpson rule for
#
# normal
# log
# brob (Brobdingnag)
#
# and compare to
# integrate
# sintegral (Bolstad2)
# integrateR (Rmpfr)

library(Brobdingnag)
library(Bolstad2)
library(Rmpfr)

source("numinteg_helper-funcs.r")
source("simpsonrule_nlb.r")

# define same fucntion for different methods
f <- function(x) sin(x)
f.log <- function(x) log(sin(x))
f.brob <- function(x) brob(log(sin(x)))
f2.brob <- function(x) as.brob(sin(x))
f.Rmpfr <- function(x) sin(mpfr(x,precBits=256))

# add lower and upper limits and Nsteps
lower <- 0
upper <- pi
Nsteps <- 1e3

x <- seq(lower,upper,length=Nsteps)
fx.sintegral <- f(x)

# calculate and review accuracy

# normal
sprintf("%.16e", simpsonrule.nlb(fx=f, lower=lower, upper=upper, type="normal",Nsteps=Nsteps) )
# log
sprintf("%.16e", exp(simpsonrule.nlb(fx=f.log, lower=lower, upper=upper, type="log",Nsteps=Nsteps)) )
# log parallel
sprintf("%.16e", exp(simpsonrule.nlb(fx=f.log, lower=lower, upper=upper, type="log",Nsteps=Nsteps, parallel=TRUE)) )
# via Brobdingnag alt 1
sprintf("%.16e", as.numeric(simpsonrule.nlb(fx=f.brob, lower=lower, upper=upper, type="brob",Nsteps=Nsteps)) )
# via Brobdingnag alt 2
sprintf("%.16e", as.numeric(simpsonrule.nlb(fx=f2.brob, lower=lower, upper=upper, type="brob",Nsteps=Nsteps)) )
# via integrate
sprintf("%.16e", integrate(f,lower,upper)$v )
# via sintegral (Bolstad2)
sprintf("%.16e", sintegral(x,fx.sintegral)$int )
# via Rmpfr
integrateR(f.Rmpfr,lower,upper)$value

# benchmark
library(microbenchmark)

micro.res <- microbenchmark( simpsonrule.nlb(fx=f, lower=lower, upper=upper, type="normal",Nsteps=Nsteps) ,
			   exp(simpsonrule.nlb(fx=f.log, lower=lower, upper=upper, type="log",Nsteps=Nsteps)) ,
			   exp(simpsonrule.nlb(fx=f.log, lower=lower, upper=upper, type="log",Nsteps=Nsteps, parallel=TRUE)) ,
			   as.numeric(simpsonrule.nlb(fx=f.brob, lower=lower, upper=upper, type="brob",Nsteps=Nsteps)) ,
			   as.numeric(simpsonrule.nlb(fx=f2.brob, lower=lower, upper=upper, type="brob",Nsteps=Nsteps)) ,
			   integrate(f,lower,upper)$v,
			   sintegral(x,fx.sintegral)$int,
			   integrateR(f.Rmpfr,lower,upper)$value
			 )
print(micro.res)
# see differences between the two brob versions although the same (!) numerical integration routine is used
