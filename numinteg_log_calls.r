library(cmna)
library(microbenchmark)

source("numinteg_helper-funcs.r")
source("numinteg_log.r")

# trapez

# define functions
f <- function(x) { sin(x)^2 + cos(x)^2 }
f.log <- function(x) { log( sin(x)^2 + cos(x)^2 ) }
f.log2 <- function(x) .llog2add.short( c(2*log(sin(x)), 2*log(cos(x))) )

lower <- -pi
upper <- pi
Nsteps <- 100

# calculate
cmna:::trap(f, -pi, pi, m=Nsteps)
trap.nl(f, lower, upper, Nsteps, method="normal")

log( cmna:::trap(f, -pi, pi, m=Nsteps) )
trap.nl(f.log, lower, upper, Nsteps, method="log")
trap.nl(f.log2, lower, upper, Nsteps, method="log")

# benchmarking
mbench.res1 <- microbenchmark( 
  cmna:::trap(f, -pi, pi, m=Nsteps),
  log( cmna:::trap(f, -pi, pi, m=Nsteps) ),
  trap.nl(f, lower, upper, Nsteps, method="normal"),
  trap.nl(f.log, lower, upper, Nsteps, method="log"),
  trap.nl(f.log2, lower, upper, Nsteps, method="log")
)

mbench.res1

# romberg

# define functions
f <- function(x) { sin(x)^2 + log(x) }
f.log <- function(x) { log( sin(x)^2 + log(x) ) }
f.log2 <- function(x) { .llog.2add.short( c( 2*log(sin(x)), log(log(x)) ) )}

lower <- 1  # a
upper <- 10 # b
Nsteps <- 3 # m

# calculate
log(cmna:::romberg(f, lower, upper, m=Nsteps, tab=FALSE))
log(romberg.nl(f, lower, upper, Nsteps, tab=FALSE, method="normal"))
romberg.nl(f.log, lower, upper, Nsteps, tab=FALSE, method="log")
romberg.nl(f.log2, lower, upper, Nsteps, tab=FALSE, method="log")

log(cmna:::romberg(f, lower, upper, m=Nsteps, tab=TRUE))
log(romberg.nl(f, lower, upper, Nsteps, tab=TRUE, method="normal"))
romberg.nl(f.log, lower, upper, Nsteps, tab=TRUE, method="log")
romberg.nl(f.log2, lower, upper, Nsteps, tab=TRUE, method="log")

# benchmarking
mbench.res2 <- microbenchmark( 
         log( cmna:::romberg(f, lower, upper, m=Nsteps, tab=FALSE) ),
				 log( romberg.nl(f, lower, upper, Nsteps, tab=FALSE, method="normal") ),
				 romberg.nl( f.log, lower, upper, Nsteps, tab=FALSE, method="log" ),
				 romberg.nl( f.log2, lower, upper, Nsteps, tab=FALSE, method="log" )
			   )

mbench.res2

