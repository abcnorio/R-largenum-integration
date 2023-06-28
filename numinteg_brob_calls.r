require(Brobdingnag)
require(pracma)
require(Bolstad2)
require(parallel)

source("numinteg_helper-funcs.r")
source("numinteg_brob.r")


# define functions
f <- function(x) sin(x)
f.brob <- function(x) as.brob(sin(x))
f2.brob <- function(x) brob(log(sin(x)))


# helper function to convert a brob list to a brob vector
#
# call:
list2vec.brob( lapply(c(1:10),as.brob))


# helper function to replace '%*%' scalarproduct
# that does not work for brob objects
#
# call:
log( c(1:10) %*% c(11:20) )
as.brob(c(1:10) %*% c(11:20))
# does not work:
as.brob(c(1:10) %*% as.brob(c(11:20)))
# therefor:
scalarprod.brob(as.brob(1:10), c(11:20))
scalarprod.brob(as.brob(1:10), as.brob(11:20))


# Bolstad2:::sintegral
#
# call:
lower <- 0
upper <- pi
Nsteps <- 1e3
sek <- seq(lower, upper, length=Nsteps)
funx <- f(sek)
as.brob(sintegral(sek, funx)$int)
sintegral.brob.parallel(fx=f.brob, sL=lower, sH=upper, Nsteps=Nsteps)
sintegral.brob.parallel(fx=f2.brob, sL=lower, sH=upper, Nsteps=Nsteps)


# pracma:::romberg
# call:
f
f.brob
a <- 0
b <- pi
maxit <- 25
tol <- 1e-15
as.brob(romberg(f=f, a, b, maxit=maxit, tol=1e-15)$value)
romberg.brob(f.brob, a, b, maxit=maxit, tol=1e-15)$value


# pracma:::cotes
#
# tweaked 'cotes' from 'pracma' to work with brob
# TODO parallel computing version
# call:
f
f.brob
#
# example from 'cotes' manpage
for(i in 2:8)
{
  cat("\ni = ",i,"\n")
  print( as.brob( cotes(f, 0.1, 1.2, 20, i) ))
  print( cotes.brob(f.brob, 0.1, 1.2, 20, i) )
}


# pracma:::integral
#
# call:
f
f.brob
no_intervals <- 8
random <- FALSE
reltol <- 1e-08
#
# pracma 'integral' 'Kronrod'
as.brob(pracma:::integral(fun=f, xmin=0, xmax=pi, method="Kronrod"))
integral.brob(fun=f.brob, xmin=0, xmax=pi, method="Kronrod", no_intervals=no_intervals)

# pracma 'integral' 'Simpson'
as.brob(pracma:::integral(fun=f, xmin=0, xmax=pi, method="Simpson"))
integral.brob(fun=f.brob, xmin=0, xmax=pi, method="Simpson", no_intervals=no_intervals)

# pracma 'integral' 'Simpson'
as.brob(pracma:::integral(fun=f, xmin=0, xmax=pi, method="Simpson", random=TRUE))
integral.brob(fun=f.brob, xmin=0, xmax=pi, method="Simpson", no_intervals=no_intervals, random=TRUE)


# pracma:::quadgk
#
# call:
f
f.brob
a
b
tol <- 1e-15
as.brob(quadgk(f, a, b))
quadgk.brob(f.brob, a, b)


# pracma:::simpadpt
#
# call:
f
f.brob
a
b
as.brob(simpadpt(f,a,b))
simpadpt.brob(f.brob ,a,b, tol=1e-08)


# pracma:::quadv
#
# seems to be just a vectorized version, no need to add this
# if one wants a vectorized version that works on multiple functions
# in parallel, one can write a short wrapper around an integral
# function
#
# otherwise it would require to simulate matrices for brob objects
#


# pracma:::quadinf
#
# seems to work without modification for brob
# call:
# from manpage
f1 <- function(x) exp(-x^2)
f1.brob <- function(x) as.brob(exp(-x^2))  # sqrt(pi)/2
#
as.brob(quadinf(f1, 0, Inf)$Q)        # 0.8862269254527413
quadinf.brob(f1.brob, 0, Inf)$Q


# pracma:::quad
#
# call:
f
f.brob
xa <- a
xb <- b
xa
xb
f2 <- function(x) x * cos(0.1*exp(x)) * sin(0.1*pi*exp(x))
f2.brob <- function(x) as.brob(x * cos(0.1*exp(x)) * sin(0.1*pi*exp(x)))
#
as.brob(pracma:::quad(f2, 0, 4))
quad.brob(f2.brob,0,4)
#
as.brob(pracma:::quad(f2, xa=a, xb=b))
quad.brob(f2.brob, xa=a, xb=b, TRACE=FALSE)
quad.brob(f2.brob, xa=a, xb=b, TRACE=TRUE, digs=7)


# pracma:::trapz
#
# call:
f
f.brob
a
b
Nsteps <- 1e2
sek <- seq(a,b, length.out=Nsteps+1)
f.sek.brob <- lapply(seq_along(sek), function(x) f.brob(sek[x]))
head(f.sek.brob)
f.sek.brob.c <- list2vec.brob(f.sek.brob)
f.sek.brob.c
#
as.brob(trapz(x=sek, y=f(sek)))
trapz.brob(x=sek, y=f.sek.brob.c)


# pracma:::quadgr
# pracma:::.rich
#
# associated to quadgr - just works, only renamed here...
#
# call:
f
f.brob
a
b
#
as.brob(pracma:::quadgr(f,a,b)$value)
quadgr.brob(f.brob,a,b)$value

# example from quadgr
f.log <- function(t) log(1-t)/t
f.log.brob <- function(t) as.brob(log(1-t)/t)
as.brob(quadgr(f.log, 1, 0, tol = 1e-15)$value)
quadgr.brob(f.log.brob, 1, 0, tol = 1e-15)$value

