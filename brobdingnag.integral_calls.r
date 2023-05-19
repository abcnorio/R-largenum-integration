source("brob.integral.r")

fx <- function(x) sin(x)
f <- function(x) as.brob(sin(x))


#########################################################################
# helper function to convert a brob list to a brob vector
#
# call:
list2vec.brob( lapply(c(1:10),as.brob))
#########################################################################


#########################################################################
# helper function to replace '%*%' scalarproduct
# that does not work for brob objects
#
# call:
as.brob(c(1:10) %*% c(11:20))
scalarprod.brob(as.brob(1:10), as.brob(11:20))
scalarprod.brob(as.brob(1:10), c(11:20))
#########################################################################


#########################################################################
# Bolstad2:::sintegral
#
# call:
sek <- seq(0,pi, length=1000)
funx <- fx(sek)
as.brob(sintegral(sek, funx)$int)
sintegral.brob.parallel(fx=f, sL=0, sH=pi, Nsteps=500)
#########################################################################


#########################################################################
# pracma:::romberg
# call:
f
fx
a <- 0
b <- pi
maxit <- 25
tol <- 1e-15
romberg.brob(f, a, b, tol=1e-15)
as.brob(romberg(f=fx, a, b, tol=1e-15)$value)
#########################################################################


#########################################################################
# pracma:::cotes
#
# tweaked 'cotes' from 'pracma' to work with brob
# TODO parallel computing version
# call:
f
fx
#
for(i in 2:8)
{
  cat("\ni = ",i,"\n")
  print( as.brob( cotes(fx, 0.1, 1.2, 20, i) ))
  print( cotes.brob(f, 0.1, 1.2, 20, i) )
}
# example from 'cotes' manpage
for(i in 2:8)
{
  cat("\ni = ",i,"\n")  
  print( as.brob(cotes(fx, 0, pi/2, 20, i) ))
  print( cotes.brob(f, 0, pi/2, 20, i) )
}
#########################################################################


######################################################################### 
# pracma:::integral
#
# call:
f
fx
no_intervals <- 8
random <- FALSE
reltol <- 1e-08
# pracma 'integral' 'Kronrod'
integral.brob(fun=f, xmin=0, xmax=pi, method="Kronrod", no_intervals=no_intervals)
as.brob(pracma:::integral(fun=fx, xmin=0, xmax=pi, method="Kronrod"))
# pracma 'integral' 'Simpson'
integral.brob(fun=f, xmin=0, xmax=pi, method="Simpson", no_intervals=no_intervals)
as.brob(pracma:::integral(fun=fx, xmin=0, xmax=pi, method="Simpson"))
# pracma 'integral' 'Simpson'
integral.brob(fun=f, xmin=0, xmax=pi, method="Simpson", no_intervals=no_intervals, random=TRUE)
as.brob(pracma:::integral(fun=fx, xmin=0, xmax=pi, method="Simpson", random=TRUE))
#########################################################################


#########################################################################
# pracma:::quadgk
#
# call:
f
fx
a
b
tol <- 1e-15
quadgk.brob(f, a, b)
as.brob(quadgk(fx, a, b))
#########################################################################


#########################################################################
# pracma:::simpadpt
#
# call:
f
fx
a
b
simpadpt.brob(f,a,b, tol=1e-08)
as.brob(simpadpt(fx,a,b))
#########################################################################


#########################################################################
# pracma:::quadv
#
# seems to be just a vectorized version, no need to add this
# if one wants a vectorized version that works on multiple functions
# in parallel, one can write a short wrapper around an integral
# function
#
# otherwise it would require to simulate matrices for brob objects
#
#########################################################################


#########################################################################
# pracma:::quadinf
#
# seems to work without modification for brob
# call:
# from manpage
f1 <- function(x) as.brob(exp(-x^2))  # sqrt(pi)/2         theory
f2 <- function(x) exp(-x^2)
as.brob(quadinf(f2, 0, Inf)$Q)        # 0.8862269254527413
quadinf.brob(f1, 0, Inf)$Q
as.brob(quadinf(f2, -Inf, 0)$Q)       # 0.8862269254527413
quadinf.brob(f1, -Inf, 0)$Q
#########################################################################


#########################################################################
# pracma:::quad
#
# call:
f
fx
xa <- a
xb <- b
xa
xb
f3 <- function(x) x * cos(0.1*exp(x)) * sin(0.1*pi*exp(x))
f4 <- function(x) as.brob(x * cos(0.1*exp(x)) * sin(0.1*pi*exp(x)))
as.brob(pracma:::quad(f3, 0, 4))
quad.brob(f4,0,4)
as.brob(pracma:::quad(fx,xa=a,xb=b))
quad.brob(f,xa=a,xb=b, TRACE=FALSE)
quad.brob(f,xa=a,xb=b, TRACE=TRUE, digs=7)
#########################################################################


#########################################################################
# pracma:::trapz
#
# call:
f
fx
a
b
sek <- seq(a,b, length.out=101)
f.sek.brob <- lapply(seq_along(sek), function(x) f(sek[x]))
head(f.sek.brob)
f.sek.brob.c <- list2vec.brob(f.sek.brob)
f.sek.brob.c
trapz.brob(x=sek, y=f.sek.brob.c)
as.brob(trapz(x=sek, y=fx(sek)))
#########################################################################


#########################################################################
# pracma:::quadgr
# pracma:::.rich
#
# associated to quadgr - just works, only renamed here...
#
# call:
f
a
b
quadgr.brob(f,a,b)
as.brob(pracma:::quadgr(fx,a,b)$value)
# example from quadgr
flog <- function(t) log(1-t)/t
flog.brob <- function(t) as.brob(log(1-t)/t)
as.brob(quadgr(flog, 1, 0, tol = 1e-15)$value)
quadgr.brob(flog.brob, 1, 0, tol = 1e-15)
#########################################################################

