# example to measure accuracy of the integration methods

# calculate accuracy
ck.accuracy <- function(comp, numinteg, tol=1e-8, digits=16, methods=NA)
{
  stopifnot(!is.na(methods))
  acdigs <- paste("%.",digits,"f",sep="")
  cat("\nChecking Numerical Integration Methods\n",rep("-",44),sep="")
  cat("\ncomparative method\t= ",methods[1],sep="")
  cat("\nnum. int. method\t= ",methods[2],sep="")
  cat("\ncomparative value\t= ",comp,sep="")
  cat("\nnum. integ. value\t= ",sprintf(acdigs, numinteg),sep="")
  diff <- comp-numinteg
  cat("\ndiff comp-num.integ.\t= ",sprintf(acdigs, diff),sep="")
  cat("\ncomp. > num.integ.\t= ",comp > numinteg,sep="")
  cat("\nabs(diff) < tol = ",tol,"\t= ",abs(diff) < tol,sep="")
  cat("\nfactor comp/num.integ.\t= ",sprintf(acdigs, comp/numinteg),"\n\n",sep="")
  cat("---\nNote: 'comparative' value can be a real value or from a different method")
}

# get resources and helper functions
library(Rmpfr)
source("numinteg_brob.r")
source("simpsonrule_nlb.r")

# define the same functions in different ways
f <- function(x) sin(x)
f.log <- function(x) log(sin(x))
f.brob <- function(x) as.brob(sin(x))
f.Rmpfr <- function(x) sin(mpfr(x,precBits=256))

# integration boundaries
lower <- 0
upper <- pi

# plot
seku <- seq(-pi,pi,length=101)
par(mfrow=c(1,2))
# 2nd half-period from 0 to pi
plot(seku[51:100], sin(seku[51:100]), col="darkred", type="l", pre.plot=grid(), bty="n", xlab="x", ylab="f(x)")
# full-period from -pi to pi
plot(seku, sin(seku), col="darkred", type="l", pre.plot=grid(), bty="n", xlab="x", ylab="f(x)")

# integral analytically solved
# sin(x) dx -> -cos(x) + constant
f.integ <- function(lower,upper) -cos(upper) -(-cos(lower))

# example -> define limits to cover half-period of the integral
int.real.v <- f.integ(lower=lower, upper=upper)

# integrate
int.integrate.v <- integrate(f, lower, upper)$v

# integrateR (Rmpfr)
# numerical integration
int.integrateR.v <- Rmpfr:::integrateR(f.Rmpfr, lower, upper)$value

# trapez (pracma)
trapez.v <- pracma:::trapz(x=sek, y=f(sek))

Nsteps <- 1e4

# log
simpson.l.v <- exp(simpsonrule.nlb(f.log, lower, upper, type="log", Nsteps=Nsteps))

# brob
sek <- seq(lower,upper, length=Nsteps+1)
f.sek.brob <- lapply(seq_along(sek), function(x) f.brob(sek[x]))
f.sek.brob.c <- list2vec.brob(f.sek.brob)
trapez.brob.v <- as.numeric( trapz.brob(x=sek, y=f.sek.brob.c) )

int.real.v
int.integate.v
int.integrateR.v
trapez.v
simpson.l.v
trapez.brob.v

# check accuracy
ck.accuracy(comp=int.real.v, numinteg=int.integrate.v, methods=c("real value","integrate()"))
ck.accuracy(comp=int.real.v, numinteg=as.numeric(int.integrateR.v), methods=c("real value","integrateR()"))
ck.accuracy(comp=int.real.v, numinteg=trapez.v,methods=c("real value","pracma:::trapz()"))
ck.accuracy(comp=int.real.v, numinteg=simpson.l.v, methods=c("real value","simpson log"))
ck.accuracy(comp=int.real.v, numinteg=trapez.brob.v, methods=c("real value","trapez brob"))

