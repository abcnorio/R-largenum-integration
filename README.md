# R code to enable integration functions to work with very large numbers

## Overview

The [Brobdingnag](https://github.com/RobinHankin/Brobdingnag) R package provides functions to enable the work with very large numbers. Most integration functions fail if large numbers temporarily occur. The R code here implements functions of the R package Brobdingnag to enable integration functions from R packages [pracma](https://github.com/cran/pracma) and [Bolstad2](https://github.com/cran/Bolstad2) to work with very large numbers. This is useful if someone tries to integrate e.g. Gamma or other functions that easily reach values beyond normal tolerance ie. near zero or up to infinite. R's tolerance is

```
> 10^(307:309)
[1] 1e+307 1e+308    Inf
> 10^(-323:-324)
[1] 9.881313e-324  0.000000e+00
```

## License

see license file included - in short it is GPL >=2 - in accordance to the packages used.

## Files

- **brobdingnag.integral.r** = tweaked functions to work with very large numbers
- **brobdingnag.integral_calls** = example code that shows that the tweaked functions produce the same results as the original functions

## Functions used - original names just got a *.brob added:

- list2vec.brob = convert a Brobdingnag list to a Brobdingnag vector
- scalarprod.brob = replace '%*%' scalarproduct that does not work for Brobdingnag objects

from [Bolstad2](https://github.com/cran/Bolstad2):

- sintegral.brob.parallel

from [pracma](https://github.com/cran/pracma):

- romberg.brob
- cotes.brob
- integral.brob ("Kronrod", "Simpson" but not "Clenshaw")
- quadgk.brob
- simpadpt.brob
- quadinf.brob
- quad.brob
- trapz.brob
- quadgr.brob

Not all integral functions of pracma are covered. What is missing are:

- Clenshaw from integral
- integral2
- integral3

## Function calls

Functions are called identical to the original version.

```
source("brobdingnag.integral.r")
[...use as always but add .brob to the name of the function call - see above]
```

e.g.

```
fx <- function(x) sin(x)
f <- function(x) as.brob(sin(x))
a <- 0
b <- pi
maxit <- 25
tol <- 1e-15
romberg.brob(f, a, b, tol=1e-15)
as.brob(romberg(f=fx, a, b, tol=1e-15)$value)
```

## R version

All R scripts should work under R >=v3.

## Disclaimer

R code was tested but is provided "as is".


