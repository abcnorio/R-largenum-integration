# R code to enable integration functions to work with very large numbers

## Overview

The [Brobdingnag](https://github.com/RobinHankin/Brobdingnag) R package provides functions to enable the work with very large numbers. Most integration functions fail if large numbers temporarily occur. The R code here implements functions of the R package Brobdingnag to enable integration functions from R packages [pracma](https://github.com/cran/pracma) and [Bolstad2](https://github.com/cran/Bolstad2) to work with very large numbers. This is helpful if one uses functions like Gamma or others within integrals that easily reach values beyond normal tolerance ie. near zero or up to infinity. R's tolerance is

```
> 10^(307:309)
[1] 1e+307 1e+308    Inf
> 10^(-323:-324)
[1] 9.881313e-324  0.000000e+00
```

## Licenses

See license file included - in short it is GPL >=2 - in accordance to the packages used and the license associated with them.

## Files

| File | Description |
| --- | --- |
| `brobdingnag.integral.r` | tweaked functions to work with very large numbers |
| `brobdingnag.integral_calls.r` | example code that shows that the tweaked functions produce the same results as the original functions |

## Functions used

Original names just got a `*.brob` added at the end of the filename.

| Function | Description |
| --- | --- |
| `list2vec.brob` | convert a Brobdingnag list to a Brobdingnag vector |
| `scalarprod.brob` | replace `%*%` scalarproduct that does not work for Brobdingnag objects |
| `sintegral.brob.parallel` | from [Bolstad2](https://github.com/cran/Bolstad2) |
| | from [pracma](https://github.com/cran/pracma): |
| `romberg.brob` | Romberg |
| `cotes.brob` | Cotes |
| `integral.brob` | Kronrod, Simpson but not Clenshaw |
| `quadgk.brob` | internal function |
| `simpadpt.brob` | Simpson |
| `quadinf.brob` | internal function |
| `quad.brob` | internal function |
| `trapz.brob` | Trapez |
| `quadgr.brob` | internal function |

Not all integration functions of pracma are covered. What is missing are:

- Clenshaw from `integral`
- `integral2`
- `integral3`

## Function calls

Functions are called identical to the original version. Use as normal but add .brob to the name of the function call and use a function that gives out a Brobdingnag object that acts as input for the integration function.

```
source("brobdingnag.integral.r")
```

e.g.

```
fx <- function(x) sin(x)
f <- function(x) as.brob(sin(x))
a <- 0
b <- pi
maxit <- 25
tol <- 1e-15
romberg.brob(f, a, b, maxit=maxit, tol=tol)
as.brob(romberg(f=fx, a, b, maxit=maxit, tol=tol)$value)
```

## TODO

Cover all integration functions.

## R version

All R scripts should work under R >=v3.

## Disclaimer


The R code was tested carefully and cross-checked against various public available results and manpages to ensure proper results.  The example calls can be used to test for accuracy. However, it is provided "as is". Use common sense to compare results with expectations. NO WARRANTY of any kind is involved here. There is no guarantee that the software is free of error or consistent with any standards or even meets your requirements. Do not use the software or rely on it to solve problems if incorrect results may lead to hurting or injurying living beings of any kind or if it can lead to loss of property or any other possible damage. If you use the software in such a manner, you are on your own and it is your own risk.


