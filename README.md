# R code to enable numerical integration functions to work with very large numbers

## Overview

The [Brobdingnag](https://github.com/RobinHankin/Brobdingnag) R package provides functions to enable the work with very large numbers. Most integration functions fail if large numbers temporarily occur. The R code here implements functions of the R package Brobdingnag to enable integration functions from R packages [pracma](https://github.com/cran/pracma) and [Bolstad2](https://github.com/cran/Bolstad2) to work with very large numbers. This is helpful if one uses functions like Gamma or others within integrals that easily reach values beyond normal tolerance ie. near zero or up to infinity. R's tolerance on a machine is

```
> 10^(307:309)
[1] 1e+307 1e+308    Inf
> 10^(-323:-324)
[1] 9.881313e-324  0.000000e+00
```

or more precisely

```
> .Machine$double.xmin
[1] 2.225074e-308
> .Machine$double.xmax
[1] 1.797693e+308
```

Another option is to transform all functions for which numerical integration is required to the log scale and then use a numerical integration method that is written for functions with log output. Normally this speeds up the process drastically compared to the usage of Brobdingnag while attaining the same accuracy, because the underlying R code is more or less exactly the same. However, it is slower than using non-log methods like `integrate`. Concretely, this means

```
fx <- function(x) sin(x)/cos(x)
```

becomes

```
fx.log <- function(x) log(sin(x)) - log(cos(x))
```

This conversion is no automatic process and for complicated integrals it requires a certain amount of effort and testing, esp. to make it efficient.

The third option is to use the package [Rmpfr](https://rmpfr.r-forge.r-project.org/) for very large numbers. It has a command called `integrateR` that works pretty fast.

## Files

| File | Description |
| --- | --- |
| `numinteg_brob.r` | numerical integration functions for Brobdingnag (from Bolstad2: `sintegral`, from pracma: `romberg`, `cotes`, `kronrod`, `quadgk`, `simpadpt`, `quadv`, `quadinf`, `quad`, `trapz`, `quadgr`) |
| `numinteg_brob_calls.r` | example calls |
| `numinteg_log.r` | numerical integration functions for log functions (rom cmna: `trapez`, `romberg`) |
| `numinteg_log_calls.r` | example calls |
| `simpsonrule_nlb.r` | general simpson rule for normal, log, and Brobdingnag numbers |
| `simpsonrule_nlb_calls.r` | example calls |

## Functions used

Original names just got a `*.brob` added at the end of the filename.

| Function | Description |
| --- | --- |
| `numinteg_helper-funcs.r` | helper functions for brob objects and log calculations |
| `list2vec.brob` | convert a brob list to a brob vector |
| `scalarprod.brob` | replace `%*%` scalarproduct that does not work for brob objects |
| `.llog2sub.short` | `log(x-y)` |
| `.llog.2add.shor` | `log(x+y)` |
|  |  |
| `numinteg_brob.r` | for brob objects |
| `sintegral.brob.parallel` | taken from [Bolstad2](https://github.com/cran/Bolstad2) |
|  | taken from [pracma](https://github.com/cran/pracma): |
| `romberg.brob` | Romberg |
| `cotes.brob` | Cotes |
| `integral.brob` | Kronrod, Simpson - but not Clenshaw |
| `quadgk.brob` | internal function |
| `simpadpt.brob` | Simpson |
| `quadinf.brob` | internal function |
| `quad.brob` | internal function |
| `trapz.brob` | Trapez |
| `quadgr.brob` | internal function |
|  |  |
| `numinteg_log.r` | taken from [cmna](https://jameshoward.us/books/computational-methods-numerical-analysis-r) for log functions: |
| `trap.nl` | trapez for log |
| `romberg.nl` | romberg for log |
|  |  |
| `simpsonrule_nlb.r` | Simpson rule |
| `simpsonrule.nlb` | for normal, log functions, and brob objects |

Note: Not all integration functions of pracma are covered. What is missing are:

- Clenshaw from `integral`
- `integral2`
- `integral3`

The same is true for cmna.

## Function calls

Functions are called identical to the original version. Use as normal but add `.brob` to the name of the function call or as outlined above, and use a function that gives out a Brobdingnag object that acts as input for the integration function. For the log version you need a function that gives out a log value.

```
source("numinteg_helper-funcs.r")
source("numinteg_brob.r")
```

e.g.

```
fx <- function(x) sin(x)
f <- function(x) as.brob(sin(x))
lower <- 0
upper <- pi
maxit <- 25
tol <- 1e-15
romberg.brob(f, lower, upper, maxit=maxit, tol=tol)
as.brob(romberg(f=fx, a, b, maxit=maxit, tol=tol)$value)
```

The *_calls.r files contain examples for each function.

## Accurcacy

There is no special treatment to ensure accuracy embedded beyond what is already in the role models from pracma, cmna, and Bolstad2. If one is interested in that, you have to write your own tolerance error function based on known error functions of numerical integration methods. Another option and more pragmatic is to use functions for which the integration result can be obtained algebraically ie. via an analytical solution beyond any doubt. This can be compared to the output of the various functions. Beyond the numerical precision of the machine one should not trust numbers anyway.

## Notes

Using parallel for multi-threading does not lead to expected results - see for yourself. It may be helpful if a function to be integrated numerically is highly complicated and it requires a lot of repitiions. Otherwise, trials showed no advantage of using parallel over normal computing. Rather, the opposite was the case - the initial overhead of parallel required more time than normal. If this is repeated, the same story starts again - so no gain at all at this stage.

## Benchmarking

For benchmarking, the R package [microbenchmark](https://github.com/joshuaulrich/microbenchmark/) is a good choice and easy to use. See the manpage of how to use it and some example calls here.

## TODO

Cover all integration methods and add them as multi-threaded nevertheless to find out how to integrate parallel computing better and more efficiently to speed up the process.

## Licenses

See license file included - in short it is GPL >=2 - in accordance to the packages used and the license associated with them.

## R version

All R scripts should work under R >=v3.

## Disclaimer

The R code was tested carefully and cross-checked against various public available results and manpages to ensure proper results.  The example calls can be used to test for accuracy. However, it is provided "as is". Use common sense to compare results with expectations. NO WARRANTY of any kind is involved here. There is no guarantee that the software is free of error or consistent with any standards or even meets your requirements. Do not use the software or rely on it to solve problems if incorrect results may lead to hurting or injurying living beings of any kind or if it can lead to loss of property or any other possible damage to the world. If you use the software in such a manner, you are on your own and it is your own risk.


