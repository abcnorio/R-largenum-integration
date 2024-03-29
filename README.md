﻿# R code to enable numerical integration methods to work with very large numbers

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
| `numinteg_accuracy.r` | check accuracy + worked example |
| `numinteg_brob.r` | numerical integration methods for brob objects (from Bolstad2: `sintegral`, from pracma: `romberg`, `cotes`, `kronrod`, `quadgk`, `simpadpt`, `quadv`, `quadinf`, `quad`, `trapz`, `quadgr`) |
| `numinteg_brob_calls.r` | example calls |
| `numinteg_log.r` | numerical integration methods for log functions (from cmna: `trapez`, `romberg`) |
| `numinteg_log_calls.r` | example calls |
| `simpsonrule_nlb.r` | general simpson rule for normal, log, and brob objects |
| `simpsonrule_nlb_calls.r` | example calls |

## Methods used

Original names just got a `*.brob` added at the end of the filename.

| Function | Description |
| --- | --- |
| `numinteg_accuracy.r` | **Check accuracy + worked example** |
| `ck.accuracy` | check accuracy |
| `numinteg_helper-funcs.r` | **Helper functions for brob objects and log calculations** |
| `list2vec.brob` | convert a brob list to a brob vector |
| `scalarprod.brob` | replace `%*%` scalarproduct that does not work for brob objects |
| `.llog2sub.short` | `log(x-y)` |
| `.llog.2add.shor` | `log(x+y)` |
|  |  |
| `numinteg_brob.r` | **Numerical integration methods for brob objects** |
| `sintegral.brob.parallel` | taken from [Bolstad2](https://github.com/cran/Bolstad2) |
|  | taken from [pracma](https://github.com/cran/pracma): |
| `romberg.brob` | Romberg |
| `cotes.brob` | Cotes |
| `integral.brob` | Kronrod, Simpson - but not Clenshaw |
| `quadgk.brob` | adaptive Gauss-Kronrod quadrature |
| `simpadpt.brob` | Simpson |
| `quadinf.brob` | infinite integrals |
| `quad.brob` | adaptive Simpson quadrature |
| `trapz.brob` | Trapez |
| `quadgr.brob` | Gaussian quadrature with Richardson extrapolation |
|  |  |
| `numinteg_log.r` | **Numerical integration methods for log functions:** |
|  | taken from [cmna](https://jameshoward.us/books/computational-methods-numerical-analysis-r): |
| `trap.nl` | Trapez |
| `romberg.nl` | Romberg |
|  |  |
| `simpsonrule_nlb.r` | **Simpson rule** |
| `simpsonrule.nlb` | Simpson for normal, log functions, and brob objects |

Note: Not all integration methods of pracma are covered. What is missing are:

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

## Accuracy

There is no special treatment to ensure accuracy embedded beyond what is already in the role models from pracma, cmna, and Bolstad2. If one is interested in that, you have to write your own tolerance error function based on known error functions of numerical integration methods. Another option and more pragmatic is to use functions for which the integration result can be obtained algebraically ie. via an analytical solution beyond any doubt. This can be compared to the output of the various nethods. Beyond the numerical precision of the machine one should not trust numbers anyway unless accuracy is a dedicated feature of the software and emulates accuracy. The file `numinteg_accuracy.r` contains a simple worked example.

```
# function to integrate algebraically and numerically
f <- function(x) sin(x)
# integral analytically solved
# sin(x) dx -> -cos(x) + constant
# integration boundaries
lower <- 0
upper <- pi
f.integ <- function(lower,upper) -cos(upper) -(-cos(lower))
# real value
int.real.v <- f.integ(lower=lower, upper=upper)
# numerical integration using method on log scale
f.log <- function(x) log(sin(x))
simpson.l.v <- exp(simpsonrule.nlb(f.log, lower, upper, method="log", Nsteps=Nsteps))
# check accuracy
ck.accuracy(comp=int.real.v, numinteg=simpson.l.v, methods=c("real value","simpson log"))
```

This results in:

```
Checking Numerical Integration Methods
--------------------------------------------
comparative method	= real value
num. int. method	= simpson log
comparative value	= 2
num. integ. value	= 2.0000000000000000
diff comp-num.integ.	= 0.0000000000000000
comp. > num.integ.	= FALSE
abs(diff) < tol = 1e-08	= TRUE
factor comp/num.integ.	= 1.0000000000000000

---
Note: 'comparative' value can be a real value or from a different method
```

And now you can compare output results to a definite value to measure accuracy within the machine tolerance boundaries. Or you can compare two numerical integration methods. If you don't know whether your function can be solved algebraically, you can try with the online integration solver of wolframalpha.

## Benchmarking (and accuracy)

For benchmarking, the R package [microbenchmark](https://github.com/joshuaulrich/microbenchmark/) is a good choice and easy to use. See the manpage of how to use it and some example calls here. As a short example we use the code from the accuracy check above:

```
> # benchmarking
> library(microbenchmark)
> Nsteps <- 1e3
> numinteg.benchm <- microbenchmark(
+           integrate(f, lower, upper)$v,
+           Rmpfr:::integrateR(f.Rmpfr, lower, upper)$value,
+           pracma:::trapz(x=sek, y=f(sek)),
+           simpsonrule.nlb(f, lower, upper, method="normal", Nsteps=Nsteps),
+           exp(simpsonrule.nlb(f.log, lower, upper, method="log", Nsteps=Nsteps)),
+           as.numeric( trapz.brob(x=sek, y=f.sek.brob.c) )
+           )
> numinteg.benchm
Unit: microseconds
                                                                       expr       min        lq        mean    median        uq      max neval  cld
                                               integrate(f, lower, upper)$v    17.950    26.165    42.94091    46.615    55.640    71.14   100 a   
                            Rmpfr:::integrateR(f.Rmpfr, lower, upper)$value 25575.446 26318.942 28044.54153 26948.197 28003.902 39713.58   100  b  
                                        pracma:::trapz(x = sek, y = f(sek))   421.650   436.935   465.50396   445.700   454.520   829.30   100 a   
       simpsonrule.nlb(f, lower, upper, method = "normal", Nsteps = Nsteps)    62.140    71.510    90.17502    91.560    98.245   269.67   100 a   
 exp(simpsonrule.nlb(f.log, lower, upper, method = "log", Nsteps = Nsteps))  1403.640  1477.560  1687.81533  1521.095  1572.291 15313.48   100   c 
						  as.numeric(trapz.brob(x = sek, y = f.sek.brob.c))  7607.822  7882.082  8339.16269  8047.017  8205.902 21620.37   100    d 
```

**Note:**
Different methods are mixed here (e.g. Simpson rule, Trapez, Romberg (`integrateR`), adaptive quadrature (`integrate`)). Additionally, the number of steps `Nsteps` differs between `simpsonrule.nlb`, `trapz`, `integrate`, and `integrateR` esp. if a cut-off criterium is chosen based on some tolerance value. Thus, causes of any differences in computer time are difficult to identify precisely. For serious comparison one should standardize all influences so that a difference is indeed caused by either the input method (like `log)()` or `as.brob()`/ `brob()`) and not by the algorithm itself (like Simpson, Trapez, Romberg, etc.) or just by the number of steps `Nsteps` used to calculate the integral. This is valid for all other comparisons. As an ideal all influences should (must) be held constant except for what you want to compare. Beyond that, benchmarking is not everything -- accuracy must be considered too if you choose a method. Both together allow for a compromise if you work with very large numbers.

Thus, we can repeat the benchmarking above but we will use only one method `simpsonrule.nlb` but with different input: normal, log, brob. First, we check for accuracy:

```
> Nsteps <- 1e3
> # repeat for one method
> # add an alternative function
> f.brob2 <- function(x) brob(log(sin(x)))
>
> simpsonrule.nlb.n <- simpsonrule.nlb(f, lower, upper, method="normal", Nsteps=Nsteps)
> simpsonrule.nlb.l <- exp(simpsonrule.nlb(f.log, lower, upper, method="log", Nsteps=Nsteps))
> simpsonrule.nlb.b <- as.numeric(simpsonrule.nlb(f.brob, lower, upper, method="brob", Nsteps=Nsteps))
> simpsonrule.nlb.b2 <- as.numeric(simpsonrule.nlb(f.brob2, lower, upper, method="brob", Nsteps=Nsteps))
>
> ck.accuracy(comp=int.real.v, numinteg=simpsonrule.nlb.n, methods=c("real","simpsonrule.nlb - normal"))

Checking Numerical Integration Methods
--------------------------------------------
comparative method	= real
num. int. method	= simpsonrule.nlb - normal
comparative value	= 2
num. integ. value	= 1.9999991775331358
diff comp-num.integ.	= 0.0000008224668642
comp. > num.integ.	= TRUE
abs(diff) < tol = 1e-08	= FALSE
factor comp/num.integ.	= 1.0000004112336012

---
Note: 'comparative' value can be a real value or from a different method

> ck.accuracy(comp=int.real.v, numinteg=simpsonrule.nlb.l, methods=c("real","simpsonrule.nlb - log"))

Checking Numerical Integration Methods
--------------------------------------------
comparative method	= real
num. int. method	= simpsonrule.nlb - log
comparative value	= 2
num. integ. value	= 2.0000000000000675
diff comp-num.integ.	= -0.0000000000000675
comp. > num.integ.	= FALSE
abs(diff) < tol = 1e-08	= TRUE
factor comp/num.integ.	= 0.9999999999999662

---
Note: 'comparative' value can be a real value or from a different method

> ck.accuracy(comp=int.real.v, numinteg=simpsonrule.nlb.b, methods=c("real","simpsonrule.nlb - brob"))

Checking Numerical Integration Methods
--------------------------------------------
comparative method	= real
num. int. method	= simpsonrule.nlb - brob
comparative value	= 2
num. integ. value	= 2.0000000000000657
diff comp-num.integ.	= -0.0000000000000657
comp. > num.integ.	= FALSE
abs(diff) < tol = 1e-08	= TRUE
factor comp/num.integ.	= 0.9999999999999671

---
Note: 'comparative' value can be a real value or from a different method

> ck.accuracy(comp=int.real.v, numinteg=simpsonrule.nlb.b2, methods=c("real","simpsonrule.nlb - brob2"))

Checking Numerical Integration Methods
--------------------------------------------
comparative method	= real
num. int. method	= simpsonrule.nlb - brob2
comparative value	= 2
num. integ. value	= 2.0000000000000657
diff comp-num.integ.	= -0.0000000000000657
comp. > num.integ.	= FALSE
abs(diff) < tol = 1e-08	= TRUE
factor comp/num.integ.	= 0.9999999999999671

---
Note: 'comparative' value can be a real value or from a different method
```
 
and then we compare the speed:
 
```
> # just compare within one method but with different input
> numinteg.benchm2 <- microbenchmark(
+           simpsonrule.nlb(f, lower, upper, method="normal", Nsteps=Nsteps),
+           simpsonrule.nlb(f.log, lower, upper, method="log", Nsteps=Nsteps),
+           simpsonrule.nlb(f.brob, lower, upper, method="brob", Nsteps=Nsteps),
+           simpsonrule.nlb(f.brob2, lower, upper, method="brob", Nsteps=Nsteps)
+ )
> numinteg.benchm2
Unit: microseconds
                                                                   expr       min        lq         mean     median         uq       max neval cld
     simpsonrule.nlb(f, lower, upper, method = "normal", Nsteps = Nsteps)     61.64     66.79     81.06962     78.810     92.600    112.64   100 a  
    simpsonrule.nlb(f.log, lower, upper, method = "log", Nsteps = Nsteps)   1398.98   1485.38   1693.13454   1525.745   1609.541  12096.86   100 a  
  simpsonrule.nlb(f.brob, lower, upper, method = "brob", Nsteps = Nsteps) 231282.95 240144.51 252877.14441 251100.391 258555.633 515170.27   100  b 
 simpsonrule.nlb(f.brob2, lower, upper, method = "brob", Nsteps = Nsteps) 157890.80 167820.18 176322.08066 175793.131 183799.018 211417.95   100   c
 ```

**Note:**
As we can see the latter two ie. `f.brob` and `f.brob2` lead to exactly the same results, use the same calculation procedure but differ greatly in time. The only difference is `as.brob(x)` vs. `brob(log(x))` whereas the first one is much slower than the second one.

## Further notes

Using parallel for multi-threading does not lead to expected results -- see for yourself. It may be helpful if a function to be integrated numerically is highly complicated and it requires a lot of repetiions. Otherwise, trials showed no advantage of using parallel over normal computing. Rather, the opposite was the case: The initial overhead of parallel computing required more time than the normal time required to compute everything. If this is repeated, the same story starts again. So, there is no gain at all at this stage. Future efforts will show whether this can be changed.

## TODO

Cover all integration methods and add them as multi-threaded nevertheless to find out how to integrate parallel computing better and more efficiently to speed up the process.

## Licenses

See license file included - in short it is GPL >=2 - in accordance to the packages used and the license associated with them.

## R version

All R scripts should work under R >=v3.

## Disclaimer

The R code was tested carefully and cross-checked against various public available results and manpages to ensure proper results.  The example calls can be used to test for accuracy. However, it is provided "as is". Use common sense to compare results with expectations. NO WARRANTY of any kind is involved here. There is no guarantee that the software is free of error or consistent with any standards or even meets your requirements. Do not use the software or rely on it to solve problems if incorrect results may lead to hurting or injurying living beings of any kind or if it can lead to loss of property or any other possible damage to the world. If you use the software in such a manner, you are on your own and it is your own risk.


