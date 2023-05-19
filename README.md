# R code to enable integration functions for very large number

## Overview

The Brobdingnag R package provides functions to enable the work with very large numbers. Most integration functions fail if large numbers temporarily occur. The R code here implements functions of the R package Brobdingnag to enable integration functions from R packages pracma and Bolstad2 to work with very large numbers.

## License

see license file included - in short it is GPL - in accordance to the packages used.

## Files

- **brobdingnag.integral.r** = tweaked functions to work with very large numbers
- **brobdingnag.integral_calls** = example code that shows that the tweaked functions produce the same results as the original functions

## Functions used - original names just got a *.brob added:

- list2vec.brob = convert a Brobdingnag list to a Brobdingnag vector
- scalarprod.brob = replace '%*%' scalarproduct that does not work for Brobdingnag objects

from Bolstad2:

- sintegral.brob.parallel

from pracma:

- romberg.brob
- cotes.brob
- integral.brob ("Kronrod", "Simpson" but not "Clenshaw")
- quadgk.brob
- simpadpt.brob
- quadinf.brob
- quad.brob
- trapz.brob
- quadgr.brob

Not all integral functions of pracma are tweaked. What is missing are:

- Clenshaw from integral
- integal2
- integral3

## Function calls

Functions are called identical to the original version.

```
source("brobdingnag.integral.r")
[...use as always but add .brob to the function call]
```

## R version

All R scripts should work under R >=v3.

## Disclaimer

R code was tested but is provided "as is".


