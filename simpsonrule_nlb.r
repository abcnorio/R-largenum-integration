
################################################################################
# Simpson rule for normal, log, and brob
simpsonrule.nlb <-  function(fx,
                             lower, upper,
                             method="normal", # log, brob
                             eps=.Machine$double.xmin,
                             #eps=.Machine$double.eps,
                             log1p.crit=1e-09, # 1e-15
                             Nsteps=100,
                             parallel=FALSE,
                             fixINF=TRUE, ...)
{
  
  if(lower == upper)
  {
    cat("\nLower and upper limits of the integral are identical.\nNo calculations done.\n")
    return(0)
  } else if(lower > upper)
  {
    cat("\nBe aware of the limits chosen:\tlower > upper\n.Negative value of the integral will result.\n")
    return( -1 * simpsonrule.nlb(fx,
                                 upper, lower,
                                 method=method,
                                 eps=eps,
                                 log1p.crit=log1p.crit,
                                 Nsteps=Nsteps,
                                 parallel=parallel,
                                 fixINF=fixINF) )
  }    
  
  sek.l <- 2*Nsteps+1
  hml.s <- sign(upper-lower)
  
  # avoid infinite values (zero with log(0), etc.)
  if(fixINF == TRUE)
  {
    if(method %in% c("normal","log"))
    {
      if(is.infinite(fx(lower))) lower <- lower + eps * hml.s
      if(is.infinite(fx(upper))) upper <- upper - eps * hml.s      
    } else if(method == "brob")
    {
      if(is.infinite(fx(lower)@x)) lower <- lower + eps * hml.s
      if(is.infinite(fx(upper)@x)) upper <- upper - eps * hml.s
    }
  }
  sek <- seq(lower,upper,length=sek.l)
  
  if(method == "normal")
  {
    xfxap <- fx(sek)
    h <- xfxap[2] - xfxap[1]
    res <- sum(h*( xfxap[2 * (1:Nsteps) - 1] +
                     4*xfxap[2 * (1:Nsteps)] +
                     xfxap[2 * (1:Nsteps) + 1] )
               /3)
  } else if(method == "log")
  {
    if(parallel)
    {
      require(parallel)
      mccores <- detectCores(all.tests=FALSE, logical=TRUE)
      fx.log <- simplify2array(mclapply(seq_along(1:sek.l),
                                        function(i) fx(sek[i]),
                                        mc.cores=mccores
      ) )
    } else
    {
      fx.log <- Vectorize(fx)(sek)
    }
    h.log <- log(sek[2] - sek[1])
    s.log <- c( fx.log[2 * (1:Nsteps) - 1],
                log(4) + fx.log[2 * (1:Nsteps)],
                fx.log[2 * (1:Nsteps) + 1]
    )
    max.s.log <- max(s.log)
    sum.log <- sum(exp(s.log - max.s.log))
    #res <- h.log - log(3) + max.s.log + log1p(sum(exp(s.log[-which.max(s.log)] - max.s.log)))
    #
    # when do log1p(x) vs. log(x)
    # https://scicomp.stackexchange.com/questions/20629/when-should-log1p-and-expm1-be-used
    # IF 1+x=1 in floating point accuracy
    # in case of x = 1e-15 which is -34.53878 on log scale
    # ie. usual double precision arithmetic
    # if is already useful for
    # x = 1e−09 which is -20.72327 on log scale
    ifelse(sum.log < log1p.crit,
           res <- h.log - log(3) +  max.s.log + log1p(sum.log),
           res <- h.log - log(3) +  max.s.log + log(sum.log)
    )
    
  } else if(method=="brob")
  {
    fx.brob <- list2vec.brob(lapply(seq_along(1:sek.l), function(i) fx(sek[i])))
    h <- as.brob(sek[2] - sek[1]) # diff
    res <- sum(h/as.brob(3) * ( fx.brob[2 * (1:Nsteps) - 1] +
                                  4 * fx.brob[2 * (1:Nsteps)] +
                                  fx.brob[2 * (1:Nsteps) + 1]
    ) )
  }
  return(res)
}
################################################################################
