M0cov <- function(dt_xci, full_dt, a0 = NULL, optimizer = "nlminb", method = NULL,
               limits = FALSE, debug = FALSE){
  dt <- copy(full_dt)
  samples <- unique(dt_xci$sample)
  if(is.null(a0)){
    a0  <- c(1,1)
  }
  for(sample_i in samples){
    if(debug)
      print(sample_i)
    dp1 <- dt_xci[sample == sample_i, dp1]
    dp  <- dt_xci[sample == sample_i, tot]
    Nxcig <- dt_xci[sample == sample_i, .N]
    if(limits){
      lowb <- c(0, 0)
      upb <- c(Inf, Inf)
      upb <- c(699, 700)
      method <- "L-BFGS-B"
    } else{
      lowb <- -Inf
      upb <- Inf
    }
    if(optimizer == "nlminb"){
      res_optim <- nlminb(a0, .logL_M0_cov, dp1 = dp1, dp = dp,
                          lower = lowb, upper = upb)
      negLogLname <- "objective"
    } else if(optimizer == "optim"){
      res_optim <- optim(par = a0, fn = .logL_M0_cov, method =  method,
                         dp1 = dp1, dp = dp,
                         lower = lowb, upper = upb)
      negLogLname <- "value"
    } else{
      stop("Unknown optimizer. Should be one of 'nlminb', 'optim'")
    }
    message <- res_optim$message
    message <- ifelse(is.null(message), "", message)
    dt[sample == sample_i, a_est := exp(res_optim$par[1]) + 1]
    dt[sample == sample_i, b_est := exp(res_optim$par[2]) + 1]
    dt[sample == sample_i, model := "M0cov"]
    dt[sample == sample_i, k := length(res_optim$par)]
    dt[sample == sample_i, logL := res_optim[[negLogLname]]] #-logL
    dt[sample == sample_i, convergence := res_optim$convergence]
    dt[sample == sample_i, flag := message]
    dt[sample == sample_i, AIC := 2*k + 2*logL]
    dt[sample == sample_i, Ntrain := Nxcig]
    dt[sample == sample_i, BIC := log(Ntrain)*k + 2*logL]
  }
  return(dt)
}
.logL_M0_cov <- function(param, dp1, dp){
  # Handling covariates is about how we obtain a1,b1
  # From a binomial distribution, I can obtain mu using mu = 1/(1+exp(-b0 +b1 X))
  # Can I use var = mu(1-mu)? I think I lose information compared to the BB since I lose 1 param.
  # Alternative is to obtain a measure of variance or the overdispersion parameter.
  # Common code
  a1 <- exp(a[1]) + 1
  b1 <- exp(a[2]) + 1
  logLbb <- sum(ldbb(dp1,dp,a1,b1)*(-1), na.rm = TRUE)
  return(logLbb)
}
