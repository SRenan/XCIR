# This file for functions under development before they get added to the files where they belong

# From XCIR output (skewing estimates) make calls for new SNPs
newCalls <- function(xcirout, newdata){
  clean <- sample_clean(xcirout)
  sscols <- c("sample", "CHROM", "POS", "GENE", "ANNO", "AD_hap1", "AD_hap2")
  newdata <- merge(clean, newdata[, sscols, with = FALSE], by = "sample")
  newdata[, dp1 := pmin(AD_hap1, AD_hap2)]
  newdata[, tot := AD_hap1 + AD_hap2]

  # Estimates
  newdata[, fg := dp1/tot]
  newdata[, var_fg := (tot * a_est * b_est * (a_est + b_est + tot))]
  newdata[, var_fg := var_fg/((a_est + b_est)^2 * (a_est + b_est + 1))]
  newdata[, var_fg := var_fg/tot^2] # Because we want the variance for the fraction, not the variance for the counts
  newdata[, t := (fg-f)/sqrt(var_fg)] #Test statistic
  newdata[, p_value := pnorm(t, lower.tail = FALSE)]
  
  newdata[, status := ifelse(p_value < 0.05, "E", "S")]

  return(newdata)
}

#'
#' Combine p-values from multiple SNPs
#'
#' @param newcalls The \code{data.table} output by \code{newCalls}.
#' @param maxSNP A \code{numeric}. The maximum number of SNPs to use (i.e:
#'  p-values to combine).
#' @param method A \code{character}. One of "fisher", "stouffer", "zscore".
#'  zscore uses the total read count to weight each SNP.
#'
#' @return The \code{newcalls} \code{data.table} with a new f_value column.
#'
fisherCombine <- function(newcalls, maxSNP = 4, method = "fisher"){
  out <- copy(newcalls)
  out[, CNT := .N, by = c("sample", "GENE")]
  out[, IDX := seq_len(.N), by = c("sample", "GENE")]
  out <- out[IDX <= maxSNP]
  out[, kcomb := pmin(CNT, maxSNP)]
  if(method == "fisher"){
    out[, X := -2*sum(log(p_value)), by = c("sample", "GENE")]
    out[, f_value := pchisq(X, df = 2*kcomb, lower.tail = FALSE)]
  } else if(method == "stouffer"){
    out[, Z := sum(qnorm(1 - p_value))/sqrt(kcomb), by = c("sample", "GENE")]
    out[, z_value := pnorm(Z, 0, kcomb, lower.tail = FALSE)]
  } else if(method == "zscore"){
    # Use total read count to weight
    out[, Z := sum(sqrt(tot)*qnorm(1 - p_value))/sqrt(sum(tot)), by = c("sample", "GENE")]
    out[, z_value := pnorm(Z, lower.tail = FALSE)]
  }
  return(out)
}


getZ <- function(xcirout){
  Zdt <- copy(xcirout)
  Zdt[, mup := mean(p_value), by = "GENE"]
  Zdt[, sdp := sd(p_value), by = "GENE"]
  Zdt[, Z := (p_value - mup)/sdp]
}

# PHASING LL
# Limitations
# - Only works for two variants
# - Assumes most expressed variant is correct (s.t there are only two combinations)
.logL_BB_phasing <- function(a, dp1, dp2, dp){
  a1 <- exp(a[1]) + 1
  b1 <- exp(a[2]) + 1
  pp1 <- exp(a[3])/(1+exp(a[3])) #Probability that phasing is correct for the top two genes
  lpp1 <- XCIR:::ldbb(dp1, dp, a1, b1)
  lpp2 <- XCIR:::ldbb(dp2, dp, a1, b1)
  lik <- pp1 * exp(lpp1) + (1-pp1) * exp(lpp2)
  logL <- -sum(lik)
  return(logL)
}
# -> BB  needs to be changed to prep dp2 and pp1 (a[3])

# MAKING DATA TO PROCESS
# i.e: Top2 SNPs per gene


BB_phasing <- function(dt_xci, full_dt, a0 = NULL, optimizer = "nlminb", method = NULL,
                       limits = FALSE, debug = FALSE){
  dt <- copy(full_dt)
  samples <- unique(dt_xci$sample)
  if(is.null(a0)){
    a0  <- c(1,1, log(0.80/0.20))
  }
  for(sample_i in samples){
    if(debug)
      print(sample_i)
    dp1 <- dt_xci[sample == sample_i, dp1]
    dp2 <- dt_xci[sample == sample_i, dp2]
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
    # Refular optim
      res_optim <- nlminb(a0, XCIR:::.logL_BB, dp1 = dp1, dp = dp)
    # Refular optim
      res_optim2 <- nlminb(a0, .logL_BB_phasing, dp1 = dp1, dp2 = dp2, dp = dp)
    negLogLname <- "objective"
    
    message <- res_optim$message
    message <- ifelse(is.null(message), "", message)
    dt[sample == sample_i, a_est := exp(res_optim$par[1]) + 1]
    dt[sample == sample_i, b_est := exp(res_optim$par[2]) + 1]
    dt[sample == sample_i, model := "BB"]
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

# EM: Try using an EM instead of the likelihood based approach
M1 <- function(x, n, as, bs){
  # E step
  # The conditional probability to which subpopulation each observation x belongs is computed as the ratio
  # of the density computed for category C over the sum of all densities.
  
  # - Compute likelihood for each mixture
  # - Estimate the fraction of each mixture
  d0 <- pi0 * dbb(x, n, as, bs)
  d1 <- pi1 * dbinom(x, n, pi_err)
  d <- d0 + d1
  
  p0 <- d0/d
  p1 <- d1/d
  # M step
  # - Compute expected value of log likelihood
  
  pi0 <- sum(p0)/n
  pi1 <- sum(p1)/n
  return() #return the parameter values and possibly the number of iteration and the likelihood difference.
}