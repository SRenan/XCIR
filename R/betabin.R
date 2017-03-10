#' @export
betaBinomXI <- function(genic_dt, singlecell = F, rm_inac = 2, method = "MM", plot = FALSE, xci = 1, hist=FALSE, type = "midp"){
  dt <- copy(genic_dt)
  dt[, dp1 := pmin(AD_hap1, AD_hap2)]
  dt[, tot := AD_hap1 + AD_hap2]

  # Use XCI to estimate parameters
  if(xci > 0){
    xcig <- readXCI()
  }
  #if(xci > 1){
  #  xcig <- xcig[! xcig %in% c(gNrm, "SEPT6")]
  #  xcig <- xcig[! xcig %in% c(gNrm)]
  #}
  #if(xci > 2){
  # dt_xci <- rmInac(dt_xci, rm_inac)
  #}
  #if(tolower(xci) == "j"){
  # xcig <- jcss[statusJ == "S", GENE]
  #}
  inactivated_genes <- xcig
  dt_xci <- dt[GENE %in% inactivated_genes]

  samples <- unique(dt_xci$sample)
  if(singlecell){
    balanced <- samples[!samples %in% getSkewedSamples(dt)]
  } else{
    balanced <- samples
  }

  if(method == "BB"){
    dt <- BB(dt_xci, dt, balanced)
    #for(sample_i in balanced){
    #  a0  <- c(1,1);
    #  dp1 <- dt_xci[sample == sample_i, dp1]
    #  dp  <- dt_xci[sample == sample_i, tot]
    #  res_optim <- nlminb(a0, logL_BB, dp1 = dp1, dp = dp)
    #  dt[sample == sample_i, a_est := exp(res_optim$par[1])]
    #  dt[sample == sample_i, b_est := exp(res_optim$par[2])]
    #}
  } else if(method == "MM"){
    for(sample_i in balanced){
      a0  <- c(1, 1, log(0.98/0.02), log(0.02/0.98))
      dp1 <- dt_xci[sample == sample_i, dp1]
      dp  <- dt_xci[sample == sample_i, tot]
      res_optim <- nlminb(a0, logL_MM, dp1 = dp1, dp = dp)
      dt[sample == sample_i, a_est := exp(res_optim$par[1])]
      dt[sample == sample_i, b_est := exp(res_optim$par[2])]
      dt[sample == sample_i, p_het := exp(res_optim$par[3])/(1 + exp(res_optim$par[3]))]
    }
  } else if(method == "MM2"){
    for(sample_i in balanced){
      a0 <- c(1, 1, #Beta for the inactivated genes mean = a/(a+b) = .5
              2, 2, #Beta for the escape  #2 has more density around 0.5
              log(0.9/0.1)) #Initial proba that any XCIg is actually escaped in this sample
      dp1 <- dt_xci[sample == sample_i, dp1]
      dp  <- dt_xci[sample == sample_i, tot]
      res_optim <- nlminb(a0, logL_MM2, dp1 = dp1, dp = dp)
      p_inac_i <- exp(res_optim$par[5])/(1 + exp(res_optim$par[5]))
      if(p_inac_i >= .5){ # We assume that the list should always have a majority of inactivated genes
        a_est_i <- exp(res_optim$par[1])
        b_est_i <- exp(res_optim$par[2])
      } else{
        p_inac_i <- 1-p_inac_i
        a_est_i <- exp(res_optim$par[3])
        b_est_i <- exp(res_optim$par[4])
      }
      dt[sample == sample_i, a_est := a_est_i]
      dt[sample == sample_i, b_est := b_est_i]
      dt[sample == sample_i, p_inac := p_inac_i]
    }
  } else if(method == "MM3"){
    for(sample_i in balanced){
      a0 <- c(1, 1, #Beta for the inactivated genes mean = a/(a+b) = .5
              2, 2, #Beta for the escape  #2 has more density around 0.5
              log(0.9/0.1), #Initial proba that any XCIg is actually escaped in this sample
              log(0.98/0.02), #Proba of SNP being actually het
              log(0.02/0.98))
      dp1 <- dt_xci[sample == sample_i, dp1]
      dp  <- dt_xci[sample == sample_i, tot]
      res_optim <- nlminb(a0, logL_MM3, dp1 = dp1, dp = dp)
      p_inac_i <- exp(res_optim$par[5])/(1 + exp(res_optim$par[5]))
      if(p_inac_i >= .5){ # We assume that the list should always have a majority of inactivated genes
        a_est_i <- exp(res_optim$par[1])
        b_est_i <- exp(res_optim$par[2])
      } else{
        p_inac_i <- 1-p_inac_i
        a_est_i <- exp(res_optim$par[3])
        b_est_i <- exp(res_optim$par[4])
      }
      dt[sample == sample_i, a_est := a_est_i]
      dt[sample == sample_i, b_est := b_est_i]
      dt[sample == sample_i, p_inac := p_inac_i]
      dt[sample == sample_i, p_het := exp(res_optim$par[6])/(1+exp(res_optim$par[6]))]
    }
  } else{
    stop("Invalid method selected")
  }

  dt[, f := a_est/(a_est + b_est)]
  dt[, fg := dp1/tot]

  # Use the estimated a and b to calculate var_fg and the test statistic
  dt[, var_fg := (tot * a_est * b_est * (a_est + b_est + tot))]
  dt[, var_fg := var_fg/((a_est + b_est)^2 * (a_est + b_est + 1))]
  dt[, var_fg := var_fg/tot^2] # Because we want the variance for the fraction, not the variance for the counts
  dt[, t := (fg-f)/sqrt(var_fg)]
  dt[, pbb := pbb(dp1, tot, a_est, b_est, type = type), by = c("GENE", "sample")]
  dt[, p_value := pnorm(t, lower.tail = F)]

  if(hist){
    hist(unique(dt[, f]), breaks = seq(0, .5, .01), main = "Cell fraction", xlab="Cell fraction", ylab = "samples")
  }
  if(plot){
    plotBBCellFrac(dt[GENE %in% xcig])
  }
  return(dt[])
}

################################################################################
dbb <- function(x, n, a, b){#Beta binomial
    #y <- choose(n, x) * beta(x+a, n-x+b)/beta(a,b)
    c <- choose(n, x)
    bn <- beta(x+a, n-x+b)
    bd <- beta(a, b)
    y <- c * bn / bd
    y[is.nan(y)] <- 0
    return(y);
}

# P-values for exact inference
pbb <- function(x, n, a, b, type = "gt"){ #p-value for exact inference
  if(is.na(x) | is.na(n)){
    return(NA)
  }
  if(x >= n)
    return(0)
  if(type == "gt"){
    from <- x+1
    sump <- sum(dbb(from:n, n, a, b), na.rm = T)
  } else if(type == "geq"){
    from <- x
    sump <- sum(dbb(from:n, n, a, b), na.rm = T)
  } else if(type == "midp"){
    from <- x+1
    t0 <- dbb(x, n, a, b)
    t1 <- sum(dbb(from:n, n, a, b), na.rm = T)
    sump <- .5*t0 + t1
  } else{
    stop("Type should be one of 'gt', 'geq', 'midp'")
  }
  return(sump)
}

################################################################################
# Beta-binomial model:

BB <- function(dt_xci, full_dt, balanced){
  dt <- copy(full_dt)
  for(sample_i in balanced){
    a0  <- c(1,1);
    dp1 <- dt_xci[sample == sample_i, dp1]
    dp  <- dt_xci[sample == sample_i, tot]
    res_optim <- nlminb(a0, logL_BB, dp1 = dp1, dp = dp)
    dt[sample == sample_i, a_est := exp(res_optim$par[1])]
    dt[sample == sample_i, b_est := exp(res_optim$par[2])]
  }
  return(dt)
}
logL_BB <- function(a, dp1, dp){
  a1 <- exp(a[1])
  b1 <- exp(a[2])
  logLbb <- sum(log(dbb(dp1,dp,a1,b1))*(-1));
  return(logLbb);
}
################################################################################
# Mixture models:

# 1 Binom for sequencing errors and 1 BB for inactivated heterozygous SNP
logL_MM <- function(a, dp1, dp) {
  a1 <- exp(a[1]);
  b1 <- exp(a[2]);
  p_het <- exp(a[3])/(1+exp(a[3])); # Proba that the SNP is indeed bi-allelic (initial proba = .98)
  pi_err <- 0.01+0.99*exp(a[4])/(1+exp(a[4])); # (initial proba = .02)
  pbb <- dbb(dp1, dp, a1, b1)
  pb <- dbinom(dp1, dp, pi_err)
  p_tot <- p_het * pbb + (1-p_het) * pb
  logL  <- -sum(log(p_tot))
  return(logL);
}

# 1 BB for inactivated SNP and 1 BB for escaped SNP
logL_MM2 <- function(a, dp1, dp) {
  # Inactivated gene from XCI
  ai <- exp(a[1]);
  bi <- exp(a[2]);
  # Escape gene in XCI list
  ae <- exp(a[3])
  be <- exp(a[4])

  p_inac <- exp(a[5])/(1+exp(a[5])) #Proba that the team is indeed inactivated
  pbbi <- dbb(dp1, dp, ai, bi)
  pbbe <- dbb(dp1, dp, ae, be)
  p_tot <- p_inac * pbbi + (1-p_inac) * pbbe
  logL  <- -sum(log(p_tot))
  return(logL);
}

# 1 BB for inac SNPs 1 BB for escaped SNP and 1 Binom for sequencing err
logL_MM3 <- function(a, dp1, dp){
  # Inactivated gene from XCI
  ai <- exp(a[1])
  bi <- exp(a[2])
  # Escape gene in XCI list
  ae <- exp(a[3])
  be <- exp(a[4])

  p_inac <- exp(a[5])/(1+exp(a[5])) #Proba that the team is indeed inactivated
  pbbi <- dbb(dp1, dp, ai, bi)
  pbbe <- dbb(dp1, dp, ae, be)

  p_het <- exp(a[6])/(1+exp(a[6]))
  pi_err  <- .01 + .99*exp(a[7])/(1+exp(a[7]))

  pbbi <- dbb(dp1, dp, ai, bi)
  pbbe <- dbb(dp1, dp, ae, be)
  pb <- dbinom(dp1, dp, pi_err)
  p_tot <- p_het * (p_inac * pbb_i + # Reads from component1
                   (1-p_inac) * pbb_e) + # Reads from component2
           (1-p_het) * pb_err # Reads from sequencing error

  logL <- -sum(log(p_tot))
  return(logL)
}

################################################################################

# Can also be used to visualize the cell frac from the calls
#' Plot cell fraction estimates
#'
#' Plot cell fraction estimates from list of known XCI genes
#'
#' @details
#' This function is mostly used in \code{betaBinomXI} to ensure that the cell
#' fraction is estimated properly. However, it can be used from the output
#' of \code{betaBinomXI} to troubleshoot estimation issues.
#'
#' @export
plotBBCellFrac <- function(xci_dt, all_names = F){
  plotfrac <- xci_dt[order(sample)]
  Nt <- plotfrac[, .N, by = "sample"]
  plotfrac <- merge(plotfrac, Nt, by = "sample")
  plotfrac[, index := unlist(sapply(Nt[, N], function(x){1:x}))]
  plotfrac[, label := ""]
  if(all_names){
    plotfrac[, label := GENE]
  }else{
    plotfrac[abs(fg - f) > .2, label := GENE]
  }
  p <- ggplot(plotfrac, aes(x = index, y = fg)) + geom_point() +
    geom_hline(aes(yintercept = f, colour = "mean"))
  p <- p + geom_text(aes(label = label))
  p <- p + geom_text(aes(x = N - 5, y= max(fg) + .01, label = paste("N =", N)))
  p <- p + facet_wrap(~sample, scales = "fixed")
  p <- p + scale_colour_manual(name = "fraction", values = c("red"), labels = c("mean")) + theme(legend.position = c(.8, 0.2))
  print(p)
  return(NULL)
}
