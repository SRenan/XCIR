# This file for functions under development before they get added to the files where they belong

# From XCIR output (skewing estimates) make calls for new SNPs
newCalls <- function(xcirout, newdata){
  clean <- XCIR:::sample_clean(xcirout)
  sscols <- c("sample", "CHROM", "POS", "GENE", "ANNO", "AD_hap1", "AD_hap2")
  newdata <- merge(clean, newdata[, sscols, with = F], by = "sample")
  newdata[, dp1 := pmin(AD_hap1, AD_hap2)]
  newdata[, tot := AD_hap1 + AD_hap2]

  # Estimates
  newdata[, fg := dp1/tot]
  newdata[, var_fg := (tot * a_est * b_est * (a_est + b_est + tot))]
  newdata[, var_fg := var_fg/((a_est + b_est)^2 * (a_est + b_est + 1))]
  newdata[, var_fg := var_fg/tot^2] # Because we want the variance for the fraction, not the variance for the counts
  newdata[, t := (fg-f)/sqrt(var_fg)] #Test statistic
  newdata[, p_value := pnorm(t, lower.tail = F)]

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
    out[, f_value := pchisq(X, df = 2*kcomb, lower.tail = F)]
  } else if(method == "stouffer"){
    out[, Z := sum(qnorm(1 - p_value))/sqrt(kcomb), by = c("sample", "GENE")]
    out[, z_value := pnorm(Z, 0, kcomb, lower.tail = F)]
  } else if(method == "zscore"){
    # Use total read count to weight
    out[, Z := sum(sqrt(tot)*qnorm(1 - p_value))/sqrt(sum(tot)), by = c("sample", "GENE")]
    out[, z_value := pnorm(Z, lower.tail = F)]
  }
  return(out)
}


getZ <- function(xcirout){
  Zdt <- copy(xcirout)
  Zdt[, mup := mean(p_value), by = "GENE"]
  Zdt[, sdp := sd(p_value), by = "GENE"]
  Zdt[, Z := (p_value - mup)/sdp]
}


### Unconstrained optimization with new parameter transformation
.logL_BB_unc <- function(a, dp1, dp){
  a1 <- exp(a[1])+1
  b1 <- exp(a[2])+1
  logLbb <- sum(ldbb(dp1,dp,a1,b1)*(-1), na.rm = T)
  return(logLbb)
}
.logL_MM_unc <- function(a, dp1, dp) {
  a1 <- exp(a[1])+1
  b1 <- exp(a[2])+1;
  # Proportions have to be kept as an exponent ratio to avoid issues with boundaries
  p_het <- exp(a[3])/(1+exp(a[3])); # Proba that the SNP is indeed bi-allelic (initial proba = .98)
  pi_err <- 0.001+0.999*exp(a[4])/(1+exp(a[4])); # (initial proba = .02)
  pb <- dbinom(dp1, dp, pi_err)

  # Operations on the log scale to avoid infinite values
  lpbb <- ldbb(dp1, dp, a1, b1)
  p_tot <- p_het * exp(lpbb) + (1-p_het)*pb

  logL  <- -sum(log(p_tot))
  return(logL);
}
.logL_MM2_unc <- function(a, dp1, dp) {
  # Inactivated gene from XCI
  ai <- exp(a[1])+1
  bi <- exp(a[2])+1
  # Escape gene in XCI list
  ae <- exp(a[3])+1
  be <- exp(a[4])+1
  # Proba that the gene is indeed inactivated (0.9)
  p_inac <- exp(a[5])/(1+exp(a[5]))

  lpbbi <- ldbb(dp1, dp, ai, bi)
  lpbbe <- ldbb(dp1, dp, ae, be)
  p_tot <- p_inac * exp(lpbbi) + (1-p_inac) * exp(lpbbe)
  logL  <- -sum(log(p_tot))
  return(logL);
}
.logL_MM3_unc <- function(a, dp1, dp){
  # Inactivated gene from XCI
  ai <- exp(a[1])+1
  bi <- exp(a[2])+1
  # Escape gene in XCI list
  ae <- exp(a[3])+1
  be <- exp(a[4])+1

  p_inac <- exp(a[5])/(1+exp(a[5]))
  p_het <- exp(a[6])/(1+exp(a[6]))
  pi_err  <- .01 + .99*exp(a[7])/(1+exp(a[7]))

  lpbb_i <- ldbb(dp1, dp, ai, bi)
  lpbb_e <- ldbb(dp1, dp, ae, be)
  pb_err <- dbinom(dp1, dp, pi_err)

  p_tot <- p_het * (p_inac * exp(lpbb_i) + # Reads from component1
                    (1-p_inac) * exp(lpbb_e)) + # Reads from component2
    (1-p_het) * pb_err # Reads from sequencing error
  logL <- -sum(log(p_tot))
  return(logL)
}



