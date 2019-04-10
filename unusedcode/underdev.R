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




