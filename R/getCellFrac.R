# dt should be addAnno for binomial and genic for naive
#' Get Cell fraction
#'
#' No longer used
#'
#' @param dt A \code{data.table}. The output of \code{getGenicDP}.
#' @export
getCellFrac2 <- function(dt, xciGenes = NULL, method = "naive", xist_only = TRUE, rm_inac = 2, plot = FALSE){
  method <- tolower(method)
  xci <- readXCI(xciGenes)
  dt_xci <- dt[GENE %in% xci]
  if(!is.null(rm_inac)){ #Remove genes that are not perfectly inactivated in the 2 most skewed samples
    dt_xci <- rmInac(dt_xci, rm_inac)
  }
  if(method == "binomial"){
    dt_frac <- getBinomCellFrac(dt_xci)
  } else if(method == "naive"){
    dt_frac <- getNaiveCellFrac2(dt_xci)
  }
  dt_xist <- getXistFrac(dt, xist_only)
  xci_frac <- unique(merge(dt_frac[, list(sample, f, var_f, frac_median)], dt_xist[, list(sample, frac_XIST, var_frac_XIST)], by = "sample"))
  if(plot){
    plotCellFrac(dt_xci, dt_xist)
  }
  all_frac <- merge(dt, xci_frac, by = "sample")
  # Now get fraction for genes that are not perfectly inactivated
  all_frac[, fg := pmin(AD_hap1, AD_hap2)/(AD_hap1 + AD_hap2)]
  all_frac[, var_fg := fg*(1-fg)/tot]
  return(all_frac[])
}

getNaiveCellFrac <- function(dt_xci){
  dt_xci[, frac := pmin(AD_hap1, AD_hap2)/(AD_hap1 + AD_hap2)]
  dt_xci <- dt_xci[, frac_mean := mean(frac, na.rm = TRUE), by = "sample"]
  dt_xci <- dt_xci[, frac_median := median(frac, na.rm = TRUE), by = "sample"]
  dt_frac <- unique(dt_xci[, list(sample, frac_mean, frac_median)])
  return(dt_frac)
}

# Used for two sample t-test
getNaiveCellFrac2 <- function(dt_xci){
  dt_xci[, fg := pmin(AD_hap1, AD_hap2)/(AD_hap1 + AD_hap2)]
  dt_xci[, var_fg := fg*(1-fg)/tot]
  dt_xci[, Nxci := .N, by = "sample"]
  dt_xci[, f := mean(fg, na.rm = TRUE), by = "sample"]
  dt_xci[, var_f := sum(var_fg)/(Nxci^2), by = "sample"] #Only xci in table
  dt_xci[, frac_median := median(fg, na.rm = TRUE), by = "sample"]
  return(dt_xci)
}

# Find the fraction that maximizes the likelihood of observing the data, assuming the data follows a binomial distribution
# Here, the data is y = AD_hap1, n = tot and what we observe is the observed fractions for each SNP
getBinomCellFrac <- function(anno, trim = TRUE){
  # The annotation table, filtered for XCI genes
  anno[, frac := AD_hap1/(AD_hap1 + AD_hap2)]
  if(trim){
    anno[, `:=`(c("lo", "hi"), as.list(quantile(frac, probs = c(.05, .95)))), by = "sample"]
    anno <- anno[frac > lo & frac < hi]
  }
  anno[, frac_naive := mean((AD_hap1)/(AD_hap1+AD_hap2)), by = "sample"]
  anno[, mu := mean(frac > .5, na.rm = TRUE), by = "sample"] #mu is the proba that AD_hap1 is the most expressed (i.e: xa)
  for(sample_i in unique(anno$sample)){
    y <- anno[sample == sample_i, AD_hap1]
    n <- anno[sample == sample_i, AD_hap1 + AD_hap2]
    p <- unique(anno[sample == sample_i, frac_naive])
    mu <- unique(anno[sample == sample_i, mu])
    pars0 <- c(p, mu) #First in Newton Raphson
    res <- nlminb(pars0, log.likelihood.binom, dp1 = y, dp = n) #Find MLE
    anno[sample == sample_i, p_est := res$par[1]]
  }
  anno[p_est > .5, p_est := 1-p_est]
  anno[frac > .5, frac := 1-frac]
  anno <- anno[, frac_median := median(frac, na.rm = TRUE), by = "sample"]
  setnames(anno, "p_est", "frac_mean")
  dt_frac <- unique(anno[, list(sample, frac_mean, frac_median)])
  return(dt_frac)
}

log.likelihood.binom <- function(pars, dp1, dp){
  p <- pars[1]
  mu <- pars[2]
  l <- (-1)*sum(log(mu*p^dp1*(1-p)^(dp-dp1)+(1-mu)*p^(dp-dp1)*(1-p)^dp1))
  return(l)
}


#' @importFrom ggplot2 ggplot aes theme geom_point
#' @importFrom ggplot2 geom_text geom_hline scale_colour_manual facet_wrap
plotCellFrac <- function(xci_dt, xist_dt){
  plotfrac <-  unique(merge(xci_dt, xist_dt[, list(sample, frac_XIST)], by = "sample"))
  plotfrac <- plotfrac[order(sample)]
  Nt <- plotfrac[, .N, by = "sample"]
  plotfrac <- merge(plotfrac, Nt, by = "sample")
  plotfrac[, index := unlist(sapply(Nt[, N], function(x){1:x}))]
  plotfrac[, label := ""]
  plotfrac[abs(fg - f) > .2, label := GENE]
  p <- ggplot(plotfrac, aes(x = index, y = fg)) + geom_point() +
    geom_hline(aes(yintercept = f, colour = "mean")) +
    geom_hline(aes(yintercept = frac_median, colour = "median")) +
    geom_hline(aes(yintercept = frac_XIST, colour = "XIST"))
  p <- p + geom_text(aes(label = label))
  p <- p + geom_text(aes(x = N - 5, y= max(fg) + .01, label = paste("N =", N)))
  p <- p + facet_wrap(~sample, scales = "fixed")
  p <- p + scale_colour_manual(name = "fraction", values = c("red", "blue", "green"), labels = c("mean", "median", "xist")) + theme(legend.position = c(.8, 0.2))
  print(p)
  return(NULL)
}

#
getXistFrac <- function(dt, xist_only = FALSE){
  if(xist_only){
    dt_xist <- dt[GENE == "XIST"]
  } else{ #We add reads from all genes
    dt_xist <- dt[grep("XIST", GENE),]
    dt_xist[, `:=`(c("n_snps", "AD_hap1", "AD_hap2", "tot"), list(sum(n_snps), sum(AD_hap1), sum(AD_hap2), sum(tot))), by = "sample"]
    dt_xist[, `:=`(c("REF", "ALT", "GENE"), list("X", "X", "XIST"))]
    dt_xist <- unique(dt_xist[, POS := NULL])
  }
  dt_xist[, frac_XIST := sum(pmin(AD_hap1, AD_hap2))/sum((AD_hap1 + AD_hap2)), by = "sample"]
  dt_xist[, var_frac_XIST := fracXIST*(1-frac_XIST)/tot]
  return(dt_xist)
}
