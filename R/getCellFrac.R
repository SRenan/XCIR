# dt should be addAnno for binomial and genic for naive
#' @importFrom ggplot2 geom_text geom_hline scale_colour_manual facet_wrap
getCellFrac2 <- function(dt, xciGenes = NULL, method = "binomial", xist_only = TRUE, rm_inac = 2, plot = F){
  xci <- readXCI(xciGenes)
  dt_xci <- dt[GENE %in% xci]
  if(!is.null(rm_inac)){
    if(!is.numeric(rm_inac))
      stop("rm_inac should be null or numeric")
    notinactivated <- unique(c(dt_xci[grep("LC2$", sample)][AD_hap1 > rm_inac & AD_hap2 > rm_inac, GENE],
                               dt_xci[grep("LC3$", sample)][AD_hap1 > rm_inac & AD_hap2 > rm_inac, GENE]))
    dt_xci <- dt_xci[!GENE %in% notinactivated]
    print(paste("Using", length(unique(dt_xci[, GENE])), "inactivated genes to calculate cell fraction."))
  }
  if(tolower(method) == "binomial"){
    dt_frac <- getNaiveCellFrac(dt_xci)
  } else if(tolower(method) == "naive"){
    dt_frac <- getNaiveCellFrac(dt_xci)
  }
  if(xist_only){
    dt_xist <- dt[GENE == "XIST"]
  } else{ #We add reads from all genes
    dt_xist <- dt[grep("XIST", GENE),]
    dt_xist[, `:=`(c("n_snps", "AD_hap1", "AD_hap2", "tot"), list(sum(n_snps), sum(AD_hap1), sum(AD_hap2), sum(tot))), by = "sample"]
    dt_xist[, `:=`(c("REF", "ALT", "GENE"), list("X", "X", "XIST"))]
    dt_xist <- unique(dt_xist[, POS := NULL])
  }
  dt_xist[, frac_XIST := pmin(AD_hap1, AD_hap2)/tot, by = "sample"]
  cellfrac <- merge(dt_frac, dt_xist[, list(sample, frac_XIST)], by = "sample")
  if(plot){
    plotfrac <-  merge(dt_xci, dt_xist[, list(sample, frac_XIST)], by = "sample")
    plotfrac <- plotfrac[order(sample)]
    Nt <- plotfrac[, .N, by = "sample"]
    plotfrac <- merge(plotfrac, Nt, by = "sample")
    plotfrac[, index := unlist(sapply(Nt[, N], function(x){1:x}))]
    p <- ggplot(plotfrac, aes(x = index, y = frac)) + geom_point() +
      geom_hline(aes(yintercept = frac_mean, colour = "red")) +
      geom_hline(aes(yintercept = frac_median, colour = "blue")) +
      geom_hline(aes(yintercept = frac_XIST, colour = "green"))

    p <- p + geom_text(aes(x = N - 5, y= max(frac) + .01, label = paste("N =", N)))
    p <- p + facet_wrap(~sample, scales = "fixed")
    p <- p + scale_colour_manual(name = "fraction", values = c("red", "blue", "green"), labels = c("mean", "median", "xist")) + theme(legend.position = c(.8, 0.2))
    print(p)
  }
  return(cellfrac)
}

getNaiveCellFrac <- function(dt_xci){
  dt_xci[, frac := pmin(AD_hap1, AD_hap2)/tot]
  dt_xci <- dt_xci[, frac_mean := mean(frac, na.rm = TRUE), by = "sample"]
  dt_xci <- dt_xci[, frac_median := median(frac, na.rm = TRUE), by = "sample"]
  dt_frac <- unique(dt_xci[, list(sample, frac_median, frac_mean)])
  return(dt_frac)
}

# Pars: c(logit of naive fraction, logit of mu)
# dp1:  AD_hap1
# dp:   tot
log.likelihood.binom <- function(pars, dp1, dp){
  logit.p <- pars[1]
  logit.mu <- pars[2]
  p  <- exp(pars[1])/(1+exp(pars[1]))
  mu <- exp(pars[2])/(1+exp(pars[2]))
  l  <- (-1)*sum(log(mu*p^dp1*(1-p)^(dp-dp1)+(1-mu)*p^(dp-dp1)*(1-p)^dp1))
  return(l)
}

log.likelihood.betabinom <- function(pars, dp1, dp){
  a <- exp(pars[1])
  b <- exp(pars[2])
  mu <- exp(pars[3])/(1+exp(pars[3]))
  l <- (-1)*sum(log(mu*dbb(dp1,dp,a,b)+(1-mu)*dbb(dp-dp1,dp,a,b)))
  return(l)
}
