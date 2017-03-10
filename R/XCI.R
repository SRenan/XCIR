#' Get fraction of cell expressing the inactivated allele
#'
#' @param dt A \code{data.table}. The output of \code{getGenicDP}.
#' @param xciGenes A \code{character}. A list or a filename that contains a
#' list of known XCI genes.
#' @param fix_phasing A \code{logical}. If set to TRUE, the fraction above .5
#'  are set to 1-frac.
#' @param plot A \code{logical}. If set to TRUE, a plot of the median fraction
#'  for each sample will be generated.
#'
#' @details
#' The default xciGenes is a list of 177 genes known to be affected by X
#' chromosome inactivation from Cotton et al. 2013
#' (dx.doi.org/10.1186/gb-2013-14-11-r122).
#'
#' @return A \code{data.table} containing the average, and median fraction of
#' expression of the inactivated chromosome for each sample. frac_XIST is the
#' fraction of total XIST expression contributed by the inactive X chromosome.
#' @seealso getGenicDP
#' @importFrom ggplot2 geom_bar
#'
#' @export
getCellFrac <- function(dt, readCountThresh = 10, fix_phasing = TRUE, overall_mean = F, xist_only = F, merge_xist = TRUE, xciGenes = NULL, plot = FALSE){
  xci <- readXCI(xciGenes)
  dt_xci <- dt[GENE %in% xci]
  #dt <- copy(dt[GENE %in% c(xci, "XIST"),])

  dt_xci <- dt_xci[, frac := AD_hap1/tot]
  if(fix_phasing){
    # If an XCI gene is mostly expressed by the Xi, the phasing is likely wrong
    dt_xci <- dt_xci[frac > .5, frac := 1-frac] #Set fraction that are above .5 to 1-frac
  }
  dt_xci <- dt_xci[tot < readCountThresh, frac := NA]

  if(overall_mean){
    dt_xci[, frac_mean := sum(pmin(AD_hap1, AD_hap2))/(sum(AD_hap1) + sum(AD_hap2)), by = "sample"]
  } else{
    dt_xci <- dt_xci[, frac_mean := mean(frac, na.rm = TRUE), by = "sample"]
  }
  dt_xci <- dt_xci[, frac_median := median(frac, na.rm = TRUE), by = "sample"]

  dt_frac <- unique(dt_xci[, list(sample, frac_median, frac_mean)])
  if(xist_only){
    dt_xist <- dt[GENE == "XIST"]
  } else{
    dt_xist <- dt[grep("XIST", GENE),]
    if(merge_xist){
      dt_xist[, AD_hap1 := sum(AD_hap1), by = "sample"]
      dt_xist[, AD_hap2 := sum(AD_hap2), by = "sample"]
      dt_xist[, n_snps := sum(n_snps), by = "sample"]
      dt_xist[, tot := sum(tot), by = "sample"]
      dt_xist[, GENE := "XIST", by = "sample"]
      dt_xist <- unique(dt_xist[, list(sample, AD_hap1, tot)])
    }
  }
  dt_xist[, frac_XIST := AD_hap1/tot, by = "sample"]
  if(length(unique(dt_xist[, sample])) < length(unique(dt[,sample]))){
    message(paste0(length(unique(dt_xist[, sample])), "/", length(unique(dt[,sample])), " samples have available expression for XIST."))
    dt_frac <- dt_frac[sample %in% unique(dt_xist$sample)]
  }
  dt_frac <- merge(dt_frac, dt_xist[, list(sample, frac_XIST)], by = "sample")

  if(plot){
    mdt_frac <- melt(dt_frac, id.vars = "sample")
    p <- ggplot(mdt_frac, aes(sample, value, fill = variable)) +
      geom_bar(position = "dodge", stat = "identity") + ylab("Fraction") #stat = identity to plot y values instead of counts
    print(p)
  }
  return(dt_frac)
}

#' Estimate inactivated X chromosome expression
#'
#' Calculate inactivated X chromosome expression.
#'
#' @details
#' Select the females. calculate the proportion of expression coming from the
#' inactivated X chromosome for each gene.
#'
#' When the allelic coverage is very big for some genes, it gives more power,
#' which allow to spot smaller significant differences. When only major
#' differences are of interest, the \code{dp_max} parameter can be set to a
#' lower value.
#'
#' @param dt A \code{data.table}. The data. It must contain a gender column.
#' @param dt_frac A \code{data.table}. The result of the \code{getCellFrac}
#'  function.
#' @param dp_max A \code{numeric}. If any gene has a total coverage above this
#'  number, the coverage will be set to this number.
#' @param output A \code{logical}. If set to TRUE, output a table for each
#'  cutoff level.
#' @param cutoffs A \code{numeric} or NA.
#'
#' @export
getXIexpr <- function(dt, dt_frac, dp_max = NULL, fix_phasing = TRUE, output = FALSE, cutoff = NA){
  dt[, frac := AD_hap1/tot]
  dt <- dt[tolower(gender) == "female"]
  if(!is.null(dp_max)){
    dt[tot > dp_max, AD_hap1 := frac*dp_max]
    dt[tot > dp_max, tot := dp_max]
  }
  full_dt <- merge(dt, dt_frac, by = "sample")
  if(fix_phasing){
    full_dt[frac > .5, frac := 1-frac]
    full_dt[frac_mean > .5, frac_mean := 1-frac_mean]
  }

  dt <- copy(full_dt)
  dt[, `:=`(c("total_cutoff"), tot)]

  if(!is.na(cutoff)){
    dt[, cutoff := cutoff]
    dt <- dt[tot >= cutoff]
  }

    # frac_mean
    dt[, tau := (frac-frac_mean)/(1-frac_mean-frac)]
    dt[, var_tau := (2*frac_mean-1)^2/(frac+frac_mean-1)^4*(1-frac_mean)*frac_mean/total_cutoff]
    dt[, sd_tau := sqrt(var_tau)]
    dt[, p_value := pnorm(tau/sd_tau, lower.tail = FALSE)]

    #dt[, frac := pmin(AD_hap1, AD_hap2)/(AD_hap1 + AD_hap2)]
    dt[, tau_var := (frac_mean - frac)/(2*frac_mean - 1)]
    dt[, sd_var := sqrt((frac_mean*(1-frac_mean))/tot + (frac*(1-frac))/tot)]
    dt[, p_value_var := pnorm(tau_var/sd_var, lower.tail = FALSE)][]

    # frac_median
    dt[, tau_med := (frac-frac_median)/(1-frac_median-frac)]
    dt[, var_tau_med := (2*frac_median-1)^2/(frac+frac_median-1)^4*(1-frac_median)*frac_median/total_cutoff]
    dt[, sd_tau_med := sqrt(var_tau_med)]
    dt[, p_value_med := pnorm(tau_med/sd_tau_med, lower.tail = FALSE)]

    # frac_XIST
    dt[, tau_xist := (frac-frac_XIST)/(1-frac_XIST-frac)]
    dt[, var_tau_xist := (2*frac_XIST-1)^2/(frac+frac_XIST-1)^4*(1-frac_XIST)*frac_XIST/total_cutoff]
    dt[, sd_tau_xist := sqrt(var_tau_xist)]
    dt[, p_value_xist := pnorm(tau_xist/sd_tau_xist, lower.tail = FALSE)]

    dt[, total_cutoff := NULL][]
    if(output)
      write.table(dt[!is.na(normExp) & !is.nan(frac_XIST)], file=paste0("XIexpr_", cutoff)) # Write the output for the cutoff
  return(dt)
}

# uses the output of getCellFrac2
#' @export
getXIexpr3 <- function(cellfrac_dt){
  xi <- copy(cellfrac_dt)

  # Two sample
  xi[, tau2 := (fg - f)]#/(2*f -1)]# close to 1 when samples are skewed
  xi[, t2 := tau2/sqrt(var_f + var_fg)]
  xi[, p_value2 := pnorm(t2, lower.tail = F)]

  xi[, tau_xist2 := (fg - frac_XIST)]#/(2*f -1)]# close to 1 when samples are skewed
  xi[, t_xist2 := tau_xist2/sqrt(var_frac_XIST + var_fg)]
  xi[, p_value_xist2 := pnorm(t_xist2, lower.tail = F)]

  # One sample
  xi[, tau1 := (fg - f)]#/(2*f -1)]# close to 1 when samples are skewed
  xi[, t1 := tau1/sqrt(var_fg)]
  xi[, p_value1 := pnorm(t1, lower.tail = F)]

  xi[, tau_xist1 := (fg - frac_XIST)]#/(2*f -1)]# close to 1 when samples are skewed
  xi[, t_xist1 := tau_xist1/sqrt(var_fg)]
  xi[, p_value_xist1 := pnorm(t_xist1, lower.tail = F)]

  #mean
  xi[, tau := (fg-f)/(1-f-fg)]
  xi[, var_tau := (2*f-1)^2/(fg+f-1)^4*(1-f)*f/tot]
  xi[, sd_tau := sqrt(var_tau)]
  xi[, p_value := pnorm(tau/sd_tau, lower.tail = FALSE)]
  # frac_median
  xi[, tau_med := (fg-frac_median)/(1-frac_median-fg)]
  xi[, var_tau_med := (2*frac_median-1)^2/(fg+frac_median-1)^4*(1-frac_median)*frac_median/tot]
  xi[, sd_tau_med := sqrt(var_tau_med)]
  xi[, p_value_med := pnorm(tau_med/sd_tau_med, lower.tail = FALSE)]
  # frac_XIST
  xi[, tau_xist := (fg-frac_XIST)/(1-frac_XIST-fg)]
  xi[, var_tau_xist := (2*frac_XIST-1)^2/(fg+frac_XIST-1)^4*(1-frac_XIST)*frac_XIST/tot]
  xi[, sd_tau_xist := sqrt(var_tau_xist)]
  xi[, p_value_xist := pnorm(tau_xist/sd_tau_xist, lower.tail = FALSE)]

  return(xi[])
}
