#' Get sum of the coverage of all SNPs for each gene
#'
#' @importFrom ggplot2 ggplot geom_histogram xlab ylab
#' @export
getGenicDP <- function(dt_anno, plot = FALSE, gender_file = NULL){
  if(!is.null(gender_file)){
    sex <- fread(gender_file)
    setnames(sex, c("sample", "gender"))
    sex[, gender := tolower(gender)]
    dt_anno <- merge(dt_anno, sex, by = "sample")
  }
  else if(!"gender" %in% names(dt_anno)){
    stop("If there is no gender column in the input, a gender_file mapping
         sample to gender should be provided")
  }
  dt_anno[, n_snps := .N, by = c("CHROM", "sample", "GENE")]
  tot <- dt_anno[, lapply(.SD, sum, na.rm = TRUE), by = c("CHROM", "gender", "GENE", "sample", "n_snps"), .SDcols = c("AD_hap1", "AD_hap2")] #nrow = nGenes * nSubs
  tot[, tot := AD_hap1 + AD_hap2]
  if(plot){
    ggplot(tot[, mean(tot), by = "GENE"], aes(V1)) + geom_histogram() + xlab(label = "coverage") + ylab("gene count")
  }
  return(tot)
}

#' Get fraction of cell expressing the inactivated allele
#'
#' @param dt A \code{data.table}. The output of \code{getGenicDP}.
#' @param xciGenes A \code{character}. The name of the genes to analyze.
#' @param fix_phasing A \code{logical}. If set to TRUE, the fraction above .5
#'  are set to 1-frac.
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
#'
#' @export
getCellFrac <- function(dt, fix_phasing = FALSE, xciGenes = NULL){
  if(is.null(xciGenes))
    xciGenes <- system.file("extdata", "xciGene.txt", package = "XCIR")
  xci <- readLines(xciGenes)
  dt <- dt[GENE %in% c(xci, "XIST"),]
  dt <- dt[, frac := AD_hap1/tot]
  if(fix_phasing){
    # If an XCI gene is mostly expressed by the Xi, the phasing is likely wrong
    dt <- dt[frac > .5, frac := 1-frac] #Set fraction that are above .5 to 1-frac
  }
  dt <- dt[tot < 10, tot := NA]
  dt <- dt[, frac_mean := mean(frac, na.rm = TRUE), by = "sample"]
  dt <- dt[, frac_median := median(frac, na.rm = TRUE), by = "sample"]

  dt_frac <- unique(dt[GENE == "XIST", list(sample, frac, frac_median, frac_mean)])
  setnames(dt_frac, "frac", "frac_XIST")
  return(dt_frac) #This is fracCell.allXci.txt
}

#Input:
# The output of estimated genicDP
# A list of genes that are known to NOT espcape XCI.
# The output of getCellFrac
#Output:
# For each cutoff: The pI, pval, sd, quantile

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
#'
#' @export
getXIexpr <- function(dt, dt_frac, dp_max = 500, fix_phasing = FALSE, output = FALSE){
  #if(is.null(xciGenes))
  #  xciGenes <- system.file("extdata", "xciGene.txt", package = "XCIR")
  #xci <- readLines(xciGenes)
  #dt <- dt[GENE %in% xci]
  dt[, frac := AD_hap1/tot]
  dt <- dt[tolower(gender) == "female"]
  dt[tot > dp_max, AD_hap1 := frac*dp_max]
  dt[tot > dp_max, tot := dp_max]
  full_dt <- merge(dt, dt_frac, by = "sample")
  if(fix_phasing){
    full_dt[frac > .5, frac := 1-frac]
    full_dt[frac_mean > .5, frac_mean := 1-frac_mean]
  }
  cutoffs <- seq(50, 10, -10)
  l <- vector("list", length(cutoffs))
  i <- 1
  #At each cutoff, consider only the genes that have a total read count higher than the cutoff
  for(dp_cutoff in cutoffs){
    # Anything that depends on the coverage is affected by the cutoff
    dt <- copy(full_dt)
    dt[, cutoff := dp_cutoff]
    dt[tot < dp_cutoff, `:=`(c("tot", "AD_hap1"), list(NA, NA))]

    dt[, tau := (frac-frac_mean)/(1-frac_mean-frac)]
    dt[, var_tau := (2*frac_mean-1)^2/(frac+frac_mean-1)^4*(1-frac_mean)*frac_mean/tot]
    dt[, sd_tau := sqrt(var_tau)]
    dt[, p_value := pnorm(tau/sd_tau, lower.tail = FALSE)]

    dt[, pI := (AD_hap1/tot + frac_mean - 1)/(2*frac_mean-1), by = c("GENE", "sample")]
    dt[, sd_pI := sqrt((AD_hap1/tot)*(1-AD_hap1/tot)/(tot))/abs(2*frac_mean-1), by = c("GENE", "sample")]
    dt[, pval_pI := pnorm(pI/sd_pI, mean = 0, lower.tail = FALSE), by = c("GENE", "sample")]
    dt[, t := pI/sd_pI]
    dt[, normExp := qnorm(pval_pI, lower.tail = FALSE)]
    if(output)
      write.table(dt[!is.na(normExp) & !is.nan(frac_XIST)], file=paste0("XIexpr_", dp_cutoff)) # Write the output for the cutoff
    l[[i]] <- dt
    i <- i+1
  }
  return(rbindlist(l))
}

