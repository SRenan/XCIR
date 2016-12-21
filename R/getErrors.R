#' Get proportion of errors
#'
#' Get the proportion of type 1 and type 2 error for each sample for selected
#' cell fraction type
#'
#' @param xi A \code{data.table}. The result of \code{getXIexpr}. The predicted
#'  status for each gene and sample.
#' @param frac_type A \code{character}, can be NULL, "mean", "median" or "XIST"
#'  NULL defaults to "mean".
#' @param alpha A \code{numeric}. The significance threshold.
#' @param cutoff A \code{numeric}. Escaped genes are genes for which the
#'  inactivated copy has more reads than this cutoff.
#' @param min_read A \code{numeric}. The minimum number of reads for a gene to
#'  be considered informative.
#' @param plot A \code{logical}. Create a barplot of the Type 1 error and power.
#'
#' @importFrom ggplot2 geom_bar
#' @export
getError <- function(xi, frac_type = NULL, alpha = .05, cutoff = 1, min_read = 0, plot = TRUE){
  # Take the two most skewed samples as the truth
  samples <- as.character(unique(xi$sample))
  truth <- unique(xi[, list(sample, frac_median)])[order(frac_median), sample][1:2]
  others <- samples[!samples %in% truth]

  xi <- copy(xi)
  xi <- xi[tot > min_read]

  if(cutoff < 1){
    escaped_genes <- xi[sample %in% truth, min(AD_hap1, AD_hap2)/tot, by = "GENE,sample"][V1 > cutoff, .N, by = "GENE"][N ==2, GENE]
    inactivated_genes <- unique(xi[!GENE %in% escaped_genes, GENE])
  } else{
    escaped_genes <- xi[sample %in% truth & AD_hap1 > cutoff & AD_hap2 > cutoff, .N, by = GENE][N==2, GENE]
    inactivated_genes <- xi[(sample %in% truth) & ((AD_hap1 <= cutoff & AD_hap2 > cutoff) | (AD_hap1 > cutoff & AD_hap2 <= cutoff))
                            , .N, by = GENE][N==2, GENE]
  }

  print(paste(length(escaped_genes), "escaped genes accross the 2 most skewed samples"))
  print(paste(length(inactivated_genes), "inactivated genes accross the 2 most skewed samples"))

  if(!is.null(frac_type)){
    if(tolower(frac_type) == "median"){
      err_table <- xi[, list(GENE, sample, tot, p_value_med)]
    } else if(tolower(frac_type) == "xist"){
      err_table <- xi[, list(GENE, sample, tot, p_value_xist)]
    } else if(tolower(frac_type) == "mean"){
      err_table <- xi[, list(GENE, sample, tot, p_value)]
    } else{
      stop("frac_type should be NULL or one of median, mean or XIST")
    }
    setnames(err_table, c("GENE", "sample", "tot", "p_value"))
  } else{
    frac_type <- "mean"
    err_table <- xi[, list(GENE, AD_hap1, AD_hap2, tot, sample, p_value)]
  }

  nt1 <- err_table[GENE %in% inactivated_genes & sample %in% others & p_value < alpha][, .N, by = sample] #number of genes with type1 error
  if(!all(others %in% nt1[["sample"]])){
    ot <- others[!others %in% nt1$sample]
    nt1 <- rbind(nt1, data.table(sample = ot, N = rep(0, length(ot))))[order(sample)]
  }
  n_info <- err_table[GENE %in% inactivated_genes & sample %in% others & tot > min_read][, .N, by = sample]
  err1 <- merge(nt1, n_info, by = "sample")
  err1[, value := N.x/N.y]
  err1[, type := "typeI"]
  nt2 <- err_table[GENE %in% escaped_genes & sample %in% others & p_value > alpha][, .N, by = sample] #number of genes with type1 error
  if(!all(others %in% nt2[["sample"]])){
    ot <- others[!others %in% nt2$sample]
    nt2 <- rbind(nt2, data.table(sample = ot, N = rep(0, length(ot))))[order(sample)]
  }
  n_info <- err_table[GENE %in% escaped_genes & sample %in% others & tot > min_read][, .N, by = sample]
  err2 <- merge(nt2, n_info, by = "sample")
  err2[, N.x := N.y - N.x] #T2 -> power
  err2[, value := N.x/N.y]
  err2[, type := "power"]

  err <- rbind(err1, err2)
  err[, fraction := frac_type]

  if(plot){
    p <- ggplot(err, aes(x = sample, y = value, fill = type)) + geom_bar(position = "dodge", stat = "identity")
    p <- p + geom_text(aes(x = sample, y = value, label = paste0(N.x, "/", N.y)), position = position_dodge(0.9), vjust = -.25)
    p <- p + geom_text(aes(x = 1.5, y = 1, label = paste0(length(escaped_genes), " escaped genes"))) +
     geom_text(aes(x = 1.5, y = .95, label = paste0(length(inactivated_genes), " inactivated genes")))
    print(p)
  }

  setnames(err, c("N.x", "N.y"), c("N", "N_info"))
  err[, value := round(value, 3)][]
  return(err)
}
