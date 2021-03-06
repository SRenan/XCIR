#' Plot error and power
#'
#' Plot power and type 1 error for each sample. Using the most skewed samples
#' as a reference. Used to asses the quality of the prediction.
#'
#' @param xi A \code{data.table}. Escape calls, as outputted by \code{betaBinomXI}.
#' @param alpha A \code{numeric}. The significance level.
#' @param escape_cutoff A \code{numeric}. A gene must have more reads than the
#' cutoff on both alleles in the two skewed samples to be considered as
#' escaping XCI.
#' @param inac_cutoff A \code{numeric}. A gene must have at most that many reads
#' on both alleles in the two skewed samples to be considered as subject to
#' inactivation.
#' @param plot A \code{logical}. If set to false, return Type1 error and power
#' table without plotting barplots.
#' @param xciGenes A \code{logical}. If set to TRUE, uses genes from the
#' training set as inactivated genes (still subject to filtering via
#' \code{inac_cutoff}). If set to FALSE, a hand curated gene list will be used
#' insteda.
#'
#' @note
#' This will only work if all samples belong to the same subject and two of the
#' samples are mono-allelicly expressed.
#'
#' @importFrom ggplot2 ggplot position_dodge geom_errorbar geom_text
#' @export
getError <- function(xi, alpha = .05, escape_cutoff = 1, inac_cutoff = 2, plot = TRUE, xciGenes = TRUE){
  # Take the two most skewed samples as the truth
  samples <- as.character(unique(xi$sample))
  truth <- .getSkewedSamples(xi)
  others <- samples[!samples %in% truth]

  xi <- copy(xi)

  escaped_genes <- xi[sample %in% truth & AD_hap1 > escape_cutoff & AD_hap2 > escape_cutoff, .N, by = GENE][N==2, GENE]
  if(xciGenes){
    xci_genes <- readXCI()
    xci_genes[! xci_genes %in% gNrm] # Not needed with the mixture-model
    dt_xci <- xi[GENE %in% xci_genes]
    notinactivated <- unique(c(dt_xci[grep(truth[1], sample)][AD_hap1 > inac_cutoff & AD_hap2 > inac_cutoff, GENE],
                               dt_xci[grep(truth[2], sample)][AD_hap1 > inac_cutoff & AD_hap2 > inac_cutoff, GENE]))
    inactivated_genes <- unique(dt_xci[!GENE %in% notinactivated, GENE])
  } else{
    inactivated_genes <- jcss[statusJ == "S", GENE]
    escaped_genes <- jcss[statusJ == "E", GENE]
  }

  message(paste(length(escaped_genes), "escaped genes accross the 2 most skewed samples"))
  message(paste(length(inactivated_genes), "inactivated genes accross the 2 most skewed samples"))

  err_table <- xi[, list(GENE, AD_hap1, AD_hap2, tot, sample, p_value)]

  nt1 <- err_table[GENE %in% inactivated_genes & sample %in% others & p_value < alpha][, .N, by = sample] #number of genes with type1 error
  if(!all(others %in% nt1[["sample"]])){
    ot <- others[!others %in% nt1$sample]
    nt1 <- rbind(nt1, data.table(sample = ot, N = rep(0, length(ot))))[order(sample)]
  }
  n_info <- err_table[GENE %in% inactivated_genes & sample %in% others][, .N, by = sample]
  err1 <- merge(nt1, n_info, by = "sample")
  err1[, value := N.x/N.y]
  err1[, type := "typeI"]
  nt2 <- err_table[GENE %in% escaped_genes & sample %in% others & p_value > alpha][, .N, by = sample] #number of genes with type1 error
  if(!all(others %in% nt2[["sample"]])){
    ot <- others[!others %in% nt2$sample]
    nt2 <- rbind(nt2, data.table(sample = ot, N = rep(0, length(ot))))[order(sample)]
  }
  n_info <- err_table[GENE %in% escaped_genes & sample %in% others][, .N, by = sample]
  err2 <- merge(nt2, n_info, by = "sample")
  err2[, N.x := N.y - N.x] #T2 -> power
  err2[, value := N.x/N.y]
  err2[, type := "power"]

  err <- rbind(err1, err2)

  setnames(err, c("N.x", "N.y"), c("N", "N_info"))
  err[type == "typeI", ymin := value -1.96*sqrt((N/N_info)*(1-N/N_info)/N_info)]
  err[type == "typeI", ymax := value +1.96*sqrt((N/N_info)*(1-N/N_info)/N_info)]
  if(plot){
    .plot_errors(err, escaped_genes, inactivated_genes)
  }

  err[, value := round(value, 3)][]
  return(err)
}

.plot_errors <- function(err, escaped_genes, inactivated_genes){
  p <- ggplot(err, aes(x = sample, y = value, fill = type)) + geom_bar(position = "dodge", stat = "identity")
  p <- p + geom_text(aes(x = sample, y = value, label = paste0(N, "/", N_info)), position = position_dodge(0.9), vjust = -.25)
  p <- p + geom_text(aes(x = 1.5, y = 1, label = paste0(length(escaped_genes), " escaped genes"))) +
    geom_text(aes(x = 1.5, y = .95, label = paste0(length(inactivated_genes), " inactivated genes")))
  p <- p + geom_errorbar(aes(ymax = ymax, ymin = ymin), position = position_dodge(.9), width = .25, na.rm = TRUE)
  print(p)
  return(p)
}
