# For now assume no NA (=no cutoff in getXT). Later, change getXIexpr
predMap <- function(xi, cutoff = 1, alpha = .05, nGenes = NULL){
  xi_table <- copy(xi)
  # Organize samples
  xi_table[, sample := gsub("^.*_", "", gsub("^.*-", "", sample))]
  samples <- as.character(unique(xi_table$sample))
  truth <- as.character(unique(xi_table[, list(sample, frac_median)])[order(frac_median), sample][1:2])
  others <- samples[!samples %in% truth]

  # status for the skewed samples
  if(cutoff < 1){
    escaped_genes <- xi_table[sample %in% truth, min(AD_hap1, AD_hap2)/tot, by = "GENE,sample"][V1 > .1, .N, by = "GENE"][N ==2, GENE]
    inactivated_genes <- unique(xi_table[!GENE %in% escaped_genes, GENE])
  } else{
    escaped_genes     <- xi_table[sample %in% truth & AD_hap1 > cutoff & AD_hap2 > cutoff, .N, by = GENE][N == 2, GENE]
    inactivated_genes <- xi_table[sample %in% truth & ((AD_hap1 <= cutoff & AD_hap2 > cutoff) | (AD_hap1 > cutoff & AD_hap2 <= cutoff))
                                  , .N, by = GENE][N == 2, GENE]
  }
  xi_table <- xi_table[GENE %in% c(escaped_genes, inactivated_genes)]
  xi_table[ GENE %in% escaped_genes & sample %in% truth, `:=`(c("status_mean", "status_median", "status_xist"), list("escaped", "escaped", "escaped"))]
  xi_table[ GENE %in% inactivated_genes & sample %in% truth, `:=`(c("status_mean", "status_median", "status_xist"), list("inactivated", "inactivated", "inactivated"))]

  # Find predicted status of each gene for the unknowns
  xi_table[sample %in% others, status_mean := ifelse(p_value < alpha, "escaped", "inactivated")]
  xi_table[sample %in% others, status_median := ifelse(p_value_med < alpha, "escaped", "inactivated")]
  xi_table[sample %in% others, status_xist := ifelse(p_value_xist < alpha, "escaped", "inactivated")]

  # Check whether predictions are correct
  idt1_mean <- xi_table[, status_mean != status_mean[sample == truth[1]] & status_mean == "escaped"]
  idt1_med  <- xi_table[, status_median != status_median[sample == truth[1]] & status_median == "escaped"]
  idt1_xist <- xi_table[, status_xist != status_xist[sample == truth[1]] & status_xist == "escaped"]
  idt2_mean <- xi_table[, status_mean != status_mean[sample == truth[1]] & status_mean == "inactivated"]
  idt2_med  <- xi_table[, status_median != status_median[sample == truth[1]] & status_median == "inactivated"]
  idt2_xist <- xi_table[, status_xist != status_xist[sample == truth[1]] & status_xist == "inactivated"]

  xi_table[sample %in% others, `:=`(c("status_mean", "status_median", "status_xist"), list("correct", "correct", "correct"))]
  xi_table[idt1_mean, status_mean := "typeI"]
  xi_table[idt1_med, status_median := "typeI"]
  xi_table[idt1_xist, status_xist := "typeI"]
  xi_table[idt2_mean, status_mean := "typeII"]
  xi_table[idt2_med, status_median := "typeII"]
  xi_table[idt2_xist, status_xist := "typeII"]

  # Plot and output
  palette <- RColorBrewer::brewer.pal(n = 6, "Set1") [c(3, 2, 6, 1, 5)]
  mxi <- melt(xi_table, measure.vars = c("status_mean", "status_median", "status_xist"), value.name ="status", variable.name = "fraction")
  if(!is.null(nGenes)){
    genes <- c(escaped_genes, inactivated_genes)
    mxi <- mxi[GENE %in% genes[1:nGenes]]
  }
  HMtheme <- theme(panel.grid = element_blank())
  if(length(unique(mxi$GENE)) > 100){
    HMtheme <- HMtheme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  p <- ggplot(mxi, aes(sample, GENE, colour = "black")) + geom_raster(aes(fill = status)) + scale_fill_manual(values = palette) + facet_wrap(~fraction) + HMtheme
  print(p)
  return(mxi)
}
