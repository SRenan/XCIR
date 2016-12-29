getError <- function(xi, frac_type = NULL, alpha = .05, escape_cutoff = 1, inac_cutoff = 2, plot = TRUE, xciGenes = NULL){
  # Take the two most skewed samples as the truth
  samples <- as.character(unique(xi$sample))
  truth <- unique(xi[, list(sample, frac_median)])[order(frac_median), sample][1:2]
  others <- samples[!samples %in% truth]


  xi <- copy(xi)
  xci_genes <- XCIR:::readXCI()

  if(cutoff < 1){
    escaped_genes <- xi[sample %in% truth, min(AD_hap1, AD_hap2)/tot, by = "GENE,sample"][V1 > cutoff, .N, by = "GENE"][N ==2, GENE]
    inactivated_genes <- unique(xi[!GENE %in% escaped_genes, GENE])
  } else{
    escaped_genes <- xi[sample %in% truth & AD_hap1 > escape_cutoff & AD_hap2 > escape_cutoff, .N, by = GENE][N==2, GENE]
    dt_xci <- xi[GENE %in% xci_genes]
    notinactivated <- unique(c(dt_xci[grep("LC2$", sample)][AD_hap1 > inac_cutoff & AD_hap2 > inac_cutoff, GENE],
                               dt_xci[grep("LC3$", sample)][AD_hap1 > inac_cutoff & AD_hap2 > inac_cutoff, GENE]))
    inactivated_genes <- unique(dt_xci[!GENE %in% notinactivated, GENE])
    #inactivated_genes <- xi[(sample %in% truth) & ((AD_hap1 <= cutoff & AD_hap2 > cutoff) | (AD_hap1 > cutoff & AD_hap2 <= cutoff))
    #                        , .N, by = GENE][N==2, GENE]
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
