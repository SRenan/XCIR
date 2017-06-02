# Read a list of known inactivated genes
#' @export
readXCI <- function(xciGenes = NULL){
  if(!is.null(xciGenes)){
    if(length(xciGenes) > 1){
      # If the input is a character vector, assume it's gene symbols
      xci <- xciGenes
      return(xciGenes)
    } else if(xciGenes == "cotton"){
      xci <- system.file("extdata", "xciGenes_cotton.txt", package = "XCIR")
      xci <- readLines(xci)
    } else if(xciGenes == "intersect"){
      xci177 <- readLines(system.file("extdata", "xciGene.txt", package = "XCIR"))
      xcic <- readLines(system.file("extdata", "xciGenes_cotton.txt", package = "XCIR"))
      xci <- intersect(xci177, xcic)
    } else if(file.exists(xciGenes)){
      xci <- readLines(xciGenes)
    } else{
      stop("The file does not exist")
    }
  } else{
    xci <- system.file("extdata", "xciGene.txt", package = "XCIR")
    xci <- readLines(xci)
  }
}

# Find the most skewed samples within an XCIR data.table
getSkewedSamples <- function(data, n = 2){
  nsamples <- length(unique(data$sample))
  if(n < 1 | n > nsamples)
    stop(paste("n should be a number between 1 and", nsamples))
  skewed_samples <- as.character(data[, median(pmin(AD_hap1, AD_hap2)/(AD_hap1+AD_hap2)), by = sample][order(V1)][1:n, sample])
  return(skewed_samples)
}


# Remove genes that are not inactivated in the 2 most skewed samples
rmInac <- function(dt_xci, rm_inac){
  if(!is.numeric(rm_inac)){
    stop("rm_inac should be null or numeric")
  }
  xci_dt <- copy(dt_xci)
  skewed <- getSkewedSamples(xci_dt)
  notinactivated <- unique(c(xci_dt[grep(skewed[1], sample)][AD_hap1 > rm_inac & AD_hap2 > rm_inac, GENE],
                             xci_dt[grep(skewed[2], sample)][AD_hap1 > rm_inac & AD_hap2 > rm_inac, GENE]))
  xci_dt <- xci_dt[!GENE %in% notinactivated]
  print(paste("Using", length(unique(xci_dt[, GENE])), "inactivated genes to calculate cell fraction."))
  return(xci_dt)
}

# Find the allelic imbalance for each sample
# Based on formula in Cotton et al. Genome Biology 2013
getAI <- function(calls, labels = FALSE){
  ai_dt <- unique(calls[, list(sample, GENE, POS, AD_hap1, AD_hap2, tot, f)])
  ai_dt <- unique(ai_dt[, list(sum(pmin(AD_hap1, AD_hap2)), sum(tot), f), by = sample]) # Assume that lowest expressed allele is Xi
  ai_dt[, pxa := (V2-V1)/V2]
  ai_dt[, pxi := V1/V2]
  ai_dt[, num := (pxa * (1-f)) + (pxi * f)]
  ai_dt[, denom := f * (pxi + pxa) + (1-f) * (pxi + pxa)]
  ai_dt[, ai := abs(num/denom -0.5)]
  p <- ggplot(ai_dt[!is.na(f)], aes(ai, pxi)) + geom_point() + geom_smooth(method="loess") +
    ggtitle("%Xi expression vs. AI") + theme(plot.title = element_text(hjust = .5))
  if(labels)
    p <- p + geom_text(aes(label = sample), vjust = 1.2)
  print(p)
  return(ai_dt)
}
