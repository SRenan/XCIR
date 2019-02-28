#' Read a list of known inactivated genes
#'
#' Read a list of gene symbols of known inactivated genes
#' to be used as training set in \code{betaBinomXI}.
#'
#' @param xciGenes A \code{character} or code{NULL}. By defaults, return a
#' vector of 177 genes. Other available choices include "cotton" and "intersect".
#' If a file path is given, the genes will be read from
#' the file.
#'
#' @details
#' The default gene list is from Carrel et al. Nature 2005. doi:10.1038/nature03479.
#' Cotton, is a list of 294 inactivated genes compiled by Cotton et al. Genome
#' Biology (2013). doi:10.1186/gb-2013-14-11-r122. Intersect is a more
#' stringent list that takes the intersection of both Carrel and Cotton's lists.
#'
#' @seealso \code{betaBinomXI}
#' @return A \code{character} vector of gene names.
#'
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
.getSkewedSamples <- function(data, n = 2){
  nsamples <- length(unique(data$sample))
  if(n < 1 | n > nsamples)
    stop(paste("n should be a number between 1 and", nsamples))
  skewed_samples <- as.character(data[, median(pmin(AD_hap1, AD_hap2)/(AD_hap1+AD_hap2), na.rm = T),
                                      by = sample][order(V1)][1:n, sample])
  return(skewed_samples)
}


# Find the allelic imbalance for each sample
# Based on formula in Cotton et al. Genome Biology 2013
.getAI <- function(calls, labels = FALSE){
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
.getAI <- function(calls){
  ai_dt <- unique(calls[, list(sample, GENE, AD_hap1, AD_hap2, tot, f, tau, p_value)])
  ai_dt[, xiexpr := pmin(AD_hap1, AD_hap2)/pmax(AD_hap1, AD_hap2)] #Assuming xaexpr is always 100%
  ai_dt[, num := (1-f) + xiexpr * f]
  ai_dt[, denom := ((1-f) * (xiexpr + 1)) + (f * (xiexpr + 1))]
  ai_dt[, ai := abs(num/denom - 0.5)]
}
# Samples with ai < cutoffs are subject to XCI
plot_escape_fraction <- function(ai, cutoff = .1){
  fai <- ai[, Nsamples := .N, by = "GENE"]
  fai <- fai[ai < cutoff, Ninac := .N, ]
}

#' Genes of interest
#'
#' These genes are particularly relevant to the biology of XCI.
#'
#' @return A \code{character} vector of gene symbols.
#' @export
#'
GOI <- function(){
  # Set of X/XCI relavant genes and aliases
  goi <- character()
  # DXZ4 is at the hinge that separates X into two superdomains
  goi <- c(goi, "DXZ4", "DANT2")
  # Genes that were found to be tumor suppressors in Male biased cancers (hypothesized
  # to be escaping, thus requiring 2 mutations in female to disrupt suppressing activity)
  goi <- c(goi, "ATRX", "KDM5C", "KDM6A", "DDX3X", "MAGEC3", "CNKSR2") # ASHG17 A.A. Lane, ATRX has never been found to escape so far. All others have.

  # Genes involved in XCI
  goi <- c(goi, "XIST", "TSIX", "TSIX|XIST", "FIRRE", "CTCF") # FIRRE helps maintain methylation on the Xi

  #SLE associated genes
  goi <- c(goi, "CXorf21", "IRAK1", "MECP2", "PRPS2")

  return(goi)
}

#' @export
sample_clean <- function(bb_table){
  clean_cols <- c("sample", "model", "f", "a_est", "b_est")
  ret <- unique(bb_table[, clean_cols, with = F])
  return(ret)
}
xcir_clean <- function(bb_table){
  clean_cols <- c("sample", "GENE", "AD_hap1", "AD_hap2", "f", "p_value", "pbb")
  ret <- unique(bb_table[, clean_cols, with = F])
  return(ret)
}

#' Classify X-genes
#'
#' Classify X-linked genes between Escape (E), Variable Escape (VE) and Silenced (S)
#'
#' @param xciObj A \code{data.table}. The table returned by \code{betaBinomXI}
#'
#' @return A \code{data.table} with genes and their XCI-state.
#'
#' @export
getXCIstate <- function(xciObj){
  if(!"status" %in% names(xciObj))
    xciObj[, status := ifelse(p_value < 0.05, "E", "S")]
  out <- setkey(xciObj, GENE, status)[, .N, by = c("GENE", "status")][CJ(GENE, status, unique = T), allow.cartesian = T][is.na(N), N := 0L]
  out[, Ntot := sum(N), by = "GENE"]
  outE <- out[status == "E"]#
  outE[, pe := N/Ntot]
  ret <- outE[, list(GENE, Ntot, pe)]
  ret[, XCIstate := ifelse(pe <= .25, "S", "VE")]
  ret[, XCIstate := ifelse(pe >= .75, "E", XCIstate)]
  return(ret)
}


#' Converting beta distribution parameters
#'
#' These functions convert parameter values between different beta distribution parametrization
#'
#' @param alpha First shape parameter
#' @param beta Second shape parameter
#' @param m Mode
#' @param theta Concentration
#' @param mean Mean
#' @param var Variance
#'
#' @return A \code{list} with the two new parameters.
#'
#' @examples
#'
#' betaMV(1, 1)
#' betaMV(1, 10)
#' betaMV(10, 10)
#' # In a data.table (such as XCIR output)
#' dt <- data.table(it = LETTERS[1:3], a_est = c(1, 1, 10), b_est = c(1, 10, 10))
#' dt[, `:=`(c("mu", "s2"), betaMV(a_est, b_est))]
#'
#'
#' @rdname BetaConversion
#' @export
betaAB <- function(mean, var){
  alpha <- (mu^2*(1-mu))/var - 1
  beta <- (mu*(1-mu)^2)/var - 1
  return(list(alpha, beta))
}
#' @rdname BetaConversion
#' @export
betaMC <- function(alpha, beta){
  theta <- alpha+beta
  m <- (alpha-1)/(theta-2)
  return(list(m, theta))
}
#' @rdname BetaConversion
#' @export
betaMV <- function(alpha, beta){
  mu <- alpha/(alpha + beta)
  var <- (alpha*beta)/((alpha+beta)^2 * (alpha+beta+1))
  return(list(mu, var))
}
