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
#' @return A \code{character} vector of gene names.
#'
#' @examples
#' xcig <- readXCI()
#' xcig <- readXCI("cotton")
#'
#' @seealso \code{betaBinomXI}
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
  skewed_samples <- as.character(data[, median(pmin(AD_hap1, AD_hap2)/(AD_hap1+AD_hap2), na.rm = TRUE),
                                      by = sample][order(V1)][seq_len(n), sample])
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


#' Sample estimates
#'
#' Return sample specific information from XCIR results
#'
#' @param bb_table A \code{data.table}. The table returned by \code{betaBinomXI}.
#'
#' @return A \code{data.table} with one entry per sample and information
#'  regarding skewing and model fitting.
#'
#' @example inst/examples/betaBinomXI.R
#'
#' @export
sample_clean <- function(bb_table){
  clean_cols <- c("sample", "model", "f", "a_est", "b_est")
  ret <- unique(bb_table[, clean_cols, with = FALSE])
  return(ret)
}
xcir_clean <- function(bb_table){
  clean_cols <- c("sample", "GENE", "AD_hap1", "AD_hap2", "f", "p_value", "pbb")
  ret <- unique(bb_table[, clean_cols, with = FALSE])
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
#' @example inst/examples/betaBinomXI.R
#'
#' @export
getXCIstate <- function(xciObj){
  if(!"status" %in% names(xciObj))
    xciObj[, status := ifelse(p_value < 0.05, "E", "S")]
  out <- setkey(xciObj, GENE, status)[, .N, by = c("GENE", "status")][CJ(GENE, status, unique = TRUE), allow.cartesian = TRUE][is.na(N), N := 0L]
  out[, Ntot := sum(N), by = "GENE"]
  outE <- out[status == "E"]#
  outE[, pe := N/Ntot]
  ret <- outE[, list(GENE, Ntot, pe)]
  ret[, XCIstate := ifelse(pe <= .25, "S", "VE")]
  ret[, XCIstate := ifelse(pe >= .75, "E", XCIstate)]
  return(ret)
}

.betaAB <- function(m, theta, mu, sigma2){
  if(!is.null(m) & !is.null(theta)){
    alpha <- m*(theta-2)+1
    beta <- (1-m)*(theta-2)+1
  } else if(!is.null(mu) & !is.null(sigma2)){
    v <- (mu * (1-mu))/sigma2 - 1
    alpha <- mu*v
    beta <- (1-mu)*v
  } else{
    stop("At least one pair of m/theta or mu/sigma2 must be specified")
  }
  return(c(alpha, beta))
}
.betaMT <- function(alpha, beta, mu, sigma2){
  if(!is.null(alpha) & !is.null(beta)){
    theta <- alpha + beta
    m <- (alpha-1)/(theta-2)
  } else if(!is.null(mu) & !is.null(sigma2)){
  } else{
    stop("At least one pair of alpha/beta or mu/sigma2 must be specified")
  }
  return(c(m, theta))
}
.betaMV <- function(alpha, beta, m, theta){
  if(!is.null(alpha) & !is.null(beta)){
    mu <- alpha/(alpha + beta)
    sigma2 <- (alpha*beta)/((alpha+beta)^2 * (alpha+beta+1))
  } else if(!is.null(m) & !is.null(theta)){
    mu <- m*(theta-2)+1/theta
    # TODO: Simplify this
    sigma2 <- ((m*(theta-2)+1) * ((1-m)*(theta-2)+1)) / # a*b
      (( m*(theta-2)+1 + (1-m)*(theta-2)+1)^2 *  # (a+b)^2
         ( m*(theta-2)+1 + (1-m)*(theta-2)+1 + 1))  # (a+b+1)
  } else{
    stop("At least one pair of alpha/beta or m/theta must be specified")
  }
  return(c(mu, sigma2))
}


#' Converting beta distribution parameters
#'
#' Convert parameter values between different beta distribution parametrization
#'
#' @param alpha First shape parameter
#' @param beta Second shape parameter
#' @param m Mode
#' @param theta Concentration
#' @param mu Mean
#' @param sigma2 Variance
#'
#' @return A \code{list} with the two new parameters.
#'
#' @examples
#'
#' betaParam(alpha = 5, beta = 5)
#' betaParam(m = 0.5, theta = 10)
#' betaParam(mu = 0.5, sigma2 = 0.02272727)
#'
#' @rdname BetaConversion
#' @export
betaParam <- function(alpha = NULL, beta = NULL, m = NULL, theta = NULL, mu = NULL, sigma2 = NULL){
  if(is.null(alpha)){
    ab <- .betaAB(m, theta, mu, sigma2)
    alpha <- ab[1]
    beta <- ab[2]
  }
  if(is.null(m)){
    mt <- .betaMT(alpha, beta, mu, sigma2)
    m <- mt[1]
    theta <- mt[2]
  }
  if(is.null(mu)){
    mv <- .betaMV(alpha, beta, m, theta)
    mu <- mv[1]
    sigma2 <- mv[2]
  }
  return(c(alpha, beta, m, theta, mu, sigma2))
}
