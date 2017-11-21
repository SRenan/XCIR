#' Perform a permutation test
#'
#' Permutation test to assess the effect of the estimated Xi expression
#' (represented as tau, the fration of the total expression) on a phenotype.
#'
#' @param table A \code{data.table}. The output of \code{betaBinomXI}.
#' @param y A \code{character}. A column of \code{table}.
#' @param nperm A \code{numeric}. The number of permutations to perform.
#' @param tau_cutoff A \code{numeric}. Optional argument, if set to a numeric
#' value, all taus less than the cutoff will be set tot 0.
#' @param varweight A \code{logical}. If set to TRUE, taus are weighted by the
#' inverse of their variance.
#' @param verbose A \code{logical}. If set to TRUE, will post additional information
#' about the permutations.
#'
#' @return A \code{data.table} with the following information
#' \item{GENE}{The gene symbol}
#' \item{nsamples}{The number of samples used for the test}
#' \item{Tobs}{The test statistic}
#' \item{perm_pvalue}{The p-value}
#'
#' @export
perm_test <- function(table, y, nperm = 100, tau_cutoff = NULL, varweight = T, verbose = F){
  # Check outcome
  if(!y %in% colnames(table))
    stop(paste(y, "is not a valid column of the dataset"))
  data <- unique(table[, c("sample", "GENE", "tau", "var_tau", y), with = F])
  outcomes <- unique(data[, get(y)])
  if(length(outcomes) > 2)
    stop(paste0(y, " is not binary (length: ", length(outcomes), ")"))
  if(verbose)
    print(paste("Permutation test for", outcomes[1], "vs.", outcomes[2]))

  if(!is.null(tau_cutoff)){
    data[, tau := ifelse(tau < tau_cutoff, 0, tau)]
  }

  # Cannot do permutations if all samples belong to the same class
  gN2rm <- unique(data[, c("GENE", y), with = F])[, .N, by = GENE][N == 1, GENE]
  gN2keep <- unique(data[, c("GENE", y), with = F])[, .N, by = GENE][N == 2, GENE]
  message(paste(length(gN2rm), "genes cannot be processed as they only appear in one condition"))
  message(paste(length(gN2keep), "genes will be run through the permutation test"))
  data <- data[GENE %in% gN2keep]

  data[, nsamples := .N, by = "GENE"]
  if(varweight){
    data[, tau := tau/var_tau] #inverse variance weighted taus
  }

  genes <- unique(data$GENE)

  for(gene in genes){
    if(verbose){
      print(paste0("Processing ", gene, " (", which(gene == genes), "/", length(genes), ")"))
    }
    subdat <- data[GENE == gene]
    subout <- subdat[, get(y)]
    n <- subdat[, .N]
    pteststat <- numeric(nperm)
    for(perm in 1:nperm){
      permout <- sample(subout, n, replace = F)
      oc1idx <- which(permout == outcomes[1])
      oc1mean <- subdat[oc1idx, mean(tau)]
      oc2mean <- subdat[!oc1idx, mean(tau)]
      pteststat[perm] <- abs(oc1mean - oc2mean) #Test stat for this perm
    }
    obsteststat <- abs(subdat[GENE == gene &  get(y)  == outcomes[1], mean(tau)] -
                         subdat[GENE == gene & get(y) == outcomes[2], mean(tau)])
    pteststat <- sort(pteststat)
    pval <- 1 - sum(obsteststat > pteststat)/nperm
    data[GENE == gene, Tobs := obsteststat]
    data[GENE == gene, perm_pvalue := pval]
  }
  ret <- unique(data[, list(GENE, nsamples, Tobs, perm_pvalue)])
  return(ret)
}

perm2_test <- function(table, y, nperm = 1000, verbose = F){
  # Here we permute the dataset first
  # Check outcome
  if(!y %in% colnames(table))
    stop(paste(y, "is not a valid column of the dataset"))
  data <- unique(table[, c("sample", "GENE", "tau", "var_tau", y), with = F])
  outcomes <- unique(data[, get(y)])
  if(length(outcomes) > 2)
    stop(paste0(y, " is not binary (length: ", length(outcomes), ")"))
  if(verbose)
    print(paste("Permutation test for", outcomes[1], "vs.", outcomes[2]))

  # For each perm
  # - Get the number of samples in each categories
  outs <- as.character(data[, get(y)])
  n <- data[, .N]
  pteststat <- numeric(nperm)
  for(perm in 1:nperm){
    if(verbose){
      print(paste0("Permutation (", perm, "/", nperm, ")"))
    }
    permutated <- sample(outs, n, replace = F)
    data[, permout := permutated]
    data[, tau := ifelse(tau < 0, 0, tau)]
    data[, ivw_tau := tau/var_tau]

    gN2rm <- unique(data[, c("GENE", "permout"), with = F])[, .N, by = GENE][N == 1, GENE]
    gN2keep <- unique(data[, c("GENE", "permout"), with = F])[, .N, by = GENE][N == 2, GENE]
    subdat <- data[GENE %in% gN2keep]
    genes <- unique(subdat$GENE)
    gteststat <- numeric(length(genes))
    for(i in seq_along(genes)){
      #num <- (subdat[GENE == genes[i] & permout == outcomes[1], sum(ivw_tau)] - subdat[GENE == genes[i] & permout == outcomes[2], sum(ivw_tau)])^2
      #den <- subdat[GENE == genes[i] & permout == outcomes[1], sum(1/var_tau)]  + subdat[GENE == genes[i] & permout == outcomes[2], sum(1/var_tau)]
      #gteststat[i] <- num/den
      gteststat[i] <- abs(subdat[GENE == genes[i] & permout == outcomes[1], mean(ivw_tau, na.rm = T)] -
                          subdat[GENE == genes[i] & permout == outcomes[2], mean(ivw_tau, na.rm = T)])
    }
    pteststat[perm] <- max(gteststat, na.rm = T)
  }

  gN2rm <- unique(data[, c("GENE", y), with = F])[, .N, by = GENE][N == 1, GENE]
  gN2keep <- unique(data[, c("GENE", y), with = F])[, .N, by = GENE][N == 2, GENE]
  message(paste(length(gN2rm), "genes cannot be processed as they only appear in one condition"))
  message(paste(length(gN2keep), "genes will be run through the permutation test"))
  data <- data[GENE %in% gN2keep]

  data[, nsamples := .N, by = "GENE"]
  genes <- unique(data$GENE)

  for(gene in genes){
    # num <- (data[GENE == gene & get(y) == outcomes[1], sum(ivw_tau)] - data[GENE == gene & get(y) == outcomes[2], sum(ivw_tau)])^2
    # den <- data[GENE == gene & get(y) == outcomes[1], sum(1/var_tau)]  + data[GENE == gene & get(y) == outcomes[2], sum(1/var_tau)]
    # obsteststat <- num/den
    obsteststat <- abs(data[GENE == gene & get(y) == outcomes[1], mean(ivw_tau, na.rm = T)] -
                       data[GENE == gene & get(y) == outcomes[2], mean(ivw_tau, na.rm = T)])
    pval <- 1 - sum(obsteststat > pteststat)/nperm
    data[GENE == gene, Tobs := obsteststat]
    data[GENE == gene, perm_pvalue := pval]
  }
  ret <- unique(data[, list(GENE, Tobs, perm_pvalue)])
  return(ret)
}


# Use ESCAPE/XCI as a predictor
perm3_test <- function(table, y, nperm = 1000, alpha = .05, verbose = FALSE){
  # Check outcome
  if(!y %in% colnames(table))
    stop(paste(y, "is not a valid column of the dataset"))
  data <- unique(table[, c("sample", "GENE", "pbb", "p_value", y), with = F])
  outcomes <- unique(data[, get(y)])
  if(length(outcomes) > 2)
    stop(paste0(y, " is not binary (length: ", length(outcomes), ")"))
  if(verbose)
    print(paste("Permutation test for", outcomes[1], "vs.", outcomes[2]))

  # Check
  gN2rm <- unique(data[, c("GENE", y), with = F])[, .N, by = GENE][N == 1, GENE]
  gN2keep <- unique(data[, c("GENE", y), with = F])[, .N, by = GENE][N == 2, GENE]
  message(paste(length(gN2rm), "genes cannot be processed as they only appear in one condition"))
  message(paste(length(gN2keep), "genes will be run through the permutation test"))
  data <- data[GENE %in% gN2keep]

  genes <- unique(data$GENE)
  data[, escape := ifelse(p_value < alpha, 1, 0)]
  data[, nsamples := .N, by = "GENE"]

  for(gene in genes){
    if(verbose){
      print(paste0("Processing ", gene, " (", which(gene == genes), "/", length(genes), ")"))
    }
    subdat <- data[GENE == gene]
    subout <- subdat[, get(y)]
    n <- subdat[, .N]
    pteststat <- numeric(nperm)
    for(perm in 1:nperm){
      permout <- sample(subout, n, replace = F)
      oc1idx <- which(permout == outcomes[1])
      oc1mean <- subdat[oc1idx, mean(escape)]
      oc2mean <- subdat[!oc1idx, mean(escape)]
      pteststat[perm] <- abs(oc1mean - oc2mean) #Test stat for this perm
    }
    obsteststat <- abs(subdat[GENE == gene &  get(y)  == outcomes[1], mean(escape)] -
                         subdat[GENE == gene & get(y) == outcomes[2], mean(escape)])
    pteststat <- sort(pteststat)
    pval <- 1 - sum(obsteststat > pteststat)/nperm
    data[GENE == gene, Tobs := obsteststat]
    data[GENE == gene, perm_pvalue := pval]
  }
  ret <- unique(data[, list(GENE, nsamples, Tobs, perm_pvalue)])
  return(ret)
}
