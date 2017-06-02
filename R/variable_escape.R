# # When the truth is unknown, plot the % of escape for every gene
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 scale_y_continuous scale_fill_manual ggtitle element_text
#' @export
plot_status <- function(sub_xi, alpha = .05, min_sup = 0, rownames = NULL,
                        inference = "asymptotic", nsamples =F, threshold = F){
  ns <- sub_xi[, .N, by = GENE] #Number of samples per gene
  if(inference == "asymptotic"){
    pe <- sub_xi[p_value < alpha, .N, by = GENE] #Samples where it escapes
  } else if(inference == "exact"){
    pe <- sub_xi[pbb < alpha, .N, by = GENE] #Samples where it escapes
  } else{
    stop("Inference should be 'asymptotic' or 'exact'")
  }
  dat <- merge(ns, pe, by = "GENE", all.x=T)
  dat <- dat[is.na(N.y), N.y := 0]
  dat <- merge(dat, unique(sub_xi[, list(GENE, POS)]), by = "GENE", all.x=T) #POS are needed to order the X axis
  setnames(dat, c("N.x", "N.y"), c("N_support", "N_escape"))
  dat[, pe := N_escape/N_support]

  Nseq <- seq(0, max(dat$N_support) + 10, 10) #Bin the number of support genes
  palette <- brewer.pal(3, "Set1")
  colfunc <- colorRampPalette(c(palette[2], palette[1]))
  ncolors <- length(unique(cut(dat$N_support, Nseq)))
  colors <- colfunc(ncolors)

  dat <- dat[N_support >= min_sup] #Also removes the genes with %e equal to 0 that dont have much support, as intended.
  dat <- dat[order(POS)]
  dat <- unique(dat[, list(GENE, N_support, N_escape, pe)])
  dat[, GENE := factor(GENE, levels = GENE)]
  p <- ggplot(data = dat)
  p <- p + geom_bar(aes(x = GENE, y = pe, fill = cut(N_support, Nseq)), stat = "identity")
  p <- p + scale_fill_manual(values = colors, name = "Number of support samples")# + ylim(0, 1) #ylim throws warning if used with an extra scale
  status_theme <- theme(axis.text.x = element_text(angle = 90, size = 7, vjust = .5),
                        plot.title = element_text(hjust=0.5))
  status_scale <- scale_y_continuous(limits = c(0, 1), expand = c(0.003, 0))
  if(nsamples){
    p <- p+geom_text(aes(label=N_support, x=GENE, y=0.05))
  }
  if(is.null(rownames)){
    ngenes <- length(unique(dat$GENE))
    if(ngenes > 100){
      status_theme <- status_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
  } else if(!rownames){
    status_theme <- status_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  p <- p + status_theme + status_scale + ylab("%escape") + ggtitle("Percentage of escape accross samples")
  print(p)

  return(dat)
}


#' Plot the fraction of genes in each category
#'
#' Barplot of the fraction of gene that are 'subject', 'variable' and 'escape'
#' inactivation from the calls made in XCIR
#'
#' @param xi A \code{data.table}. The calls, output from \code{betaBinomXI}.
#' @param alpha A \code{numeric}. The significance level.
#' @param inference A \code{character}. Which p-value should be used to
#'  categorize each gene. Valid choices are 'asymptotic' and 'exact'.
#' @param threshold A \code{character} of length 2. Genes with a percent of
#'  escape lower than the first number will be considered subject to
#'  inactivation. Genes with a higher percent of escape than the second number
#'  will be categorized as escaping inactivation.
#' @param min_support A \code{numeric}. Only calculate percentage of escape for
#'  genes that have at least that many samples for which a prediction was made.
#'
#' @return An invisible \code{data.table} with one row for each gene.
#'
#' @export
plot_status_fraction <- function(xi, alpha = .05, inference = "asymptotic", threshold = c(.25, .75), min_support = 0, plot = TRUE){
  dat <- get_status(xi, alpha = alpha, inference = inference)
  dat <- dat[N_support >= min_support]

  dat[, status := ifelse(pe <  threshold[1], "subject", "variable")]
  dat[, status := ifelse(pe >= threshold[2], "escape", status)]
  dat[, count := .N, by = status]
  dat[, percent := round(count/nrow(dat), 3)*100]
  dat[, status := factor(status, levels = c("subject", "variable", "escape"))]

  # Plot
  if(plot){
    p <- ggplot(dat) + geom_bar(aes(status, fill = status))
    p <- p + geom_text(aes(status, y = count, label = paste(percent, "%",  "\n", count)))
    p <- p + ggtitle("Escape status") + theme(plot.title = element_text(hjust=.5))
    print(p)
  }
  return(invisible(dat[]))
}

get_status <- function(xi, alpha, inference){
  ns <- xi[, .N, by = GENE] #Number of samples per gene
  if(inference == "asymptotic"){
    pe <- xi[p_value < alpha, .N, by = GENE] #Samples where it escapes
  } else if(inference == "exact"){
    pe <- xi[pbb < alpha, .N, by = GENE] #Samples where it escapes
  } else{
    stop("Inference should be 'asymptotic' or 'exact'")
  }
  dat <- merge(ns, pe, by = "GENE", all.x=T)
  dat <- dat[is.na(N.y), N.y := 0]
  setnames(dat, c("N.x", "N.y"), c("N_support", "N_escape"))
  dat[, pe := N_escape/N_support]
  return(dat)
}

#' @export
concordance <- function(xi, xciGenes = NULL, inference = "exact"){
  if(is.null(xciGenes))
    xcig <- readXCI()
  else
    xcig <- xciGenes
  skews <- c(15, 25, 50)
  l <- vector('list', length(skews))
  names(l) <- paste0("s", skews)
  for(skew in skews){
    par1   <- plot_status_fraction(xi[f < skew/100 & POS < 2699520], inference = inference, plot = F)
    xci177 <- plot_status_fraction(xi[f < skew/100 & GENE %in% xcig], inference = inference, plot = F)
    all    <- plot_status_fraction(xi[f < skew/100], inference = inference, plot = F)
    par1[, region := "Par 1"]
    xci177[, region := "XCI"]
    all[, region := "chrX"]
    s <- rbindlist(list(par1, xci177, all))
    s[, skew := skew]
    l[[paste0("s", skew)]] <- s
  }
  dat <- rbindlist(l)
  p <- ggplot(dat) + geom_bar(aes(status, fill = status)) +
    geom_text(aes(status, y = count, label = paste(percent, "%"))) +
    facet_grid(region ~ skew, scales = "free_y")
  print(p)
  return(dat)
}

#' Cochran-Mantel-Haenszel test
#'
#' Runs CMH on \code{XCIR} calls to test the association between XCI status
#' and an outcome of interest.
#'
#' @param calls A \code{character}. XCI status prediction as outputed by
#'  \code{betaBinomXI}.
#' @param condition A \code{character}. The column name that contains the
#'  outcomes of interest.
#' @param value A \code{character}. The outcome of interest. All other values
#'  in the \code{condition} columns will be matched against it.
#'
#' @note
#' The CMH will only be computed using genes that have at least 1 sample for
#' each outcome.
#'
#' @seealso \code{\link{betaBinomXI}}
#'
#' @export
CMH <- function(calls, condition, value){
  if(!condition %in% names(calls)){
    stop(paste(condition, "is not a column of the dataset"))
  }
  if(!value %in% calls[[condition]]){
    stop(paste(value, "does not exist in the ", condition, "column"))
  }
  ## TODO: Allow flexibility
  inference <- "asymptotic"
  if(tolower(inference) == "exact"){
    pval <- "pbb"
  } else{
    pval <- "p_value"
  }
  # sample, condition, gene, p_value
  cmh_dat <- calls[, list(t = .N,
                 m1 = sum(p_value > .05),  # Subjects accross conditions
                 m2 = sum(p_value <= .05), # Escape accross conditions
                 n1 = sum(get(condition) == value),
                 n2 = sum(get(condition) != value),
                 a = sum(p_value > .05 & source == value)),
          by="GENE"]
  gn_rm <- cmh_dat[n1 == 0 | n2 == 0, GENE]
  message(paste0(length(gn_rm), "/", nrow(cmh_dat), " genes were removed because they do not have at least one call for each condition."))
  cmh_dat <- cmh_dat[!GENE %in% gn_rm]

  cmh_dat[, num := a - (n1*m1)/t]
  cmh_dat[, den := (n1*n2*m1*m2)/(t^2*(t-1))]
  cmh_stat <- cmh_dat[, sum(num)^2/sum(den)] #Reject H0:R=1 for large values of cmh_stat
  # Compare to a chisq with 1 df
  cmh_pv <- pchisq(cmh_stat, 1, lower.tail = F)

  cmh_list <- list(cmh_stat = cmh_stat, cmh_pv = cmh_pv)
  print_cmh(cmh_list)
  return(cmh_list)
}

print_cmh <- function(cmh_list, alpha = .05){
  decision <- ifelse(cmh_list$cmh_pv < alpha, "rejects", "fails to reject")
  print(data.frame(X2 = cmh_list$cmh_stat, pvalue = cmh_list$cmh_pv), row.names = F)
  print(paste("The test", decision, "the null hypothesis that there is no association between XCI status and the outcome of interest."))
  return(NULL)
}

#' @importFrom pheatmap pheatmap
#' @export
res_heatmap <- function(XCIres, value="tau", condition=NULL, FC=FALSE){
  data <- dcast(XCIres, GENE ~ sample, value.var = value)
}
