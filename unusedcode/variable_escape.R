# # When the truth is unknown, plot the % of escape for every gene

#' Plot escape across X
#'
#' Barplot of the escape status of X-linked genes. Each bar represent the
#' proportions of samples in which a gene is found to escape.
#' Genes are ordered along the X chromosome.
#'
#' @param sub_xi A \code{data.table}. The output of \code{betaBinomXI}.
#' @param alpha A \code{numeric} between 0 and 1. The significance level for
#' escape calls.
#' @param min_sup A \code{numeric}. Only genes with at least that many samples
#' will be displayed.
#' @param rownames A \code{logical} or NULL. If set to NULL, gene names will
#' be displayed only if there are less than 100 genes on the plot.
#' @param inference A \code{character}. Which of the \code{betaBinomXI}
#' predictions should be used, "asymptotic" or "exact".
#' @param nsamples A \code{logical}. If TRUE, display the total number of samples
#' in the dataset.
#' @param cutoffs A \code{numeric} vector of length 2. The boundaries for
#' calling subject, variable and escape genes.
#' @param theme A \code{character}. Available options are "gray", "status" and
#' "samplecount".
#' @param hg A \code{numeric}. The genome version that was used to align the
#' RNA-Seq data. This is used to display annotations.
#'
#' @return An invisible \code{data.table} with the following columns
#' \item{GENE}{Gene symbol}
#' \item{N_support}{The number of samples for which a prediction was made for that gene}
#' \item{N_escape}{The number of samples where this gene was found to escape}
#' \item{pe}{The fraction of escape, N_escape/N_support}
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 scale_y_continuous scale_fill_manual ggtitle element_text
#' @importFrom ggplot2 geom_vline geom_bar
#' @importFrom ggplot2 aes theme
#' @importFrom grDevices colorRampPalette
#' @export
plot_status <- function(sub_xi, alpha = .05, min_sup = 0, rownames = NULL,
                        inference = "asymptotic", nsamples =F, cutoffs = c(25, 75),
                        theme = "gray", hg = 19){
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
  dat <- dat[order(POS)]

  plotdat <- dat[N_support >= min_sup] #Also removes the genes with %e equal to 0 that dont have much support, as intended.
  plotdat <- unique(plotdat[, list(GENE, N_support, N_escape, pe)])
  plotdat[, GENE := factor(GENE, levels = GENE)]
  p <- ggplot(data = plotdat, mapping = aes(x = GENE, y = pe))
  if(theme == "gray"){
    p <- p + geom_bar(aes(x = GENE, y = pe), stat = "identity")
  } else if(theme == "status"){
    xmax <- nrow(plotdat)
    p <- p + geom_rect(data = NULL, xmin = 0, xmax = 275, ymin = 0, ymax = .25, fill = "lightyellow") +
             geom_rect(data = NULL, xmin = 0, xmax = xmax, ymin = .25, ymax = .75, fill = "lightgreen") +
             geom_rect(data = NULL, xmin = 0, xmax = xmax, ymin = .75, ymax = 1, fill = "lightblue")
    p <- p + geom_bar(stat = "identity")
    xmaxlabel <- ifelse(xmax > 10, xmax -1, xmax)
    p <- p + geom_text(aes(label = "Escape", xmaxlabel, y = .76, hjust = 1, vjust = 0), color = "red", cex = 10) +
             geom_text(aes(label = "Variable", xmaxlabel, y = .5, hjust = 1),color = "red", cex = 10) +
             geom_text(aes(label = "X inactivated", xmaxlabel, y = .24, hjust = 1, vjust = 1), color = "red",cex = 10)
  } else{
    Nseq <- seq(0, max(dat$N_support) + 10, 10) #Bin the number of support genes
    palette <- brewer.pal(3, "Set1")
    colfunc <- colorRampPalette(c(palette[2], palette[1]))
    ncolors <- length(unique(cut(dat$N_support, Nseq)))
    colors <- colfunc(ncolors)
    p <- p + geom_bar(aes(x = GENE, y = pe, fill = cut(N_support, Nseq)), stat = "identity")
    p <- p + scale_fill_manual(values = colors, name = "Number of support samples")# + ylim(0, 1) #ylim throws warning if used with an extra scale
  }
  # Add vertical line to denote PAR1 region
  if(hg == 19){
    par1_l <- 60001
    par1_u <- 2699520
    par1_ug <- length(unique(dat[POS < par1_u, GENE]))
    if(par1_ug > 0){
      p <- p + geom_vline(xintercept = par1_ug + .5, color = "red") #+
               #geom_text(aes(label = "PAR1", x = par1_ug, y = 1.1, hjust = 1, vjust = 0), color = "red")
    }
  } else if(hg == 38){
    par1_l <- 10001
    par1_u <- 2781479
  }

  status_theme <- theme(axis.text.x = element_text(angle = 90, size = 7, vjust = .5),
                        plot.title = element_text(hjust=0.5))
  status_scale <- scale_y_continuous(limits = c(0, 1), expand = c(0.003, 0))
  if(nsamples){
    p <- p+geom_text(aes(label=N_support, x=GENE, y=0.05))
  }
  if(is.null(rownames)){
    ngenes <- length(unique(plotdat$GENE))
    if(ngenes > 100){
      status_theme <- status_theme + theme(axis.text.x = element_blank(),
                                           axis.ticks.x = element_blank(),
                                           panel.grid.major.x = element_blank())
    }
  } else if(!rownames){
    status_theme <- status_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  if(max(sub_xi$f) < .25){
    title <- "Percentage of escape across skewed samples (> 25/75)"
  } else{
    title <- "Percentage of escape across samples"
  }
  p <- p + status_theme + status_scale + ylab("%escape") + ggtitle(title)
  print(p)

  return(invisible(plotdat))
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
#' @param plot A \code{logical}. If set to FALSE, do not print the plot and only
#' return the (invisible) \code{data.table}.
#'
#' @return An invisible \code{data.table} with the following columns
#' \item{GENE}{Gene symbol}
#' \item{N_support}{The number of samples for which a prediction was made for that gene}
#' \item{N_escape}{The number of samples where this gene was found to escape}
#' \item{pe}{The fraction of escape, N_escape/N_support}
#'
#' @export
plot_status_fraction <- function(xi, alpha = .05, inference = "asymptotic",
                                 threshold = c(.25, .75), min_support = 0, plot = TRUE){
  dat <- .get_status(xi, alpha = alpha, inference = inference)
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

.get_status <- function(xi, alpha, inference){
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
#' @value A \code{list} with the test statistic and p-value for each gene.
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

#' Summarize escape frequency
#'
#'
#' @param table A \code{data.table}. The output of \code{betaBinomXI}.
#' @param fmax A \code{numeric}. The maximum skew (takes values between 0 and
#' 0.5, 0.5 corresponds to balanced samples).
#' @param alpha A \code{numeric}. The significance level of escape calls.
#'
#' @export
summescape <- function(table, fmax = 0.5, alpha = 0.05){
  summ <- copy(table[f <= fmax])
  summ[, status := ifelse(p_value < alpha, "E", "S")]
  summ <- setkey(summ, GENE, status)[CJ(GENE, status, unique = T), .N, by = .EACHI]
  summ <- dcast.data.table(summ, GENE ~ status, value.var = "N")
  summ[, N := E+S]
  summ[, esc := E/N]
  return(summ[])
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

