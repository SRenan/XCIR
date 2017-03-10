# # When the truth is unknown, plot the % of escape for every gene
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 scale_y_continuous
#' @export
plot_status <- function(sub_xi, alpha = .05, min_sup = 0, rownames = T, inference = "asymptotic"){
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
  if(!rownames){
    status_theme <- status_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  p <- p + status_theme + status_scale + ggtitle("Percentage of escape accross samples")
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
plot_status_fraction <- function(xi, alpha = .05, inference = "asymptotic", threshold = c(.25, .75), min_support = 0){
  dat <- get_status(xi, alpha = alpha, inference = inference)
  dat <- dat[N_support >= min_support]

  dat[, status := ifelse(pe <  threshold[1], "subject", "variable")]
  dat[, status := ifelse(pe >= threshold[2], "escape", status)]
  dat[, count := .N, by = status]
  dat[, percent := round(count/nrow(dat), 3)*100]
  dat[, status := factor(status, levels = c("subject", "variable", "escape"))]

    # Plot
  p <- ggplot(dat) + geom_bar(aes(status, fill = status))
  p <- p + geom_text(aes(status, y = count, label = paste(percent, "%",  "\n", count)))
  p <- p + ggtitle("Escape status") + theme(plot.title = element_text(hjust=.5))
  print(p)
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

concordance <- function(xi, xciGenes = NULL){
  if(is.null(xciGenes))
    xcig <- readXCI()
  else
    xcig <- xciGenes
  skews <- c(15, 25, 50)
  l <- vector('list', length(skews))
  names(l) <- paste0("s", skews)
  for(skew in skews){
    par1   <- plot_status_fraction(xi[f < skew/100 & POS < 2699520])
    xci177 <- plot_status_fraction(xi[f < skew/100 & GENE %in% xcig])
    all    <-   plot_status_fraction(xi[f < skew/100])
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
