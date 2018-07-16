#' Fit mixture model
#'
#' Fit a mixture model to estimate mosaicism and XCI-escape.
#'
#' @param genic_dt A \code{data.table}. The table as outputted by \code{getGenicDP}.
#' @param model A \code{character} indicating which model to use to estimate
#'  the mosaicism. Valid choices are "AUTO", "BB", "MM", "MM2", "MM3". See details.
#' @param plot A \code{logical}. If set to TRUE, information about the training
#'  set and the skewing estimate will be plotted.
#' @param hist A \code{logical}. If set to TRUE, an histogram of the skewing
#'  estimates will be displayed.
#' @param flag A \code{numeric}. Specify how to handle convergence issues. See
#'  details.
#' @param xciGenes A \code{character} or NULL. To be passed to \code{readXCI} to
#'  select the training set of inactivated genes.
#' @param limits A \code{logical}. If set to TRUE, the optimization will be
#'  constrained. Using upper bounds on the probability of sequencing error and
#'  escape in the training set ensures that the dominant mixture represents the
#'  skewing for inactivated genes.
#'
#' @details
#' The model determines the number of components used in the mixture model. By
#' default, "AUTO" tries all combinations of mixtures and the best estimate is
#' kept using backward selection based on AIC.
#' BB is a simple beta-binomial. MM adds a binomial component to model the
#' sequencing errors. MM2 jointly models the probability of misclasification
#' in the training set. MM3 include all 3 components.
#'
#' Flags in the output reports issues in convergence. If \code{flag} is set to 0,
#' nothing is done. If set to 1, the model selection will avoid flagged models
#' (will favor parcimonious models).
#' If set to 2, calls for which the best selected model had convergence issue
#' will be removed.
#'
#' @seealso getGenicDP readXCI
#' @export
betaBinomXI <- function(genic_dt,  model = "AUTO", plot = FALSE, hist = FALSE,
                        flag = 0, xciGenes = NULL, limits = TRUE){
  dt <- copy(genic_dt)
  dt[, dp1 := pmin(AD_hap1, AD_hap2)]
  dt[, tot := AD_hap1 + AD_hap2]

  # Use XCI to estimate parameters
  #if(xci > 0){
  xcig <- readXCI(xciGenes)
  inactivated_genes <- xcig
  dt_xci <- dt[GENE %in% inactivated_genes]
  if(nrow(dt_xci) == 0){
    warning("No known silenced gene found in data.")
  }

  ## TODO: This does not really need to be passed and could be handled at the subfunction level
  balanced <- unique(dt_xci$sample)

  model <- .check_model(model)
  modl <- vector("list", length(model))
  for(i in seq_along(modl)){
    modi <- model[i]
    if(modi == "BB"){
      dt <- BB(dt_xci, dt, balanced, limits)
    } else if(modi == "MM"){
      dt <- MM(dt_xci, dt, balanced, limits)
    } else if(modi == "MM2"){
      dt <- MM2(dt_xci, dt, balanced, limits)
    } else if(modi == "MM3"){
      dt <- MM3(dt_xci, dt, balanced, limits)
    }
    modl[[i]] <- dt
  }
  if(length(modl) > 1){
    dt <- .back_sel(modl, flag = flag)
  } else{
    dt <- modl[[1]]
  }
  #if(model == "BB"){
  #  dt <- BB(dt_xci, dt, balanced, limits)
  #} else if(model == "MM"){
  #  dt <- MM(dt_xci, dt, balanced, limits)
  #} else if(model == "MM2"){
  #  dt <- MM2(dt_xci, dt, balanced, limits)
  #} else if(model == "MM3"){
  #  dt <- MM3(dt_xci, dt, balanced, limits)
  #} else if(model == "AUTO"){
  #  bb <- BB(dt_xci, dt, balanced, limits)
  #  mm <- MM(dt_xci, dt, balanced, limits)
  #  mm2 <- MM2(dt_xci, dt, balanced, limits)
  #  mm3 <- MM3(dt_xci, dt, balanced, limits)

  #  aics <- .back_sel(bb, mm, mm2, mm3, flag = flag)
  #  dt <- aics
  #} else{
  #  stop("Invalid model selected")
  #}

  dt[, f := a_est/(a_est + b_est)]
  dt[, fg := dp1/tot]

  # Use the estimated a and b to calculate var_fg and the test statistic
  # Under H0, f_g ~ BB
  dt[, var_fg := (tot * a_est * b_est * (a_est + b_est + tot))]
  dt[, var_fg := var_fg/((a_est + b_est)^2 * (a_est + b_est + 1))]
  dt[, var_fg := var_fg/tot^2] # Because we want the variance for the fraction, not the variance for the counts
  dt[, t := (fg-f)/sqrt(var_fg)] #Test statistic
  dt[, pbb := pbb(dp1, tot, a_est, b_est, type = "midp"), by = c("GENE", "sample")]
  dt[, p_value := pnorm(t, lower.tail = F)]

  #tau is the Xi expression
  dt[, tau := (f-fg)/(2*f-1)]
  dt[, var_tau := (2*f-1)^-2 * var_fg]
  dt[, ivw_tau := tau/sqrt(var_tau)]

  if(hist){
    hist(unique(dt[, f]), breaks = seq(0, .5, .01), main = "Cell fraction", xlab="Cell fraction", ylab = "samples")
  }
  if(plot){
    plotBBCellFrac(xci_dt = dt[GENE %in% xcig])
  }
  return(dt[])
}

.check_model <- function(model){
  model <- unique(toupper(model))
  validmodels <- c("BB", "MM", "MM2", "MM3")
  if("AUTO" %in% model){
    model <- validmodels
  }
  if(any(!model %in% validmodels)){
     err <- paste("Invalid models:", paste(model[!model %in% validmodels], collapse=", "),
                  "should be one of", paste(validmodels, collapse=", "))
     stop(err)
  }
  return(model)
}

################################################################################
dbb <- function(x, n, a, b){#Beta binomial
    #y <- choose(n, x) * beta(x+a, n-x+b)/beta(a,b)
    c <- choose(n, x)
    bn <- beta(x+a, n-x+b)
    bd <- beta(a, b)
    y <- c * bn / bd
    #y[is.nan(y)] <- 0
    ninf <- sum(is.infinite(y))
    if(ninf > 0){
      warning(paste("dbb returns",ninf, "infinite values"))
    }
    y <- y[is.finite(y)]
    return(y);
}

# Calculate log likelihood directly to avoid computation overflow
ldbb <- function(x, n, a, b){
  lc <- lchoose(n, x)
  lbn <- lbeta(x+a, n-x+b)
  lbd <- lbeta(a, b)
  ly <- lc + lbn - lbd
  ninf <- sum(is.infinite(ly))
  if(ninf > 0){
    warning(paste("ldbb returns",ninf, "infinite values"))
  }
  ly <- ly[is.finite(ly)]
  return(ly)
}

# P-values for exact inference
pbb <- function(x, n, a, b, type = "midp"){ #p-value for exact inference
  if(is.na(x) | is.na(n)){
    return(as.numeric(NA))
  }
  if(x >= n)
    return(0)
  if(type == "gt"){
    from <- x+1
    sump <- sum(dbb(from:n, n, a, b), na.rm = T)
  } else if(type == "geq"){
    from <- x
    sump <- sum(dbb(from:n, n, a, b), na.rm = T)
  } else if(type == "midp"){
    from <- x+1
    t0 <- dbb(x, n, a, b)
    t1 <- sum(dbb(from:n, n, a, b), na.rm = T)
    sump <- .5*t0 + t1
  } else{
    stop("Type should be one of 'gt', 'geq', 'midp'")
  }
  return(sump)
}

################################################################################
# Beta-binomial model:

BB <- function(dt_xci, full_dt, balanced, limits = F){
  dt <- copy(full_dt)
  for(sample_i in balanced){
    a0  <- c(1,1);
    dp1 <- dt_xci[sample == sample_i, dp1]
    dp  <- dt_xci[sample == sample_i, tot]
    Nxcig <- dt_xci[sample == sample_i, .N]
    if(limits){
      res_optim <- nlminb(a0, .logL_BB, dp1 = dp1, dp = dp,
                          lower = c(0, 0), upper = c(Inf, Inf))
    } else{
      res_optim <- nlminb(a0, .logL_BB, dp1 = dp1, dp = dp)
    }
    dt[sample == sample_i, a_est := exp(res_optim$par[1])]
    dt[sample == sample_i, b_est := exp(res_optim$par[2])]
    dt[sample == sample_i, model := "BB"]
    dt[sample == sample_i, k := length(res_optim$par)]
    dt[sample == sample_i, logL := res_optim$objective] #-logL
    dt[sample == sample_i, flag := res_optim$message]
    dt[sample == sample_i, AIC := 2*k + 2*logL]
    dt[sample == sample_i, Ntrain := Nxcig]
    dt[sample == sample_i, BIC := log(Ntrain)*k + 2*logL]
  }
  return(dt)
}
.logL_BB <- function(a, dp1, dp){
  a1 <- exp(a[1])
  b1 <- exp(a[2])
  logLbb <- sum(ldbb(dp1,dp,a1,b1)*(-1), na.rm = T)
  # logLbb <- sum(log(dbb(dp1,dp,a1,b1))*(-1), na.rm = T)
  return(logLbb)
}
################################################################################
# Mixture models:

# 1 Binom for sequencing errors and 1 BB for inactivated heterozygous SNP
MM <- function(dt_xci, full_dt, balanced, limits = F){
  dt <- copy(full_dt)
  for(sample_i in balanced){
    a0  <- c(1, 1, log(0.98/0.02), log(0.02/0.98))
    dp1 <- dt_xci[sample == sample_i, dp1]
    dp  <- dt_xci[sample == sample_i, tot]
    Nxcig <- dt_xci[sample == sample_i, .N]
    if(limits){
      res_optim <- nlminb(a0, .logL_MM, dp1 = dp1, dp = dp,
                          lower = c(0, 0, 0, -Inf),
                          #lower = c(-Inf, -Inf, 0, -Inf),
                          upper = c(Inf, Inf, Inf, log(0.2/0.8)))
    } else{
      res_optim <- nlminb(a0, .logL_MM, dp1 = dp1, dp = dp)
    }
    dt[sample == sample_i, a_est := exp(res_optim$par[1])]
    dt[sample == sample_i, b_est := exp(res_optim$par[2])]
    dt[sample == sample_i, p_het := exp(res_optim$par[3])/(1 + exp(res_optim$par[3]))]
    dt[sample == sample_i, pi_err := exp(res_optim$par[4])/(1 + exp(res_optim$par[4]))]
    dt[sample == sample_i, model := "MM"]
    dt[sample == sample_i, k := length(res_optim$par)]
    dt[sample == sample_i, logL := res_optim$objective] #Actually -logL
    dt[sample == sample_i, flag := res_optim$message]
    dt[sample == sample_i, AIC := 2*k + 2*logL]
    dt[sample == sample_i, Ntrain := Nxcig]
  }
  return(dt)
}

.logL_MM <- function(a, dp1, dp) {
  a1 <- exp(a[1]);
  b1 <- exp(a[2]);
  # Proportions have to be kept as an exponent ratio to avoid issues with boundaries
  p_het <- exp(a[3])/(1+exp(a[3])); # Proba that the SNP is indeed bi-allelic (initial proba = .98)
  pi_err <- 0.001+0.999*exp(a[4])/(1+exp(a[4])); # (initial proba = .02)
  pb <- dbinom(dp1, dp, pi_err)

  # Operations on the log scale to avoid infinite values
  lpbb <- ldbb(dp1, dp, a1, b1)
  p_tot <- p_het * exp(lpbb) + (1-p_het)*pb

  logL  <- -sum(log(p_tot))
  return(logL);
}



# 1 BB for inactivated SNP and 1 BB for escaped SNP
MM2 <- function(dt_xci, full_dt, balanced, limits = F){
  dt <- copy(full_dt)
  for(sample_i in balanced){
    a0 <- c(1, 1, #Beta for the inactivated genes mean = a/(a+b) = .5
            2, 2, #Beta for the escape  #2 has more density around 0.5
            log(0.9/0.1)) #Initial proba that a training gene is indeed inactivated in this sample
    dp1 <- dt_xci[sample == sample_i, dp1]
    dp  <- dt_xci[sample == sample_i, tot]
    Nxcig <- dt_xci[sample == sample_i, .N]
    if(limits){
      #ai,bi,ae,be > .5 to have a mode. p_inac > log(0.75/0.25) 3/4 of the training genes should be inactivated
      res_optim <- nlminb(a0, .logL_MM2, dp1 = dp1, dp = dp,
                          lower = c(0, 0, -Inf, -Inf, log(0.75/0.25)), # 0 mins min(p_inac) > .5
                          upper = c(Inf, Inf, Inf, Inf, Inf))
    } else{
      res_optim <- nlminb(a0, .logL_MM2, dp1 = dp1, dp = dp)
    }
    #p_inac_i <- exp(res_optim$par[5])/(1 + exp(res_optim$par[5]))
    expar5 <- exp(res_optim$par[5])
    if(is.finite(expar5)){
      p_inac_i <- expar5/(1 + expar5)
    } else if(expar5 == Inf){
      p_inac_i <- 1
    } else{
      stop("Something went wrong with the estimation of the inactivated mixture proportion")
    }
    # Selecting the right component to get alpha/beta corresponding to the inactivated group
    #if(p_inac_i >= .5){ # We assume that the list should always have a majority of inactivated genes
    #  a_est_i <- exp(res_optim$par[1])
    #  b_est_i <- exp(res_optim$par[2])
    #  a_est_e <- exp(res_optim$par[3])
    #  b_est_e <- exp(res_optim$par[4])
    #} else{
    #  p_inac_i <- 1-p_inac_i
    #  a_est_i <- exp(res_optim$par[3])
    #  b_est_i <- exp(res_optim$par[4])
    #  a_est_e <- exp(res_optim$par[1])
    #  b_est_e <- exp(res_optim$par[2])
    #}
    a_est_i <- exp(res_optim$par[1])
    b_est_i <- exp(res_optim$par[2])
    a_est_e <- exp(res_optim$par[3])
    b_est_e <- exp(res_optim$par[4])
    flagi <- res_optim$message
    fi <- a_est_i/(a_est_i + b_est_i)
    fe <- a_est_e/(a_est_e + b_est_e)
    if(!is.finite(fi) | !is.finite(fe)){
      stop(paste("For sample", sample_i, "fi is", fi, "and fe is", fe))
    }
    if(fi > fe){
      flagi <- "MM2: fi>fe"
    }
    #  # If the escape mixture has a lower average than the inactive, this is probably not the correct model
    #  dt[sample == sample_i, a_est := NA]
    #  dt[sample == sample_i, b_est := NA]
    #  dt[sample == sample_i, p_inac := NA]
    #  dt[sample == sample_i, p_het := NA]
#
    #  dt[sample == sample_i, pi_escape := NA]
    #  dt[sample == sample_i, k := NA]
    #  dt[sample == sample_i, logL := NA] #Actually -logL
    #  dt[sample == sample_i, AIC := NA]
    #  dt[sample == sample_i, Ntrain := Nxcig]
    #} else{
      dt[sample == sample_i, a_est := a_est_i]
      dt[sample == sample_i, b_est := b_est_i]
      dt[sample == sample_i, p_inac := p_inac_i]

      dt[sample == sample_i, pi_escape := a_est_e/(a_est_e + b_est_e)]
      dt[sample == sample_i, model := "MM2"]
      dt[sample == sample_i, k := length(res_optim$par)]
      dt[sample == sample_i, logL := res_optim$objective] #Actually -logL
      dt[sample == sample_i, flag := flagi]
      dt[sample == sample_i, AIC := 2*k + 2*logL]
      dt[sample == sample_i, Ntrain := Nxcig]
    #}
  }
  return(dt)
}
.logL_MM2 <- function(a, dp1, dp) {
  # Inactivated gene from XCI
  ai <- exp(a[1]);
  bi <- exp(a[2]);
  # Escape gene in XCI list
  ae <- exp(a[3])
  be <- exp(a[4])

  p_inac <- exp(a[5])/(1+exp(a[5])) #Proba that the gene is indeed inactivated (0.9)
  lpbbi <- ldbb(dp1, dp, ai, bi)
  lpbbe <- ldbb(dp1, dp, ae, be)
  p_tot <- p_inac * exp(lpbbi) + (1-p_inac) * exp(lpbbe)
  logL  <- -sum(log(p_tot))
  return(logL);
}

# 1 BB for inac SNPs 1 BB for escaped SNP and 1 Binom for sequencing err
MM3 <- function(dt_xci, full_dt, balanced, limits = F){
  dt <- copy(full_dt)
  for(sample_i in balanced){
    a0 <- c(1, 1, #Beta for the inactivated genes mean = a/(a+b) = .5
            2, 2, #Beta for the escape  #2 has more density around 0.5
            log(0.9/0.1), #Initial proba that a training gene is indeed inactivated in this sample
            log(0.98/0.02), #Proba of SNP being actually het
            log(0.02/0.98))
    dp1 <- dt_xci[sample == sample_i, dp1]
    dp  <- dt_xci[sample == sample_i, tot]
    Nxcig <- dt_xci[sample == sample_i, .N]
    err <- try({
      if(limits){
        res_optim <- nlminb(a0, .logL_MM3, dp1 = dp1, dp = dp,
                            lower = c(0, 0, -Inf, -Inf, 0, 0, -Inf),
                            upper = c(Inf, Inf, Inf, Inf, Inf, Inf, log(0.2/0.8)))
      } else{
        res_optim <- nlminb(a0, .logL_MM3, dp1 = dp1, dp = dp)
      }
      expar5 <- exp(res_optim$par[5])
      if(is.finite(expar5)){
        p_inac_i <- expar5/(1 + expar5)
      } else if(expar5 == Inf){
        p_inac_i <- 1
      } else{
        stop("Something went wrong with the estimation of the inactivated mixture")
      }
      # if(p_inac_i >= .5){ # We assume that the list should always have a majority of inactivated genes
      #   a_est_i <- exp(res_optim$par[1])
      #   b_est_i <- exp(res_optim$par[2])
      #   a_est_e <- exp(res_optim$par[3])
      #   b_est_e <- exp(res_optim$par[4])
      # } else{
      #   p_inac_i <- 1-p_inac_i
      #   a_est_i <- exp(res_optim$par[3])
      #   b_est_i <- exp(res_optim$par[4])
      #   a_est_e <- exp(res_optim$par[1])
      #   b_est_e <- exp(res_optim$par[2])
      # }
      a_est_i <- exp(res_optim$par[1])
      b_est_i <- exp(res_optim$par[2])
      a_est_e <- exp(res_optim$par[3])
      b_est_e <- exp(res_optim$par[4])
      flagi <- res_optim$message
      fi <- a_est_i/(a_est_i + b_est_i)
      fe <- a_est_e/(a_est_e + b_est_e)
      ferr_i <- exp(res_optim$par[7])/(1 + exp(res_optim$par[7]))
      if(fi > fe){
        flagi <- "MM3: fi>fe"
      } else if(fi < ferr_i){
        flagi <- "MM3: ferr>fi"
      }
    })
    if(inherits(err, "try-error")){
      dt[sample == sample_i, a_est := NA]
      dt[sample == sample_i, b_est := NA]
      dt[sample == sample_i, p_inac := NA]
      dt[sample == sample_i, p_het := NA]

      dt[sample == sample_i, pi_escape := NA]
      dt[sample == sample_i, model := "MM3"]
      dt[sample == sample_i, k := NA]
      dt[sample == sample_i, logL := NA] #Actually -logL
      dt[sample == sample_i, flag := "MM3: Error"]
      dt[sample == sample_i, AIC := NA]
      dt[sample == sample_i, Ntrain := Nxcig]
    } else{
      dt[sample == sample_i, a_est := a_est_i]
      dt[sample == sample_i, b_est := b_est_i]
      dt[sample == sample_i, p_inac := p_inac_i]
      dt[sample == sample_i, p_het := exp(res_optim$par[6])/(1+exp(res_optim$par[6]))]

      dt[sample == sample_i, pi_escape := a_est_e/(a_est_e + b_est_e)]
      dt[sample == sample_i, pi_err := ferr_i]
      dt[sample == sample_i, model := "MM3"]
      dt[sample == sample_i, k := length(res_optim$par)]
      dt[sample == sample_i, logL := res_optim$objective] #Actually -logL
      dt[sample == sample_i, flag := flagi]
      dt[sample == sample_i, AIC := 2*k + 2*logL]
      dt[sample == sample_i, Ntrain := Nxcig]
    }
  }
  return(dt)
}


.logL_MM3 <- function(a, dp1, dp){
  # Inactivated gene from XCI
  ai <- exp(a[1])
  bi <- exp(a[2])
  # Escape gene in XCI list
  ae <- exp(a[3])
  be <- exp(a[4])

  if(is.finite(exp(a[5]))){
    p_inac <- exp(a[5])/(1+exp(a[5])) #Proba that the gene is not a training set error
  } else if(exp(a[5]) == Inf){
    p_inac <- 1 #Tends to 1
  } else{
    stop("Something is wrong with the a[5]")
  }
  if(is.finite(exp(a[6]))){
    p_het <- exp(a[6])/(1+exp(a[6])) #Proba that the gene is indeed heterozygous
  } else if(exp(a[6]) == Inf){
    p_het <- 1 #Tends to 1
  } else{
    stop("Something is wrong with the a[6]")
  }

  pi_err  <- .01 + .99*exp(a[7])/(1+exp(a[7]))

  lpbb_i <- ldbb(dp1, dp, ai, bi)
  lpbb_e <- ldbb(dp1, dp, ae, be)
  pb_err <- dbinom(dp1, dp, pi_err)
  ## TODO: Handle cases where exp lpbb_i/lpbb_e > 0
  ## lpbb should be between 0 and 1
  p_tot <- p_het * (p_inac * exp(lpbb_i) + # Reads from component1
                   (1-p_inac) * exp(lpbb_e)) + # Reads from component2
           (1-p_het) * pb_err # Reads from sequencing error

  logL <- -sum(log(p_tot))
  return(logL)
}
################################################################################
# Model selection
.back_sel <- function(modl, criterion = "AIC", flag = 0){
  cols <- Reduce(intersect, lapply(modl, names))
  modl <- lapply(modl, function(XX){XX[, cols, with = F]})
  aics <- rbindlist(modl)
  if(flag == 1){# Models where the sample had errors are discarded
    aics <- aics[!grep("^MM", flag)]
  }
  aics <- aics[aics[, .I[which.min(AIC)], by = "sample,GENE"]$V1] #Model that minimizes AIC
  if(flag == 2){# Remove samples where the selected model had errors
    aics <- aics[!grep("^MM", flag)]
  }
  return(aics)
}
################################################################################

# Can also be used to visualize the cell frac from the calls
#' Plot cell fraction estimates
#'
#' Plot cell fraction estimates from list of known XCI genes
#'
#' @param xci_dt A \code{data.table}. The data to be used for the estimate
#' of skewing (i.e: limited to XCI genes).
#' @param xcig A \code{logical}. If \code{xci_dt} was not subset for training
#' genes only, setting xcig to TRUE will filter the data.
#' @param gene_names A \code{character}. If left blank, only genes that are
#' further than 20% away from the estimated skewing will be annotated, if set
#' to "all", all genes will be named. Set to "none" to remove all annotations.
#' Alternately, a \code{character} vector can be passed to annotate specific
#' genes of interest.
#'
#' @return NULL
#'
#' @details
#' This function is mostly used in \code{betaBinomXI} to ensure that the cell
#' fraction is estimated properly. However, it can be used from the output
#' of \code{betaBinomXI} to troubleshoot estimation issues.
#'
#' @importFrom ggplot2 element_blank scale_colour_manual facet_wrap
#' @importFrom ggplot2 geom_point geom_hline theme_bw
#' @export
plotBBCellFrac <- function(xci_dt, xcig = NULL, gene_names = ""){
  plotfrac <- xci_dt[order(sample)]
  if(!is.null(xcig)){
    xcig <- readXCI(xcig)
    plotfrac <- plotfrac[GENE %in% xcig]
  }
  Nt <- plotfrac[, .N, by = "sample"]
  plotfrac <- merge(plotfrac, Nt, by = "sample")
  plotfrac[, index := unlist(sapply(Nt[, N], function(x){1:x}))]
  plotfrac[, label := ""]
  if(length(gene_names) > 1){
    plotfrac[GENE %in% gene_names, label := GENE]
  } else if(gene_names == "all"){
    plotfrac[, label := GENE]
  } else if(gene_names == "none"){
    plotfrac[, label := ""]
  } else if(gene_names == ""){
    plotfrac[abs(fg - f) > .2, label := GENE]
  } else{
    plotfrac[GENE %in% gene_names, label := GENE]
  }
  if("model" %in% colnames(plotfrac)){
    gp <- geom_point(aes(shape = model))
  } else{
    gp <- geom_point()
  }

  p <- ggplot(plotfrac, aes(x = index, y = fg))
  ## TODO: Add the other components when available (also add option to return components from all models to betaBinomXI)

  p <- p + geom_hline(aes(yintercept = f))
  p <- p + gp + geom_text(aes(label = label))
  #p <- p + geom_text(aes(x = N - 5, y= max(fg) + .01, label = paste("N =", N)))
  p <- p + geom_text(aes(x = 0, y= max(fg) + .01, label = paste("N =", N)), hjust = "left")
  p <- p + facet_wrap(~sample, scales = "free_x")
  p <- p + theme_bw() + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x = element_blank())
  print(p)
  return(invisible(p))
}

