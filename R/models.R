
# Summarizing options:
# - Sum across all SNPs
# - Retain only top SNP
# - Retain top <X> SNPs 
# - Sum only top <X> SNPs
# Pooling options:
# - independent samples
# - Select SNPs from multiple samples.
# - One estimate for multiple samples
# Return options
# - Same size as input. Can be useful for re-runs
# - Same size as number of calls
# - One call per SNP?
fitXCIR <- function(dt_anno, selection = c("topsnp", "sum"), nsnps = 1,
                    model = "AUTO", plot = FALSE, hist = FALSE,
                    flag = 0, xciGenes = NULL, a0 = NULL,
                    optimizer = c("nlminb", "optim"), method = NULL, limits = TRUE,
                    return = c("full", "calls"),
                    keep_params = FALSE, debug = FALSE){
  # getGenicDP
  snps <- copy(dt_anno)
  if(!"sex" %in% names(snps)){
    warning("There are no 'sex' column in the dataset. If not all subjects
            are female, remove them or add a 'sex' column to indicate 'female'
            samples and run again.")
    snps <- snps[, sex := "female"]
  }
  snps <- snps[sex == "female"]
  snps <- snps[, n_snps := .N, by = c("CHROM", "sample", "GENE")] 
  
  sel_snps <- select_snps(snps)
  
  # betaBinomXI
  dt[, xg := pmin(AD_hap1, AD_hap2)]
  dt[, ng := AD_hap1 + AD_hap2]
  
  training_genes <- readXCI(xciGenes)
  dt_train <- dt[GENE %in% training_genes]
  if(nrow(dt_train) == 0){
    warning("No known silenced gene found in data.")
  }
  
  optimizer <- match.arg(optimizer)
  model <- .check_model(model)
  modl <- vector("list", length(model))
  for(i in seq_along(modl)){
    modi <- model[i]
    if(modi == "M0"){
      dt <- M0(dt_train, dt, a0, optimizer, method, limits, debug)
    } else if(modi == "M1"){
      dt <- M1(dt_train, dt, a0, optimizer, method, limits, debug)
    } else if(modi == "M2"){
      dt <- M2(dt_train, dt, a0, optimizer, method, limits, debug)
    } else if(modi == "MF"){
      dt <- MF(dt_train, dt, a0, optimizer, method, limits, debug)
    }
    modl[[i]] <- dt
  }
  if(length(modl) > 1){
    dt <- .back_sel(modl, flag = flag, keep_params = keep_params)
  } else{
    dt <- modl[[1]]
  }
  
  dt[, f := a_est/(a_est + b_est)]
  dt[, fg := dp1/tot]
  
  # Use the estimated a and b to calculate var_fg and the test statistic
  # Under H0, f_g ~ M0
  dt[, var_fg := (tot * a_est * b_est * (a_est + b_est + tot))]
  dt[, var_fg := var_fg/((a_est + b_est)^2 * (a_est + b_est + 1))]
  dt[, var_fg := var_fg/tot^2] # Because we want the variance for the fraction, not the variance for the counts
  
  dt[, t := (fg-f)/sqrt(var_fg)] #Test statistic
  dt[, pbb := .pbb_midp(dp1, tot, a_est, b_est), by = c("GENE", "sample")]
  dt[, p_value := pnorm(t, lower.tail = FALSE)]
  dt[, status := ifelse(p_value < 0.05, "E", "S")]
  
  #tau is the Xi expression
  dt[, tau := (f-fg)/(2*f-1)]
  # tau is a fraction (0<tau<1) but estimate can be < 0 if fg < f
  # Which could only happens if a gene is tightly inactivated
  dt[, tau := ifelse(tau < 0, 0, tau)]
  dt[, var_tau := (2*f-1)^-2 * var_fg]
  dt[, ivw_tau := tau/sqrt(var_tau)]
  
  # Return
  return <- match.arg(optimizer)
  if(return == "full"){
    ret <- merge(dt_an)
  } else if(return == "calls"){
    
  }
}

.select_snps <- function(snps, selection, nsnps){
  
}