#' Get expression at the gene level
#'
#' Calculate allele specific expression for each gene in each sample, either
#' using only the most expressed SNP or using all SNPs (when phasing has been
#' performed).
#'
#' @param dt_anno A \code{data.table}. An annotated table of read counts for
#'  each SNP, as outputted by \code{addAnno}
#' @param highest_expr A \code{logical}. If FALSE, all SNPs will be
#'  summed within each gene. This should only be set to FALSE when high quality
#'  phasing information is available. If set to TRUE, the highest expressed SNP
#'  (across both alleles) will be used instead.
#' @param pool A \code{logical}. Only works when \code{highest_expr} is set to
#'  TRUE. If set to TRUE, the read counts are pooled accross all samples for
#'  each SNP. Only use this if the samples come from the same subject
#' @param sex_file A \code{character} or NULL. Leave NULL if \code{dt_anno}
#'  already contains a sex column. The file must contain at least a "sample"
#'  and "sex" column with samples matching the samples in \code{dt_anno}.
#'
#' @return A \code{data.table}. That should be used as input for
#'  \code{betaBinomXI}.
#'
#' @example inst/examples/workflow.R
#'
#' @seealso betaBinomXI, addAnno
#' @importFrom ggplot2 ggplot geom_histogram xlab ylab aes
#' @export
getGenicDP <- function(dt_anno, highest_expr = TRUE, pool = FALSE, sex_file = NULL){
  if(!is.null(sex_file)){
    sex <- fread(sex_file)
    setnames(sex, c("sample", "sex"))
    sex[, sex := tolower(sex)]
    dt_anno <- merge(dt_anno, sex[, list(sample, sex)], by = "sample")
  }
  else if(!"sex" %in% names(dt_anno)){
    warning("There are no 'sex' column in the dataset. If not all subjects
            are female, remove them or add a 'sex' column and run again.")
    dt_anno <- dt_anno[, sex := "female"]
  }
  dt_anno <- dt_anno[sex == "female"]
  dt_anno <- dt_anno[, n_snps := .N, by = c("CHROM", "sample", "GENE")]

  if(highest_expr){
    bygene <- .getGenicMostExp(dt_anno, pool)
  } else{
    bygene <- .getGenicAll(dt_anno)
  }
  return(bygene[])
}

.getGenicMostExp <- function(annotated, pool = FALSE){
  dt_anno <- copy(annotated)
  if(pool){ # All samples come from the same subject
    summed <- dt_anno[, lapply(.SD, sum, na.rm = TRUE), by = c("CHROM", "POS", "GENE", "n_snps"), .SDcols = c("AD_hap1", "AD_hap2")]
    summed <- summed[, tot := AD_hap1 + AD_hap2]
    maxpos <- summed[summed[, .I[which.max(tot)], by = "GENE"]$V1]$POS
    bygene <- dt_anno[POS %in% maxpos, list(CHROM, POS, REF, ALT, GENE, sample, sex, n_snps, AD_hap1, AD_hap2)] #nrow = nGenes * nSubs
    bygene[, tot := AD_hap1 + AD_hap2]
  } else{
    high <- copy(dt_anno)
    high <- high[, tot := AD_hap1 + AD_hap2]
    bygene <- high[high[, .I[which.max(tot)], by = "sample,GENE"]$V1]
  }
  return(bygene)
}

.getGenicAll <- function(annotated){
  bygene <- annotated[, lapply(.SD, sum, na.rm = TRUE), by = c("CHROM", "sex", "GENE", "sample", "n_snps"), .SDcols = c("AD_hap1", "AD_hap2")] #nrow = nGenes * nSubs
  bygene[, tot := AD_hap1 + AD_hap2]
  ua <- unique(annotated[, list(POS, GENE, sample)])
  ua <- ua[ua[, .I[which.min(POS)], by = c("sample","GENE")]$V1]
  bygene <- merge(bygene, ua, by = c("GENE", "sample"))
  setcolorder(bygene, c("sample", "sex", "CHROM", "POS", "GENE", "AD_hap1", "AD_hap2", "tot", "n_snps"))
  return(bygene)
}

#' @importFrom utils tail
.getTopSNPs <- function(annotated, nsnps = 1){
  annotated[, top := AD_hap1 + AD_hap2]
  setkeyv(annotated, "sample", "GENE", "tot")
  top <- annotated[, tail(.SD, nsnps), by = c("sample", "GENE")]
  return(top)
}

