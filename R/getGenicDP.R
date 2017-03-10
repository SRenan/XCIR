#' Get sum of the coverage of all SNPs for each gene
#'
#' @param pool A \code{logical}. Only works when \code{highest_expr} is set to
#'  TRUE. If set to TRUE, the read counts are pooled accross all samples for
#'  each SNP. Only use this if the samples come from the same subject
#'
#' @importFrom ggplot2 ggplot geom_histogram xlab ylab aes
#' @export
getGenicDP <- function(dt_anno, highest_expr = FALSE, pool = FALSE, gender_file = NULL){
  if(!is.null(gender_file)){
    sex <- fread(gender_file)
    setnames(sex, c("sample", "gender"))
    sex[, gender := tolower(gender)]
    dt_anno <- merge(dt_anno, sex[, list(sample, gender)], by = "sample")
  }
  else if(!"gender" %in% names(dt_anno)){
    stop("If there is no 'gender' column in the input, a gender_file mapping
         sample to gender should be provided")
  }
  dt_anno <- dt_anno[gender == "female"]
  dt_anno <- dt_anno[, n_snps := .N, by = c("CHROM", "sample", "GENE")]

  if(highest_expr){
    bygene <- getGenicMostExp(dt_anno, pool)
  } else{
    bygene <- getGenicAll(dt_anno)
  }
  return(bygene[])
}

getGenicMostExp <- function(annotated, pool = FALSE){
  dt_anno <- copy(annotated)
  if(pool){ # All samples come from the same subject
    summed <- dt_anno[, lapply(.SD, sum, na.rm = TRUE), by = c("CHROM", "POS", "GENE", "n_snps"), .SDcols = c("AD_hap1", "AD_hap2")]
    summed <- summed[, tot := AD_hap1 + AD_hap2]
    maxpos <- summed[summed[, .I[which.max(tot)], by = "GENE"]$V1]$POS
    bygene <- dt_anno[POS %in% maxpos, list(CHROM, POS, REF, ALT, GENE, sample, gender, n_snps, AD_hap1, AD_hap2)] #nrow = nGenes * nSubs
    bygene[, tot := AD_hap1 + AD_hap2]
  } else{
    high <- copy(dt_anno)
    high <- high[, tot := AD_hap1 + AD_hap2]
    bygene <- high[high[, .I[which.max(tot)], by = "sample,GENE"]$V1]
  }
  return(bygene)
}

getGenicAll <- function(annotated){
  bygene <- annotated[, lapply(.SD, sum, na.rm = TRUE), by = c("CHROM", "gender", "GENE", "sample", "n_snps"), .SDcols = c("AD_hap1", "AD_hap2")] #nrow = nGenes * nSubs
  bygene[, tot := AD_hap1 + AD_hap2]
  ua <- unique(annotated[, list(POS, ANNO, GENE, sample)])
  ua <- ua[ua[, .I[which.min(POS)], by = c("sample","GENE")]$V1]
  bygene <- merge(bygene, ua, by = c("GENE", "sample"))
  setcolorder(bygene, c("sample", "gender", "CHROM", "POS", "ANNO", "GENE", "AD_hap1", "AD_hap2", "tot", "n_snps"))
  return(bygene)
}
