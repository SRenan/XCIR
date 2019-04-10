
#' @importFrom stringr str_extract
#' @export
read_AD_format <- function(file, het_cutoff){
  fr <- fread(file)
  re <- "^\\d+,\\d+" #First two numerics separated by a comma (ignore additional ALTs which are always 0 AFAIK)
  mfr <- melt(fr, id.vars = c("CHROM", "POS"), variable.name = "sample", value.name = "AD")
  mfr[, ADs := str_extract(AD, re)] #102 seconds
  mfr <- mfr[ADs != "0,0"]
  mfr[, AD_hap2 := as.numeric(gsub("^.*,", "", ADs))]
  mfr <- mfr[AD_hap2 != 0]
  mfr[, AD_hap1 := as.numeric(gsub(",.*$", "", ADs))]
  het <- mfr[ AD_hap1 != 0 & AD_hap2 != 0]
  return(het)
}


# Read phenotype
readPhenotypes <- function(vcf_file){
  geno <- fread(vcf_file, skip = 40) #fread gets messed up if it sees lines without strings b4 lines with strings (header is actually line 54)
  geno <- geno[nchar(REF) == 1 & nchar(ALT) == 1] #Remove the indels
  geno <- geno[FILTER == "PASS"]

  mvars <- names(geno)[10:ncol(geno)]
  geno <- geno[, c("#CHROM", "POS", "REF", "ALT", mvars), with = FALSE]
  mg <- melt(geno, measure.vars = mvars)
  mg <- melt(geno, id.vars = c("#CHROM", "POS", "REF", "ALT"))
  mg <- mg[!value %in% c("0|0", "1|1")]
  mg <- mg[!value %in% c("0", "1")] # Another way of encoding homozygous SNPs
}
