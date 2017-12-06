#The file should be a vcf file with
#Chr, pos, ref, AD_hap1, AD_hap2

#' Read SNPs from RNA-Seq
#'
#' Read SNPs from RNA-Seq that have not been phased.
#'
#' @param vcf_file A \code{character}. The path to a vcf file.
#'
#' @details
#' For phased samples, use \code{readXVcf}.
#'
#' @export
readRNASNPs <- function(vcf_file){
  vcf_param <- ScanVcfParam(fixed = c("ALT"), info = NA, geno = c("AD"))
  vcf <- readVcf(vcf_file, genome = "chrX", param = vcf_param)
  x = CharacterList(alt(vcf))
  x = unstrsplit(x)
  keep = which(nchar(x)==1) # Removes empty (monomorphic) and multiple (indels).
  nomono <- vcf[keep]

  # Extract allelic depth
  AD <- geno(nomono)$AD #We use unfiltered allele depth b/c filters are usually made for DNA data
  # AD is n1,n2, where n1 is the number of copies of ref and n2 is the number of copies of alt
  dt <- .makeDTFromGeno(AD)
  dt[, ID := keep] #To filter the vcf when dt is filtered
  ref <- as.character(ref(nomono))
  alt <- unstrsplit(CharacterList(alt(nomono)))
  dt[, REF := ref]
  dt[, ALT := alt]
  dt <- dt[nchar(REF) == 1] #indels are a lot harder to analyse and will move the reading frame. Remove them.
  nomono <- vcf[dt$ID] # At this point it's vcf - monomorphic sites - indels on ALT - indels on REF
  samples <- colnames(dt)[!colnames(dt) %in% c("POS", "CHROM", "ID", "REF", "ALT") ]
  mdt <- melt(dt, measure.vars = samples, variable.name = "sample", value.name = "AD")
  mdt[, ID := NULL]
  mdt[, `:=`(c("AD_hap1", "AD_hap2"), tstrsplit(AD, split = ","))] #It doesn't actually matter which is REF and ALT
  mdt[, AD_hap1 := as.numeric(AD_hap1)]
  mdt[, AD_hap2 := as.numeric(AD_hap2)]
  return(mdt)
}


# Read phenotype
readPhenotypes <- function(vcf_file){
  geno <- fread(vcf_file, skip = 40) #fread gets messed up if it sees lines without strings b4 lines with strings (header is actually line 54)
  geno <- geno[nchar(REF) == 1 & nchar(ALT) == 1] #Remove the indels
  geno <- geno[FILTER == "PASS"]

  mvars <- names(geno)[10:ncol(geno)]
  geno <- geno[, c("#CHROM", "POS", "REF", "ALT", mvars), with = F]
  mg <- melt(geno, measure.vars = mvars)
  mg <- melt(geno, id.vars = c("#CHROM", "POS", "REF", "ALT"))
  mg <- mg[!value %in% c("0|0", "1|1")]
  mg <- mg[!value %in% c("0", "1")] # Another way of encoding homozygous SNPs
}



