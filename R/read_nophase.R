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
#' @importFrom S4Vectors unstrsplit
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
  mdt[, POS := as.numeric(POS)]
  #mdt[, `:=`(c("AD_hap1", "AD_hap2"), tstrsplit(AD, split = ","))] #It doesn't actually matter which is REF and ALT
  #mdt[, AD_hap1 := as.numeric(AD_hap1)]
  #mdt[, AD_hap2 := as.numeric(AD_hap2)]
  mdt[, AD_hap1 := as.numeric(sapply(AD, "[[", 1))]
  mdt[, AD_hap2 := as.numeric(sapply(AD, "[[", 2))]
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

#' @export
readVCF4 <- function(vcf_file, merge_alts = F){
  vcf_param <- ScanVcfParam(fixed = c("ALT"), info = NA, geno = c("AD"))
  vcf <- readVcf(vcf_file, genome = "chrX", param = vcf_param)
  alts <- as.character(unstrsplit(alt(vcf)))
  refs <- as.character(ref(vcf))
  # Remove monomorphic and indels
  keep <- which(nchar(gsub("<\\*>", "", alts)) == 1 & nchar(refs) == 1)
  nomono <- vcf[keep]

  # Make a table with genotypes and counts
  geno <- geno(nomono)
  AD <- geno$AD
  sampleNames <-  gsub("\\..*", "", basename(colnames(AD)))
  rn <- rownames(AD)
  chr <- gsub(":.*", "", rn)
  pos <- as.numeric(gsub("_.*", "", gsub("^.*:", "", rn)))
  ret <- data.table(CHROM = chr, POS = pos, REF = refs[keep], ALT = alts[keep])

  # Check for alternate alleles formats
  alt2 <- unique(substr(alts[keep], 2, 4))
  alt2 <- alt2[alt2 != ""]
  if(length(alt2) > 0){ #Are there multiple alternate
    if(alt2 != "<*>"){  #If so, is it just structural variants
      err <- "There are more than one secondary alternate allele or it is not a deletion '<*>'"
      stop(err)
    }
  }
  #if(!(length(alt2) == 1 & alt2 == "<*>")){
  #  err <- "There are more than one secondary alternate allele or it is not a deletion '<*>'"
  #  stop(err)
  #}

  nref <- sapply(AD, "[[", 1)
  nalt1 <- sapply(AD, "[[", 2)

  # REF
  refDT <- data.table(matrix(nref, ncol = ncol(AD)))
  setnames(refDT, sampleNames)
  refDT[, `:=`(c("CHROM", "POS"), list(chr, pos))]
  mref <- melt(refDT, id.vars = c("CHROM", "POS"), value.name = "AD_hap1", variable.name = "sample")

  # ALT (Consider only the first ALT)
  altDT <- data.table(matrix(nalt1, ncol = ncol(AD)))
  setnames(altDT, sampleNames)
  altDT[, `:=`(c("CHROM", "POS"), list(chr, pos))]
  malt <- melt(altDT, id.vars = c("CHROM", "POS"), value.name = "AD_hap2", variable.name = "sample")

  counts <- merge(mref, malt, by = c("CHROM", "POS", "sample"))
  ret <- merge(ret, counts, by = c("CHROM", "POS"))[order(CHROM, POS)]
  return(ret)
}


# New functions for samtools v4.3
.makeDTFromGeno <- function(geno){
  sampleNames <-  gsub("\\..*", "", basename(colnames(geno)))
  rn <- rownames(geno)
  chr <- gsub(":.*", "", rn)
  pos <- gsub("_.*", "", gsub("^.*:", "", rn))
  dt <- data.table(geno)
  setnames(dt, sampleNames)
  dt[, CHROM := chr]
  dt[, POS := pos]
  return(dt)
}

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
