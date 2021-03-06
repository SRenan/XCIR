.makeDTFromGeno <- function(geno){
  # Fix samplenames, find CHROM and POS from rownames.
  sampleNames <-  gsub("\\..*", "", basename(colnames(geno)))
  rn <- rownames(geno)
  chr <- gsub(":.*", "", rn)
  pos <- gsub("_.*", "", gsub("^.*:", "", rn))
  if(is.list(geno[1,1])){
    ul <- unlist(geno)
    vals <- paste(ul[seq(1, length(ul), by = 2)], ul[seq(2, length(ul), by = 2)], sep = ",")
    geno <- data.table(matrix(vals, ncol = ncol(geno)))
  }
  dt <- data.table(geno)
  setnames(dt, sampleNames)
  dt[, POS := as.numeric(pos)]
  dt[, CHROM := as.character(chr)]
  return(dt)
}

#' Read a vcf file to extract information relevant to XCI
#'
#' Read a given vcf file and extract minimum information required for the
#' estimation of the expression of inactivated X chromosome.
#'
#' @param vcf_file A \code{character}. The path to a vcf file.
#' @param haps_file A \code{character}. The path to a .haps file, the output of
#' the shapeit software.
#' @param rm_homo A \code{logical}. If set to TRUE, homozygous sites will be
#' removed from the output. Otherwise, they will be set to NA.
#'
#' @return A \code{data.table} object in long format of length sample*site.
#'
#' @references
#' O. Delaneau, B. Howie, A. Cox, J-F. Zagury, J. Marchini (2013) Haplotype
#' estimation using sequence reads. American Journal of Human Genetics 93 (4)
#'  787-696.
#'
#' @example inst/examples/workflow.R
#'
#'
#' @importFrom VariantAnnotation ScanVcfParam readVcf geno alt ref
#' @importFrom IRanges CharacterList
#' @importFrom Biostrings unstrsplit
#' @importFrom tools file_ext
#' @export
readXVcf <- function(vcf_file, haps_file = NULL, rm_homo = TRUE){
  vcf_param <- ScanVcfParam(fixed = c("ALT"), info = NA, geno = c("GT", "AD"))
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

  # Add reference and alternate alleles
  ref <- as.character(ref(nomono))
  alt <- unstrsplit(CharacterList(alt(nomono)))
  dt[, REF := ref]
  dt[, ALT := alt]
  dt <- dt[nchar(REF) == 1] #indels are a lot harder to analyse and will move the reading frame. Remove them.
  nomono <- vcf[dt$ID] # At this point it's vcf - monomorphic sites - indels on ALT - indels on REF
  samples <- colnames(dt)[!colnames(dt) %in% c("POS", "CHROM", "ID", "REF", "ALT") ]
  mdt <- melt(dt, measure.vars = samples, variable.name = "sample", value.name = "AD")
  mdt[, ID := NULL]

  # Extract GT
  if(!is.null(haps_file)){
    message("Specifying haps_file only works when the data contains a single haplotype.")
    dtGT <- fread(haps_file)
    ext <- file_ext(haps_file)
    if(ext == "vcf"){
      # Using shapeit to get vcf file: shapeit -convert --input-haps gwas.phased --output-vcf gwas.phased.vcf
      # Where there are two files gwas.phased.sample and gwas.sample.haps
      # Using phaser to get vcf file: use the --gw_phase_vcf 1 option to output a vcf file with modified GT information
      # python2.7 ~/programs/phaser/phaser/phaser.py --vcf combined_phased.vcf --bam hisat2-LC1.sorted.bam --paired_end 1 --gw_phase_vcf 1 --pass_only 0 --mapq 1 --baseq 1 --sample GM07345 --o combined_phasered
      setnames(dtGT, "#CHROM", "CHROM")
      dtGT[, CHROM := as.character(CHROM)]
      dtGT[, GT := gsub(":.*$", "", get(names(dtGT)[10:ncol(dtGT)]))]
      unphased <- grep("/", dtGT$GT)
      if(length(unphased) > 0){
        warning(paste0(length(unphased), "/", nrow(dtGT), "SNPs are unphased."))
      }
      dtGT[, `:=`(c("h1", "h2"), tstrsplit(GT, split = "\\|"))]
      dtGT[, a_hap1 := ifelse(h1 == 0, REF, ALT)]
      dtGT[, a_hap2 := ifelse(h2 == 0, REF, ALT)]
    } else if(ext %in% c("haps", "")){ #make it work with empty for shiny
      setnames(dtGT, c("CHROM", "ID", "POS", "a1", "a2", "h1", "h2"))
      dtGT[, CHROM := as.character(CHROM)]
      dtGT[, a_hap1 := ifelse(h1 == 0, a1, a2)]
      dtGT[, a_hap2 := ifelse(h1 == 0, a2, a1)]
    } else{
      stop("Unknown haps_file format.")
    }
    dtGT[CHROM == 23 | tolower(CHROM) == "chrx", CHROM := "X"]
    bys <- c("CHROM", "POS")
  } else{
    message("No haps file specified. The vcf file must contain phased samples.")
    GT <- geno(nomono)$GT
    dtGT <- .makeDTFromGeno(GT)
    dtGT <- melt(dtGT, measure.vars = samples, variable.name = "sample", value.name = "GT")
    dtGT[, `:=`(c("h1", "h2"), tstrsplit(GT, split = "/"))]
    #TODO: split the genotype columns into h1, h2
    bys <- c("sample", "CHROM", "POS")
  }
  dt2 <- merge(mdt, dtGT, by = bys)

  dt2[, `:=`(c("AD_ref", "AD_alt"), tstrsplit(AD, split = ","))]
  dt2[, AD_ref := as.numeric(AD_ref)]
  dt2[, AD_alt := as.numeric(AD_alt)]

  if(rm_homo){
    dt2 <- dt2[!(h1 == 0 & h2 == 0) & !(h1 == 1 & h2 == 1)]
  } else{
    dt2[!(h1 == 0 & h2 == 0) & !(h1 == 1 & h2 == 1), `:=`(c("h1", "h2"), c(0,0))]
  }

  if("REF.x" %in% names(dt2)){
    setnames(dt2, c("REF.x", "ALT.x"), c("REF", "ALT"))
  }
  if(!"a_hap1" %in% colnames(dt2)){ #Already phased samples
    dt2[, a_hap1 := ifelse(h1 == 0, REF, ALT)]
    dt2[, a_hap2 := ifelse(h2 == 0, REF, ALT)]
  }
  # Si h1 == 1. L'allele sur chrA est a1. si a1 == ref. AD_chrA <- AD_hap1
  dt2[, AD_hap1 := ifelse(a_hap1 == REF, AD_ref, AD_alt)]
  dt2[, AD_hap2 := ifelse(a_hap2 == REF, AD_ref, AD_alt)]
  dt2 <- dt2[, list(CHROM, POS, REF, ALT, sample, h1, h2, a_hap1, a_hap2, AD_hap1, AD_hap2)]
  return(dt2)
}


