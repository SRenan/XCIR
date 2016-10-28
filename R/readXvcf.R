makeDTFromGeno <- function(geno){
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
#' @importFrom VariantAnnotation ScanVcfParam readVcf geno alt ref
#' @importFrom IRanges CharacterList
#' @importFrom Biostrings unstrsplit
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
  dt <- makeDTFromGeno(AD)
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
    setnames(dtGT, c("CHROM", "ID", "POS", "a1", "a2", "h1", "h2"))
    dtGT[, CHROM := as.character(CHROM)]
    dtGT[CHROM == 23, CHROM := "X"]
    bys <- c("CHROM", "POS")
  } else{
    GT <- geno(nomono)$GT
    dtGT <- makeDTFromGeno(GT)
    dtGT <- melt(dtGT, measure.vars = samples, variable.name = "sample", value.name = "GT")
    #TODO: split the genotype columns into h1, h2
    bys <- c("sample", "CHROM", "POS")
  }

  dt2 <- merge(mdt, dtGT, by = bys)
  #if(rm_homo){
  #  dt2 <- dt2[!GT %in% c("0/0", "1/1")]
  #} else{
  #  dt2[GT %in% c("0/0", "1/1"), AD := "NA,NA"]
  #}
  dt2[, `:=`(c("AD_ref", "AD_alt"), tstrsplit(AD, split = ","))]
  #dt2[GT == "1/0",`:=`(c("AD_ref", "AD_alt"), list(AD_alt, AD_ref))]
  #dt2[, AD_ref := as.numeric(AD_ref)] #will coerce "NA" into NA with a warning
  #dt2[, AD_alt := as.numeric(AD_alt)]

  if(rm_homo){
    dt2 <- dt2[!(h1 == 0 & h2 == 0) & !(h1 == 1 & h2 == 1)]
  } else{
    dt2[!(h1 == 0 & h2 == 0) & !(h1 == 1 & h2 == 1), `:=`(c("h1", "h2"), c(0,0))]
  }
  dt2[, AD_ref := as.numeric(AD_ref)]
  dt2[, AD_alt := as.numeric(AD_alt)]
  # Si h1 == 1. L'allele sur chrA est a1. si a1 == ref. AD_chrA <- AD_ref
  dt2[, a_chrA := ifelse(h1 == 0, a1, a2)]
  dt2[, a_chrB := ifelse(h1 == 0, a2, a1)]
  dt2[, AD_chrA := ifelse(a_chrA == REF, AD_ref, AD_alt)]
  dt2[, AD_chrB := ifelse(a_chrB == REF, AD_ref, AD_alt)]
  return(dt2)
}


#' @importFrom seqminer annotatePlain makeAnnotationParameter
seqm_anno <- function(input, output_dir, reference = NULL, geneFile = NULL){
  if(is.null(reference)){
    reference <- system.file("extdata", "human.g1k.v37.fa", package = "XCIR")
    if(!file.exists(reference))
      stop("No reference file given or found.")
  }
  if(is.null(geneFile)){
    geneFile <- system.file("extdata", "refFlat.gencode.v19.gz", package = "XCIR")
    if(!file.exists(geneFile))
      stop("No gene file given or found.")
  }
  inFile <- tempfile(tmpdir = output_dir)
  write_siteRef(input, ref_file = inFile) #Write annotations into the file that will be used by seqminer
  param <- (list(reference = reference, geneFile = geneFile, inputFormat = "plain"))
  output_file <- file.path(output_dir, "anno_out")
  annotatePlain(inFile, output_file, params=param) #seqminer does not allow using an existing filename
  return(output_file)
}

#IN: the return of readXVcf, anno_file generated by dajiang's scripts. test_siteRef.tsv.anno
#OUT: dt merged with the annotations //noMono.chr*.DPR.anno//
#' Read annotation file
#'
#' Read a given annotation file and merge it with a data.table containing the
#' relevant information to exstimate inactivated X chromosome expression.
#'
#' @param dt A \code{data.table} object.
#' @param seqm_annotate A \code{logical}. If set to TRUE, the \code{seqminer}
#' package will be used to annotate \code{dt}
#' @param anno_file A \code{character}. The name of a file containing annotations.
#'
#' @return A \code{data.table} object that contains allelic coverage, genotype
#' and annotations at the covered SNPs.
#'
#' @export
addAnno <- function(dt, seqm_annotate = TRUE, anno_file = NULL){
  if(seqm_annotate){
    output_dir <- tempdir()
    anno_file <- seqm_anno(input = copy(dt), output_dir)
  } else if(is.null(anno_file)){
    stop("No annotations given, set seqm_annotate to TRUE or provide an annotation filename to anno_file")
  }
  anno <- fread(anno_file)
  setnames(anno, gsub("#", "", names(anno)))
  anno[, GENE := gsub("^.*:", "", ANNO)]
  anno[, ANNO := gsub(":.*$", "", ANNO)]
  anno[, CHROM := as.character(CHROM)]
  anno[, POS := as.numeric(POS)]
  anno[, ANNO_FULL := NULL]
  anno <- anno[GENE != "Intergenic"]
  DPR.anno <- merge(unique(anno), dt,  by = c("CHROM", "POS", "REF", "ALT"))
  if(seqm_annotate){
    unlink(output_dir, recursive = TRUE)
  }
  return(DPR.anno)
}

#' @export
write_siteRef <- function(dt, ref_file = "test_ref.anno.tsv"){
  write.table(unique(dt[, list(CHROM, POS, REF, ALT)]), file = ref_file,
              row.names = FALSE, quote = FALSE, sep = "\t")
}
