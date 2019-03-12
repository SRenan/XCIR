# This is an alternative to using seqminer for annotations
# Instead of relying on a large annotation file, and seqminer,
# we fetch gene coordinates from bioMart and find whether the given list
# of positions are covered by any of them

#' biomaRt genes
#'
#' Extract gene informations from biomaRt
#'
#' @param release A \code{character}. Genome release name. Valid releases are
#'  "hg19", "hg38".
#' @param chr A \code{character} or \code{NULL}. If specified, only the genes
#'  from the specified chromosomes will be returned.
#'
#' @return A \code{data.table} with the gene symbol, start and end position
#'  and matching ensembl transcripts.
#'
#' @examples
#' #Chromosome X, hg19
#' egX <- mart_genes()
#' #Full genome, latest release
#' eg <- mart_genes("hg38")
#'
#' @importFrom biomaRt useMart getBM
#' @export
mart_genes <- function(release = "hg19", chr = "X"){
  release <- tolower(release)
  if(release %in% c("grch38", "hg38", "latest")){
    mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    message("Fetching genes for hg38")
  } else if(release %in% c("grch37", "hg19")){
    mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    message("Fetching genes for hg19")
  } else if(release %in% c("hg17", "hg18")){
    stop(paste(release, "is too old (hg18 = 2006)"))
  } else{
    stop("Incorrect release number")
  }
  atts <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "gene_biotype")
  if(!is.null(chr)){
    p2g <- data.table(getBM(attributes = atts, mart = mart, filters = "chromosome_name", values = "X"))
  } else{
    p2g <- data.table(getBM(attributes = atts, mart = mart))
  }
  setnames(p2g, c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), c("GENE", "CHROM", "start", "end"))
  return(p2g)
}

#' @export
annotateX <- function(xciObj, read_count_cutoff = 20, het_cutoff = 3,
                      release = "hg19", verbose = F){
  # 1) Extract the genes
  egX <- mart_genes(release)[GENE != ""]
  egX <- egX[gene_biotype %in% c("protein_coding", "pseudogene", "lincRNA")] #Remove miRNA, rRNA,
  # Handle genes that have multiple start/end positions
  dupg <- egX[duplicated(egX$GENE), GENE]
  if(verbose){
    message(paste("There are", length(dupg), "genes with multiple entries."))
  }
  egX[, start:= min(start), by = "GENE"]
  egX[, end:= max(end), by = "GENE"]
  egX <- unique(egX)

  #2) Map positions to the genes
  # Filtering first to make everything faster
  xciAnno <- xciObj[AD_hap1 + AD_hap2 >= read_count_cutoff & AD_hap1 >= het_cutoff & AD_hap2 >= het_cutoff]
  xciAnno[, ID := paste(CHROM, POS, sep = ":")]
  p2gX <- pos2gene(xciAnno[, list(ID, POS)], egX)
  p2gX[, GENE_pos := (start+end)/2] #This can be useful for plots to order genes
  xciAnno <-  merge(p2gX[, list(ID, GENE, GENE_pos)], xciAnno, by = "ID") #The merge naturally removes intergenic SNPs
  return(xciAnno)
}

pos2gene <- function(target, interval){
  pos <- unique(copy(target)) #Still need to copy b/c if target is already unique we'll overwrite
  gene <- copy(interval)
  # TODO:
  # a) Handle chromosome
  # b) Fetch genes from biomaRt

  # INPUTS
  # Assume first column is an ID and second is the positions
  setnames(pos, c("ID", "POS"))
  pos[, start := POS]
  pos[, end := POS]
  setkey(pos, start, end)
  # Assume genes have columns GENE, start, end
  olps <- foverlaps(gene, pos, type = "any")
  ret <- olps[!is.na(ID), list(ID, POS, GENE, i.start, i.end)]
  setnames(ret, c("i.start", "i.end"), c("start", "end"))
  return(ret)
}

################################################################################
# Annotation using seqminer
################################################################################

#' @importFrom seqminer annotatePlain makeAnnotationParameter
.seqm_anno <- function(input, output_dir, reference = NULL, geneFile = NULL){
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
  .write_siteRef(input, ref_file = inFile) #Write annotations into the file that will be used by seqminer
  param <- (list(reference = reference, geneFile = geneFile, inputFormat = "plain"))
  output_file <- file.path(output_dir, "anno_out")
  capture.output({
    annotatePlain(inFile, output_file, params=param) #seqminer does not allow using an existing filename
  })
  return(output_file)
}


#' Read annotation file
#'
#' Read a given annotation file and merge it with a data.table containing the
#' relevant information to estimate inactivated X chromosome expression and
#' filter out SNPs with low coverage.
#'
#' @param dt A \code{data.table} object.
#' @param seqm_annotate A \code{logical}. If set to TRUE, the \code{seqminer}
#' package will be used to annotate \code{dt}. If set to FALSE, this function
#' is a simple read count filtering step.
#' @param read_count_cutoff A \code{numeric}. Keep only SNPs that have at least
#'  that many reads.
#' @param het_cutoff A \code{numeric}. Keep only SNPs that have at least that
#'  many reads on each allele.
#' @param filter_pool_cutoff A \code{numeric}. Keep only SNPs that have at
#'  least that many reads on each allele across all samples. See details for
#'  more information.
#' @param anno_file A \code{character}. The name of a file containing annotations.
#'
#' @details
#' If the samples all have the same genotype (e.g: technical replicates),
#' \code{filter_pool_cutoff} will sum counts across samples and preserve SNPs
#' that pass the cutoff on both the reference and alternate alleles. This
#' may lead to samples with 0 counts on either allele but will prevent removing
#' heterozygous sites with lower coverage (especialliy in skewed samples).
#'
#' @return A \code{data.table} object that contains allelic coverage, genotype
#' and annotations at the covered SNPs.
#'
#' @export
addAnno <- function(dt, seqm_annotate = TRUE, read_count_cutoff = 20,
                    het_cutoff = 3, filter_pool_cutoff = 3, anno_file = NULL){
  dt <- dt[AD_hap1 + AD_hap2 > read_count_cutoff]
  if(seqm_annotate){
    output_dir <- tempdir()
    anno_file <- .seqm_anno(input = copy(dt), output_dir)
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
  DPR_anno <- merge(unique(anno), dt,  by = c("CHROM", "POS", "REF", "ALT"))

  # Finding heterozygous SNPs
  if(filter_pool_cutoff > 0){
    poskp <- DPR_anno[, list(sum(AD_hap1), sum(AD_hap2)), by = c("CHROM", "POS")][V1 >= filter_pool_cutoff & V2 >= filter_pool_cutoff][, POS]
    DPR_anno <- DPR_anno[POS %in% poskp]
  }
  DPR_anno <- DPR_anno[AD_hap1 >= het_cutoff & AD_hap2 >= het_cutoff]

  return(DPR_anno)
}

#' @export
.write_siteRef <- function(dt, ref_file = "test_ref.anno.tsv"){
  write.table(unique(dt[, list(CHROM, POS, REF, ALT)]), file = ref_file,
              row.names = FALSE, quote = FALSE, sep = "\t")
}

#' XCI consensus
#'
#' Read consensus & XCIR calls for all X-linked genes
#'
#' @param redownload A \code{logical}. If set to TRUE, the original supplementary
#' file is redownloaded from PMC.
#' @param simple A \code{logical}. If set to TRUE, minimal information is returned,
#' only for genes with an available XCIR classification.
#'
#' @details
#' The consensus is as published in Supplementary table S1 of
#' Balaton et al. (Biol Sex Differ. 2015). doi: 10.1186/s13293-015-0053-7
#'
#' @importFrom readxl read_xlsx
#' @export
consensusXCI <- function(redownload = F, simple = T){
  if(redownload){
    tmp = tempfile()
    download.file("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696107/bin/13293_2015_53_MOESM1_ESM.xlsx", destfile = tmp, mode = "wb")
    cons <- data.table(read_xlsx(tmp))
  } else{
    cons <- data.table(read_xlsx(system.file("extdata/13293_2015_53_MOESM1_ESM.xlsx", package = "XCIR")))
  }
  setnames(cons, tolower(chartr(" ", "_", names(cons))))
  geu <- fread(system.file("extdata/GEU_xcir_181122_flag.tsv", package = "XCIR"))
  geuall <- getXCIstate(geu)
  geuskew <- getXCIstate(geu[f < .25])
  geus <- merge(geuall, geuskew, by = "GENE", suffixes = c("_all", "_skew"), all.x = T)
  cons <- merge(geus, cons, by.x = "GENE", by.y = "gene_name", all.y = T)
  cons <- cons[, c(names(cons)[1:8], "y_homology"), with = F]
  if(simple)
    cons <- cons[!is.na(XCIstate_skew)]
  return(cons)
}
