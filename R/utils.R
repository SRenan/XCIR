readXCI <- function(xciGenes = NULL){
  if(!is.null(xciGenes)){
    if(length(xciGenes) > 1){
      # If the input is a character vector, assume it's gene symbols
      xci <- xciGenes
      return(xciGenes)
    } else if(xciGenes == "cotton"){
      xci <- system.file("extdata", "xciGene_cotton.txt", package = "XCIR")
      xci <- readLines(xci)
    } else if(file.exists(xciGenes)){
      xci <- readLines(xciGenes)
    }
  } else{
    xci <- system.file("extdata", "xciGene.txt", package = "XCIR")
    xci <- readLines(xci)
  }
}

# Combine snps from pileup and array into a new vcf for phasing
#
# array should be GM07345.vcf
# pileup should be the result of ->
#array <- "~/workspace/XCI/snp1Kgp/GM07345.vcf" #6470 heterozyguous sites
#pileup <- "~/workspace/XCI/umap/het.site.RNAseq.txt" #401 new SNPs
#pileup <- "~/workspace/XCI/umap/het.site.RNAseq0.txt" #993 new SNPs


# This takes the output of a pileup and subset to get only the heterozyguous sites
getRNAhet <- function(pileup_vcf, thresh = 0, out = NULL){
  if(is.null(out))
    stop("out should be a character")
  system(paste("vcftools --vcf", pileup_vcf, "--get-INFO AD")) #Get allelic depth
  system(paste0("grep \",\" out.INFO | sed s/\",\"/\"\t\"/g | awk '{if(NF==6 && length($3)==1 && length($4)==1 && $5>", thresh, " && $6>", thresh, ") print $0}' > ", out)) #keep only heterozyguous SNPs that have at least thresh reads on each allele
  message(paste("The file", out, " can be used as input of combineVCF along with genotyped SNPs"))
  return(NULL)
}

#
combineVCF <- function(pileup, array, out){
  vcf_dt <- fread(array, skip = 29)
  sn <- colnames(vcf_dt)[ncol(vcf_dt)]
  setnames(vcf_dt, c("#CHROM", sn), c("CHROM", "sample"))
  vcf_dt <- vcf_dt[CHROM == 23]
  vcf_dt[, CHROM := as.character(CHROM)]
  vcf_dt[, CHROM := "X"]
  vcf_dt <- vcf_dt[! sample %in% c("1/1", "0/0", "./.")]

  rna_dt <- fread(pileup) #This is already the subset of heterozyguous SNPs w/ at least 2 SNPs on each allele
  setnames(rna_dt, c("CHROM", "POS", "REF", "ALT", "AD1", "AD2"))
  rna_dt[, ID := "."]
  rna_dt[, QUAL := "."]
  rna_dt[, FILTER := "."]
  rna_dt[, INFO := "."]
  rna_dt[, FORMAT := "GT"]
  rna_dt[, sample := "0/1"]

  new_vcf <- merge(vcf_dt, rna_dt, all = TRUE)
  new_vcf <- new_vcf[!(POS %in% new_vcf[duplicated(new_vcf$POS)]$POS & is.na(AD1))]
  new_vcf[, `:=`(c("AD1", "AD2"), list(NULL, NULL))]
  setnames(new_vcf, "CHROM", "#CHROM")

  msg <- paste0(nrow(new_vcf), " rows. ", length(which(rna_dt$POS %in% vcf_dt$POS)), "rows existed in both datasets. Kept the ones from RNA-Seq")
  message(msg)
  write.table(new_vcf, file = out,
              row.names = FALSE, quote = FALSE, sep = "\t")
  return(new_vcf)
}
