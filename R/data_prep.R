# Keep SNPs from a pileup
getRNASNPs <- function(pileup_vcf, out = NULL){
  if(is.null(out))
    stop("out should be a character")
  system(paste("vcftools --vcf", pileup_vcf, "--get-INFO AD")) #Get allelic depth
  fr <- fread("out.INFO")
  hetSNP <- fr[REF != "." & ALT != "." & nchar(REF) == 1 & nchar(ALT) == 1]
  hetSNP_file <- paste0(pileup_vcf, "_hetSNP.tsv")
  write.table(hetSNP[, list(CHROM, POS)], file = hetSNP_file,  row.names = F, col.names = F, sep = "\t", quote = F) #Create a position file
  message(paste("bi-allelic positions were written in"), hetSNP_file)
  system(paste("vcftools --vcf", pileup_vcf, "--positions", hetSNP_file, "--recode --out", out)) #Keep only the positions selected above
  return(NULL)
}
