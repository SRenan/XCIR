# Example workflow for documentation

vcff <- system.file("extdata/AD_example.vcf", package = "XCIR")
# Reading functions
vcf <- readRNASNPs(vcff)
vcf <- readVCF4(vcff)

# Annotation functions
# Using seqminer (requires additional annotation files)
\donttest{
anno <- addAnno(vcf)
}
# Using biomaRt
anno <- annotateX(vcf)
# Do not remove SNPs with 0 count on minor allele
anno0 <- annotateX(vcf, het_cutoff = 0)

# Summarise read counts per gene
# Assuming data is phased, reads can be summed across genes.
genic <- getGenicDP(anno, highest_expr = FALSE)
# Unphased data, select SNP with highest overall expression.
genic <- getGenicDP(anno, highest_expr = TRUE)
