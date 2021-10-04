
snpmart = useEnsembl(biomart = "snp", dataset="hsapiens_snp")
getBM(attributes = c('refsnp_id','allele','chrom_start','chrom_strand'), 
      filters = c('chr_name','start','end'), 
      values = list(8, 148350, 148420), 
      mart = snpmart)


ggplot(aanno, aes(x = as.factor))


nt <- rnbinom(40, size = 80, prob = 0.5)
y <- rbinom(n = 40, size = nt, prob = 0.5)
dt <- data.table(nt, y)
ggplot(dt, aes(x = nt, y = y)) + geom_point()

glm(y ~ nt)
glm(y ~ nt, family = )

#########################################
# With annotation packages
BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")

library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
chrX_snps <- snpsBySeqname(snps, "X")
dtX <- as.data.table(chrX_snps)

# Merging
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
atts <- c("hgnc_symbol", "chromosome_name", "start_position", 
          "end_position", "gene_biotype")
p2g <- data.table(getBM(attributes = atts, mart = mart, 
                        filters = "chromosome_name", values = "X"))
gs <- p2g[gene_biotype %in% c("protein_coding", "pseudogene", "lncRNA") & hgnc_symbol != ""]
gs[duplicated(hgnc_symbol)]
#########################################
snpmart <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
snpatts <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
snpfilt <- 
res <- getBM(attributes = listAttributes(snpmart)[1:8, 1], mart = snpmart,
             filters = "snp_filter", values = "rs1296676918")
res <- getBM(attributes = listAttributes(snpmart)[1:8, 1], mart = snpmart,
             filters = "start", values = "2781481")
