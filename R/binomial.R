aq <- copy(anno_q10)
aqSV <- copy(anno_q10)

# Start from the anno.
# fraction for the sample
# fraction for the gene
# variance for the sample fraction
# variance for the gene fraction

## Options
rm_inac <- 2 # If LC2 or LC3 have more than that many reads for a GENE, then the GENE is not considered inactivated


## Get Fraction for the sample
# AD1 ~ Bin(AD1 + AD2, frac)
# Which then gives us MLE for frac, then we can get its variance (AD1+AD2)*frac*(1-frac)
xci <- XCIR:::readXCI()
aq <- aq[GENE %in% xci]
aq[, frac_s := AD_hap1/(AD_hap1 + AD_hap2)]
aq[, frac_sample := mean((AD_hap1)/(AD_hap1+AD_hap2)), by = "sample"]
