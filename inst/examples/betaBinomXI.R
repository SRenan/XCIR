library(data.table)
# Simulated data
dtf <- system.file("extdata/data2_vignette.tsv", package = "XCIR")
dt <- fread(dtf)
xcigf <- system.file("extdata/xcig_vignette.txt", package = "XCIR")
xcig <- readLines(xcigf)
# Run all models on the data
all <- betaBinomXI(dt, xciGenes = xcig)
# Simple BetaBinomial model and show histogram of skewing
bb <- betaBinomXI(dt, xciGenes = xcig, model = "BB", hist = TRUE)

# Plotting fits
stoshow <- paste0("sample", c(1, 11, 31, 40)) #interesting samples
plotQC(all[sample %in% stoshow], xcig = xcig)

# Summarizing results
# Sample information
samps <- sample_clean(all)
# Gene-level predictions
xcistates <- getXCIstate(all)
