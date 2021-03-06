% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betabin.R
\name{betaBinomXI}
\alias{betaBinomXI}
\title{Fit mixture model}
\usage{
betaBinomXI(
  genic_dt,
  model = "AUTO",
  plot = FALSE,
  hist = FALSE,
  flag = 0,
  xciGenes = NULL,
  a0 = NULL,
  optimizer = c("nlminb", "optim"),
  method = NULL,
  limits = TRUE,
  keep_params = FALSE,
  debug = FALSE
)
}
\arguments{
\item{genic_dt}{A \code{data.table}. The table as outputted by \code{getGenicDP}.}

\item{model}{A \code{character} indicating which model to use to estimate
the mosaicism. Valid choices are "AUTO", "M0", "M1", "M2", "MF". See details.}

\item{plot}{A \code{logical}. If set to TRUE, information about the training
set and the skewing estimate will be plotted.}

\item{hist}{A \code{logical}. If set to TRUE, an histogram of the skewing
estimates will be displayed.}

\item{flag}{A \code{numeric}. Specify how to handle convergence issues. See
details.}

\item{xciGenes}{A \code{character} or NULL. To be passed to \code{readXCI} to
select the training set of inactivated genes.}

\item{a0}{A \code{numeric} or NULL. Starting values for the optimization. This
should not be used with more than one model as different models have
different parameters. Leave NULL unless you know what you're doing.}

\item{optimizer}{A \code{character}. The optimization function to use for minimization
of the log-likelihood. Should be one of "nlminb" or "optim".}

\item{method}{A \code{character}. The method to be passed to \code{optim}
when it is the selected \code{optimizer}.}

\item{limits}{A \code{logical}. If set to TRUE, the optimization will be
constrained. Using upper bounds on the probability of sequencing error and
escape in the training set ensures that the dominant mixture represents the
skewing for inactivated genes.}

\item{keep_params}{A \code{logical}. If set to TRUE, all parameters will be
reported instead of just the alpha and beta estimates. Can useful for
model specific analysis but clutters the table.}

\item{debug}{A \code{logical}. If set to TRUE, information about each iteration
will be printed (Useful to identify problematic samples).}
}
\value{
A \code{data.table} with an entry per sample and per gene.
}
\description{
Fit a mixture model to estimate mosaicism and XCI-escape.
}
\details{
The model determines the number of components used in the mixture model. By
default, "AUTO" tries all combinations of mixtures and the best estimate is
kept using backward selection based on AIC.
M0 is a simple beta-binomial. M1 adds a binomial component to model the
sequencing errors. M2 jointly models the probability of misclasification
in the training set. MF include all 3 components.

Flags in the output reports issues in convergence. If \code{flag} is set to 0,
nothing is done. If set to 1, the model selection will avoid flagged models
(will favor parcimonious models).
If set to 2, calls for which the best selected model had convergence issue
will be removed.
}
\examples{
library(data.table)
# Simulated data
dtf <- system.file("extdata/data2_vignette.tsv", package = "XCIR")
dt <- fread(dtf)
xcigf <- system.file("extdata/xcig_vignette.txt", package = "XCIR")
xcig <- readLines(xcigf)
# Run all models on the data
all <- betaBinomXI(dt, xciGenes = xcig)
# Simple BetaBinomial model and show histogram of skewing ~~A~~~
bb <- betaBinomXI(dt, xciGenes = xcig, model = "M0", hist = TRUE)

# Plotting fits
stoshow <- paste0("sample", c(31, 33, 35, 40)) #interesting samples
plotQC(all[sample \%in\% stoshow], xcig = xcig)

# Summarizing results
# Sample information
samps <- sample_clean(all)
# Gene-level predictions
xcistates <- getXCIstate(all)
}
\seealso{
getGenicDP readXCI
}
