#' @title Estimating inactivated X chromosome expression
#'
#' @description Tools for the analysis of X chromosome inactivation (XCI) and XCI-escape inference.
#'
#' @name XCIR-package
#' @aliases XCIR
#' @author Renan Sauteraud <rxs575@psu.edu>
#' @import data.table methods stats
#' @importFrom utils capture.output write.table globalVariables
NULL

gv <- c("AD_alt", "AD_hap1", "AD_hap2", "AD_ref", "AIC_bb", "AIC_mm",
        "AIC_mm2", "AIC_mm3", "ALT", "ANNO", "ANNO_FULL", "CHROM", "FILTER",
        "GENE", "ID", "N", "N.x", "N.y", "N_escape", "N_info", "N_support",
        "Ninac", "Nsamples", "Ntrain", "POS", "REF", "Tobs", "V1", "V2",
        "a", "a1", "a2", "a_est", "a_hap1", "a_hap2", "ai", "b_est",
        "capture.output", "count", "den", "denom", "dp1", "escape", "f",
        "facet_wrap", "fg", "gNrm", "sex", "geom_hline", "geom_point",
        "geom_rect", "h1", "h2", "index", "ivw_tau", "jcss", "k", "label",
        "logL", "m1", "m2", "model", "n1", "n2", "n_snps", "nsamples",
        "num", "p_het", "p_inac", "p_value", "pe", "percent", "perm_pvalue",
        "permout", "pi_err", "pi_escape", "status", "statusJ", "tau",
        "theme", "tot", "type", "value", "var_fg", "var_tau", "write.table",
        "xiexpr", "ymax", "ymin", "convergence", "flag", "gene_biotype",
        "GENE_pos", "Ntot", "XCIstate", "i.start", "i.end", "XCIstate_skew",
        "pbb", ".")
globalVariables(gv)

