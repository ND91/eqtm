#' eqtm
#'
#' This function was setup to perform an expression quantitative trait methylation analysis for methylation and gene expression data obtained from the same samples.
#' @param dmrs_gr A GenomicRanges-derived object where the ranges represent the coordinates of the DMR and whose metadata contains an associated gene identifier.
#' @param gene_col A column index containing of the dmrs_gr object containing the gene identifiers.
#' @param meth_data A matrix whose row names contain CpG identifiers and whose values represent the methylation signal. The column names must contain the columns names of the gene expression data.
#' @param expr_data A matrix whose row names contain gene identifiers and whose values represent the expression signal. The column names must contain the columns names of the methylation data.
#' @param meth_anno A GenomicRanges-derived object where the ranges represent the DMR and whose metadata contains an associated gene identifier.
#' @param cor_type Type of correlation ("Pearson", or "Spearman").
#' @param N The number of bootstraps and permutations to be performed when calculating the 95% confidence intervals and p-values respectively.
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return 
#' @keywords eqtm, methylation, expression
#' @export
#' @import bsseq
#' @import GenomicRanges
#' @import boot
#' @examples 
#' 

eqtm <- function(dmrs_gr, meth_data, expr_data, meth_anno){
  #Split the DMRs according to chromosome
  dmrs_split <- split(dmrs_gr, seqnames(dmrs_gr))
  
  #Generate a null-distribution for each number of CpGs per chromosome. No need to recalculate the null-distribution for
  ncpgs_chrom <- lapply(dmrs_split, function(chrom){
    as.vector(unique(table(queryHits(findOverlaps(chrom, meth_anno)))))
  })
  
  dmr_mean(dmr_gene = dmr_gene, methylation = methylation)
  
  
  anno_cpgs <- subsetByOverlaps(meth_anno, dmrs_gr[1,])
  names(anno_cpgs)
}