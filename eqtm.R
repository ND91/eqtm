#' eqtm
#'
#' This function was setup to perform an expression quantitative trait methylation analysis for methylation and gene expression data obtained from the same samples.
#' @param dmrs_gr A GenomicRanges-derived object where the ranges represent the coordinates of the DMR and whose metadata contains an associated gene identifier.
#' @param gene_col A column index containing of the dmrs_gr object containing the gene identifiers.
#' @param meth_data A matrix whose row names contain CpG identifiers and whose values represent the methylation signal. The column names must contain the columns names of the gene expression data.
#' @param expr_data A matrix whose row names contain gene identifiers and whose values represent the expression signal. The column names must contain the columns names of the methylation data.
#' @param meth_anno_gr A GenomicRanges-derived object where the ranges represent the coordinates of the assayed CpGs and whose rownames contain an identifier.
#' @param cor_type Type of correlation ("Pearson", or "Spearman").
#' @param N The number of bootstraps and permutations for calculating the 95% confidence intervals and p-values respectively (defaults to 1000).
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return 
#' @keywords eqtm, methylation, expression
#' @export
#' @import GenomicRanges
#' @import boot
#' @examples 

eqtm <- function(dmrs_gr, gene_col, meth_data, expr_data, meth_anno_gr, cor_type = "pearson", N = 1000, seed = NULL, ...){
  #Argument checking
  if(is.null(dmrs_gr)) stop("dmr_gr cannot be found")
  if(class(dmrs_gr) != "GRanges") stop("dmr_gr must be a GRanges object")
  if(is.null(gene_col)) stop("gene_col cannot be found")
  if(is.null(meth_data)) stop("meth_data cannot be found")
  if(is.null(expr_data)) stop("expr_data cannot be found")
  if(is.null(meth_anno_gr)) stop("meth_anno_gr cannot be found")
  if(class(meth_anno_gr) != "GRanges") stop("meth_anno_gr must be a GRanges object")
  
  ol_samples <- intersect(colnames(meth_data), colnames(expr_data))
  if(length(ol_samples) == 0) stop("There exists no overlap between the column names present in meth_data and expr_data")
  
  #Subsetting the datasets according to the genes
  dmr_genes <- data.frame(coordinates = gsub("[;-]", "_", as.character(dmr_genes)),
                          geneid = as.character(mcols(dmrs_gr)[, gene_col]))
  ol_genes <- intersect(dmr_genes$geneid, rownames(expr_data))
  
  dmrs_ol <- dmrs_gr[as.character(mcols(dmrs_gr)[, gene_col]) %in% ol_genes,]
  expr_ol <- expr_data[rownames(expr_data) %in% ol_genes, ol_samples]
  meth_anno_ol <- meth_anno_gr[unique(subjectHits(findOverlaps(dmrs_ol, meth_anno_gr))),]
  meth_ol <- meth_data[names(meth_anno_ol), ol_samples]
  
  #No need to regenerate the null distribution for same length DMRs
  dmrs_split <- split(dmrs_ol, seqnames(dmrs_ol))
  
  ncpgs_chrom <- lapply(dmrs_split, function(chrom){
    sort(as.vector(unique(table(queryHits(findOverlaps(chrom, meth_anno_ol))))))
  })
  
  #Methylation aggregation
  dmrs_meth <- dmr_aggregator(dmrs_gr = dmrs_ol, meth_data = meth_ol, meth_anno = meth_anno_ol)
  
  #Correlation
  lapply(dmrs_agg)
  
  if(!is.null(seed)) set.seed(seed)
  correlations <- sapply(X = overlapping_dmr_gene_anno, simplify = T, USE.NAMES = F, FUN = function(dmr_gene_entry){
    cor_df <- data.frame(meth = as.vector(unlist(mcols(subsetByOverlaps(dmr_data_culled, dmr_gene_entry)))),
                         expr = as.vector(unlist(expr_data_culled[unlist(mcols(dmr_gene_entry[, id_col_name])), ])))
    ci_bootstrapper(cor_df = cor_df, iterations = iterations, cor_type = cor_type)
  })
}

dmr_aggregator <- function(dmrs_gr, meth_data, meth_anno){
  meth_fo <- data.frame(findOverlaps(dmrs_gr, meth_anno))
  meth_fo_split <- split(meth_fo, meth_fo$queryHits)
  dmr_mean <- lapply(meth_fo_split, FUN = function(dmr){
    colMeans(meth_data[dmr$subjectHits,], na.rm = T)
  })
  dmr_mean <- do.call(rbind, dmr_mean)
  rownames(dmr_mean) <- gsub("[:-]", "_", as.character(dmrs_gr))
  
  return(dmr_mean)
}


ci_bootstrapper <- function(cor_df, iterations, cor_type){
  #Bootstraps for the CI
  me_correlator_correlator <- function(df, indices){
    return(cor(df[indices, "meth"], df[indices, "expr"], method = cor_type)) 
  }
  boot_results <- boot(cor_df, me_correlator_correlator, R = iterations, stype = "i")
  CI95 <- boot.ci(boot_results, type = "bca")
  
  return(c(corcoef = boot_results$t0, corcoefCI95_lower = CI95$bca[4], corcoefCI95_upper = CI95$bca[5]))
}