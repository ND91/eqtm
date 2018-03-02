#' dmr_mean 
#'
#' This function was setup to perform an expression quantitative trait methylation analysis for methylation and expression data obtained from the same samples.
#' @param dmr_gene A GenomicRanges or bsseq object where the coordinates represent the CpG locations and the metadata represents the methylation signal. Metadata column names must be the same as the expression data column names.
#' @param transcription A dataframe of transcripts with gene/transcript identifiers as rownames. Column names must be the same as the methylation metadata column names.
#' @param dmr_gene_annotation A GenomicRanges object where the coordinates represent the CpG locations and the metadata contains at least a column of the gene/transcript identifiers as used in the transcription matrix.
#' @param ID_column A string representing the name of the column in dmr_gene_annotation that contains the gene/transcript identifiers.
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

dmr_mean <- function(dmr_gene, methylation, ...){
  if(is.null(dmr_gene)) stop("No GRanges object with DMRs provided")
  if(is.null(methylation)) stop("No GRanges object with methylation data provided")
  
  dmr_chr <- gsub(pattern = "chr", replacement = "", x = as.character(seqnames(dmr_gene)), ignore.case = T)
  
  meth_chr <- as.character(seqnames(methylation))
  meth_pos <- start(methylation)
  
  if(class(methylation) == "BSseq"){
    betas <- bsseq::getMeth(methylation[which(meth_chr == dmr_chr & meth_pos >= start(dmr_gene) & meth_pos <= end(dmr_gene)), ])
  } else if(class(methylation) == "GenomicRanges"){
    betas <- methylation[which(meth_chr == dmr_chr & meth_pos >= start(dmr_gene) & meth_pos <= end(dmr_gene)), ]
    betas <- as.data.frame(betas)
  } else stop("The methylation matrix must be a GRanges(-derived) object")
  betas <- colMeans(betas, na.rm = T)
  return(betas)
}