#' plot_eqtm
#'
#' This function was setup to perform an expression quantitative trait methylation analysis for methylation and gene expression data obtained from the same samples.
#' @param dmrs_se A SummarizedExperiment object generated from the eqtm() function
#' @param sort.by What to sort by? ("pval", "cor_coef", "CI95_diff")
#' @param volcanoplot A boolean whether or not to generate a volcano plot (x = correlation, y = -log10(p-value))
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return 
#' @keywords eqtm, methylation, expression
#' @export
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import ggplot2
#' @import gviz
#' @examples 

plot_eqtm <- function(dmrs_se, meth_data, expr_data, index, factor_interest, colors = NULL){
  if(is.null(meth_data)) meth_data <- assays(eqtms)$meth_summarized
  if(is.null(meth_data)) expr_data <- assays(eqtms)$expr
  
  
  
}

genome_plot <- function(dmr_gr, txdb, flanks = NULL, genome_version = "hg19", title = NULL){
  
  
  
  
  start_plot = start - flanks
  end_plot = end + flanks
  
  #Preparing the tracks
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = genome_version, chromosome = chr)
  genetrack <- GeneRegionTrack(txdb, chromosome = chr, from = start_plot, to = end_plot, geneSymbol = T, name = "Ensemble", fontsize = 14)
  methtrack <- AnnotationTrack(range = dmrs.gr, genome = genome_version, name = "DMR", transcriptAnnotation = "", fontsize = 14)
  
  #Plotting the tracks
  plotTracks(list(gtrack, itrack, genetrack, methtrack), 
             chromosome = chr, 
             from = start_plot, 
             to = end_plot,
             cex.title = 1.5,
             cex.axis = 1.5,
             fontcolor = "black",
             background.title = "darkgray") 
}