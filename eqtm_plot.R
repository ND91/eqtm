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
#' @import Gviz
#' @import biomaRt
#' @import NDlib
#' @import cowplot
#' @examples 

plot_eqtm <- function(dmrs_se, meth_data, meth_groups, anno_gr, expr_data, expr_groups, united_groups, bm, index, colors = NULL){
  if(is.null(meth_data)) meth_data <- assays(eqtms)$meth_summarized
  if(is.null(expr_data)) expr_data <- assays(eqtms)$expr
  
  dmr_se <- dmrs_se[which(mcols(dmrs_se)$index == index),]
  dmr_gr <- rowRanges(dmr_se)
  
  genome_plot(dmr_gr = dmr_gr, bm = bm)
  topplot <- recordPlot()
  dmr_plot(dmr_gr = dmr_gr, anno_gr = anno_gr, meth_groups = meth_groups, color = gg_color_hue(meth_groups), flanks = NULL, genome_version = "hg19", title = NULL)
  bottomleft <- recordPlot()
  bottommid <- NDlib::transcript_strip_plot(id = dmr_gr$geneid, counts = expr_data, factor_interest = expr_groups, type = "SE", legend = F)
  bottomright <- cor_plot(dmr_se = dmr_se, united_groups = united_groups)
  
  plot_grid(topplot,
            bottomleft,
            bottommid,
            bottomright,
            labels = "AUTO", 
            hjust = 0, 
            vjust = 1,
            scale = c(1., 1., 0.9, 0.9)
            )
}

genome_plot <- function(dmr_gr, bm, flanks = NULL, genome_version = "hg19", title = NULL, ...){
  if(!is.null(bm) & class(bm) != "Mart") stop("'biomart' object is not of type 'Mart'")
  
  ensg_gr <- makeGRangesFromDataFrame(df = getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"), 
                                                 filters = "ensembl_gene_id", 
                                                 values = as.character(dmr_gr$geneid), 
                                                 mart = bm), 
                                      keep.extra.columns = T, 
                                      seqnames.field = "chromosome_name", 
                                      start.field = "start_position", 
                                      end.field = "end_position")
  
  seqlevelsStyle(ensg_gr) <- "UCSC"
  
  plotrange <- range(c(range(dmr_gr), range(ensg_gr)))
  if(is.null(flanks)) flanks <- width(plotrange)
  plotrange <- plotrange + flanks
  
  #Preparing the tracks
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = genome_version, chromosome = seqnames(plotrange))
  genetrack <- BiomartGeneRegionTrack(start = start(plotrange), 
                                    end = end(plotrange), 
                                    chromosome = seqnames(plotrange), 
                                    genome = genome_version,
                                    biomart = biomart, 
                                    name = "ENS",
                                    collapseTranscripts = "meta",
                                    transcriptAnnotation = "symbol", 
                                    #...,
                                    fontsize = 14)
  methtrack <- AnnotationTrack(range = dmr_gr, genome = genome_version, name = "DMR", transcriptAnnotation = "", fontsize = 14, col = "blue")
  tlist <- list(genetrack, methtrack)
  
  hltrack <- HighlightTrack(trackList = tlist, range = ensg_gr + width(plotrange)/100, col = "red")
  
  #Plotting the tracks
  plotTracks(list(gtrack, itrack, hltrack), 
             chromosome = as.character(seqnames(plotrange)), 
             from = start(plotrange), 
             to = end(plotrange),
             cex.title = 1.5,
             cex.axis = 1.5,
             fontcolor = "black",
             background.title = "darkgray") 
}

dmr_plot <- function(dmr_gr, anno_gr, meth_groups, color, flanks = NULL, genome_version = "hg19", title = NULL){
  if(is.null(flanks)) flanks <- width(dmr_gr)
  
  plotrange <- range(dmr_gr + flanks)
  
  cg_ids <- subsetByOverlaps(x = anno_gr, ranges = dmr_gr)
  meth_gr <- GRanges(seqnames = seqnames(cg_ids), ranges = ranges(cg_ids)) 
  mcols(meth_gr) <- meth_data[names(cg_ids),]
  
  #Preparing the tracks
  methtrack <- DataTrack(range = meth_gr, 
                         groups = meth_groups, 
                         genome = genome_version, 
                         name = "% Methylation", 
                         type = c("a", "p"),
                         ylim = c(0,1),
                         lty = sort(as.numeric(as.factor(levels(meth_groups)))),
                         #col = color,
                         cexp = 0.7,
                         fontsize = 10,
                         fontsize = 14)
  
  #Plotting the tracks
  plotTracks(methtrack, 
             chromosome = as.character(seqnames(plotrange)), 
             from = start(plotrange), 
             to = end(plotrange),
             cex.title = 1.5,
             cex.axis = 1.5,
             fontcolor = "black",
             background.title = "darkgray") 
  
  
}

cor_plot <- function(dmr_se, united_groups){
  plot_df <- data.frame(expr = as.vector(assays(dmr_se)$expr),
                        meth = as.vector(assays(dmr_se)$meth_summarized),
                        Groups = united_groups)
  plot_obj <- ggplot(plot_df, aes(x = expr, y = meth)) +
    geom_point(aes(fill = Groups), color = "black", shape = 21, size = 5) +
    geom_smooth(method = 'lm', formula = y ~ x, se = F) +
    theme_bw() +
    ylim(0, 1) +
    xlab("Log expression") +
    ylab("% methylation") +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12), 
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "bottom")
}

gg_color_hue <- function(factor_interest) {
  
  unique_factor <- length(unique(factor_interest))
  
  hues <- seq(15, 375, length = unique_factor + 1)
  hcl(h = hues, l = 65, c = 100)[1:unique_factor]
}
