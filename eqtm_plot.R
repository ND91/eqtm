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
#' @import gridExtra
#' @examples 

plot_eqtm <- function(dmrs_se, meth_data, meth_groups, anno_gr, expr_data, expr_groups, united_groups, bm, index, colors = NULL){
  if(is.null(meth_data)) meth_data <- assays(dmrs_se)$meth_summarized
  if(is.null(expr_data)) expr_data <- assays(dmrs_se)$expr
  
  dmr_se <- dmrs_se[which(mcols(dmrs_se)$index == index),]
  dmr_gr <- rowRanges(dmr_se)
  
  topplot <- grid.grabExpr(genome_plot(dmr_gr = dmr_gr, 
                                       bm = bm))
  meth_plot <- dmr_plot(dmr_gr = dmr_gr, 
                                    meth_data = meth_data, 
                                    anno_gr = anno_gr, 
                                    meth_groups = meth_groups, 
                                    flanks = NULL, 
                                    title = NULL)
  
  expr_plot <- NDlib::transcript_strip_plot(id = as.character(dmr_gr$geneid), 
                                                   counts = expr_data, 
                                                   factor_interest = expr_groups, 
                                                   type = "SE", 
                                                   legend = T, 
                                                   title = as.character(dmr_gr$Symbol), 
                                                   y_lab = "Expr") +
    labs(title = "Expression",
         subtitle = paste0(as.character(dmr_gr$Symbol), " (", as.character(dmr_gr$geneid), ")"))
  
  correl_plot <- cor_plot(dmr_se = dmr_se, 
                          united_groups = united_groups)
  
  bottomplots <- ggplotGrob(ggarrange(meth_plot, expr_plot, correl_plot, nrow = 1, ncol = 3, common.legend = T, align = "hv"))
  
  grid.arrange(topplot, bottomplots, nrow = 2, heights = c(0.75, 1.25))
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
  if(is.null(flanks)) flanks <- width(plotrange)/2
  plotrange <- plotrange + flanks
  
  #Preparing the tracks
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = genome_version, chromosome = seqnames(plotrange))
  genetrack <- BiomartGeneRegionTrack(start = start(plotrange), 
                                    end = end(plotrange), 
                                    chromosome = seqnames(plotrange), 
                                    genome = genome_version,
                                    biomart = bm, 
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
             # cex.title = 1.5,
             # cex.axis = 1.5,
             fontcolor = "black",
             background.title = "darkgray") 
}

dmr_plot <- function(dmr_gr, meth_data, anno_gr, meth_groups, color, flanks = NULL, title = NULL){
  if(is.null(flanks)) flanks <- width(dmr_gr)/2
  
  plotrange <- range(dmr_gr + flanks)
  
  cg_ids <- queryHits(findOverlaps(query = anno_gr, subject = dmr_gr))
  cg_ids <- c(min(cg_ids)-1, cg_ids, max(cg_ids)+1)
  anno_sub <- anno_gr[cg_ids,]
  
  meth_df <- data.frame(pos = start(anno_sub), 
                        meth_data[names(anno_sub),])
  
  meth_df_melt <- melt(meth_df, id = "pos")
  meth_df_melt$Group <- rep(meth_groups, each = length(cg_ids))
  
  ggplot(meth_df_melt, aes(x = pos, y = value, col = Group, group = Group)) +
    geom_point() +
    stat_summary(fun.y=mean, geom="smooth") +
    coord_cartesian(xlim = c(start(plotrange), end(plotrange))) +
    ylim(0,1) +
    ylab("Beta") +
    labs(title = "Methylation",
         subtitle = as.character(dmr_gr)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 12), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "none")
}

cor_plot <- function(dmr_se, united_groups){
  
  subtext <- paste0("r = ", 
                    round(rowRanges(dmr_se)$cor_coef, 2), 
                    " [", 
                    round(rowRanges(dmr_se)$CI95_lower, 2), 
                    ", ", 
                    round(rowRanges(dmr_se)$CI95_upper, 2),
                    "]; p-value = ",
                    rowRanges(dmr_se)$pval)
  
  
  plot_df <- data.frame(expr = as.vector(assays(dmr_se)$expr),
                        meth = as.vector(assays(dmr_se)$meth_summarized),
                        Group = united_groups)
  plot_obj <- ggplot(plot_df, aes(x = expr, y = meth)) +
    geom_point(aes(fill = Group), color = "black", shape = 21, size = 5) +
    geom_smooth(method = 'lm', formula = y ~ x, se = F) +
    theme_bw() +
    ylim(0, 1) +
    xlab("Expr") +
    ylab("Meth") +
    labs(title = "Correlation",
         subtitle = subtext) +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12), 
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "none")
}

gg_color_hue <- function(factor_interest) {
  
  unique_factor <- length(unique(factor_interest))
  
  hues <- seq(15, 375, length = unique_factor + 1)
  hcl(h = hues, l = 65, c = 100)[1:unique_factor]
}
