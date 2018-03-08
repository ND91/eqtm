#' topcor
#'
#' This function extracts the 
#' @param dmrs_se A SummarizedExperiment object generated from the eqtm() function
#' @param sort.by What to sort by? ("pval", "cor_coef", "CI95_diff")
#' @param volcanoplot A boolean whether or not to generate a volcano plot (x = correlation, y = -log10(p-value))
#' 
#' @author Andrew Y.F. Li Yim
#' 
#' @return A data.frame containing the eqtms ranked according to the 'sort.by' argument
#' @keywords eqtm, methylation, expression
#' @export
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import ggplot2
#' @examples 

topcor <- function(dmrs_se, sort.by = c("pval", "cor_coef", "CI95_diff"), volcanoplot = T){
  sort.by <- match.arg(sort.by)
  
  eqtms_df <- data.frame(rowRanges(dmrs_se))
  
  #Sort by number of CpGs
  eqtms_df <- eqtms_df[order(eqtms_df$nCpGs, decreasing = T), ]
  switch(sort.by,
         pval = eqtms_df <- eqtms_df[order(eqtms_df$pval),],
         cor = eqtms_df <- eqtms_df[order(abs(eqtms_df$cor_coef)),],
         CI95_diff = eqtms_df <- eqtms_df[order(abs(eqtms_df$CI95_upper - eqtms_df$CI95_lower), decreasing = F),])
  
  if(volcanoplot){
    print(ggplot(eqtms_df, aes(x = cor_coef, y = -log10(pval))) +
            geom_point() +
            theme_bw() +
            xlim(-1,1))
  }
  
  return(eqtms_df)
}