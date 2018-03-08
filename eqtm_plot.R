plot_eqtm <- function(dmrs_se, meth_data, expr_data, index, factor_interest, colors = NULL){
  #
  if(is.null(meth_data)) meth_data <- assays(eqtms)$meth_summarized
  if(is.null(meth_data)) expr_data <- assays(eqtms)$expr
  
}