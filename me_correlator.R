#' me_correlator 
#'
#' This function was setup to perform an expression quantitative trait methylation analysis for methylation and expression data obtained from the same samples.
#' @param dmr_data A GRanges object whose coordinates correspond to the coordinates of the DMR and whose metadata represents the aggregated methylation value per sample. Can be generated using the dmr_beta_mean.
#' @param expr_data A dataframe containing transcription with a transcript/gene identifier as row names and sample names as column names. Note that the sample names must overlap with the sample names of dmr_data.
#' @param dmr_gene_anno A GRanges object whose coordinates correspond to the coordinates of the DMR and whose metadata contains a column containing the transcript/gene identifier associated to the DMR.
#' @param id_col_name The name of the metadata column of dmr_gene_anno that contains the transcript/gene identifier.
#' @param meth_data A GRanges(-derived) object (such as BSseq) necessary for the permutation analysis 
#' @param iterations An integer representing the number of bootstraps and permutations to perform.
#' @param cor_type A string representing which type of correlation to calculate using base::cor(), choose from "pearson", "spearman", "kendall".
#' @param ncores An integer representing the number of cores to use for parallelization. Defaults to the singlecore implementation. Note that the multicore implementation has the tendency to crash if the datasets to work with are large (i.e. WGBS)
#' @param seed An integer to seed the randomization process. Useful for reproducibility purposes.
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

me_correlator <- function(dmr_data, expr_data, dmr_gene_anno, id_col_name, meth_data, iterations, cor_type = c("pearson", "kendall", "spearman"), ncores = 1, seed = NULL){
  #Argument checking
  if(is.null(dmr_data)) stop("dmr_data cannot be found")
  if(class(dmr_data) != "GRanges") stop("dmr_data must be a GRanges object")
  if(is.null(expr_data)) stop("expr_data cannot be found")
  if(!colnames(mcols(dmr_data)) %in% colnames(expr_data)) stop("There exists no overlap between the samples present in dmr_data and expr_data as evident from the column names")
  if(is.null(id_col_name)) stop("id_col_name cannot be found")
  if(is.null(meth_data)) stop("meth_data cannot be found")
  if(!class(meth_data) %in% c("BSseq", "GRanges")) stop("meth_data must be a GRanges or BSseq object")
  if(is.null(iterations)) stop("iterations cannot be found")
  
  iterations <- round(iterations)
  cor_type <- match.arg(cor_type)
  
  if(!is.numeric(ncores)) stop("ncores must be a numeric")
  if(ncores > 1){
    cat(paste0("Parallelization will be implemented using ", ncores, " cores.\n"))
    #Parallelization under UNIX environments works better if forked (better memory management). However, forking is not an option under windows
    if(Sys.info()['sysname'] != "Windows"){
      cl_type <- "FORK"
    } else{
      cl_type <- "PSOCK"
    }
  } else cat(paste0("Computation will be performed on a single core.\n"))
  
  #Basic QC expression
  expr_data <- expr_data[rowSums(expr_data) != 0, ]
             
  #Find overlapping samples
  dmr_samples <- as.character(colnames(mcols(dmr_data)))
  expr_samples <- as.character(colnames(expr_data))
  overlapping_samples <- intersect(dmr_samples, expr_samples)
             
  #Find DMR-annotations for which aggregated methylation values exist
  overlapping_dmr_gene_anno <- subsetByOverlaps(dmr_gene_anno, dmr_data)
  #Find overlapping transcription features
  overlapping_dmr_gene_anno <- overlapping_dmr_gene_anno[mcols(overlapping_dmr_gene_anno)[, id_col_name] %in% rownames(expr_data), ]
             
  dmr_data_culled <- dmr_data[, overlapping_samples]
  expr_data_culled <- expr_data[, overlapping_samples]
             
  cat("Pulling up Baron Munchausen by his bootstraps\n")
  if(ncores == 1){
    if(!is.null(seed)) set.seed(seed)
    correlations <- sapply(X = overlapping_dmr_gene_anno, simplify = T, USE.NAMES = F, FUN = function(dmr_gene_entry){
    cor_df <- data.frame(meth = as.vector(unlist(mcols(subsetByOverlaps(dmr_data_culled, dmr_gene_entry)))),
                         expr = as.vector(unlist(expr_data_culled[unlist(mcols(dmr_gene_entry[, id_col_name])), ])))
    ci_bootstrapper(cor_df = cor_df, iterations = iterations, cor_type = cor_type)
  })} else{
    cl <- makeCluster(spec = cores, type = cl_type)
    clusterExport(cl, c("dmr_data_culled", "dmr_gene_entry", "id_col_name", "iterations", "cor_type", "ci_bootstrapper"), envir = environment())
    clusterEvalQ(cl, library(boot))
    if(!is.null(seed)) clusterSetRNGStream(cl = cl, iseed = seed)
    
    correlations <- parSapply(cl = cl, X = overlapping_dmr_gene_anno, simplify = T, USE.NAMES = F, FUN = function(dmr_gene_entry){
      cor_df <- data.frame(meth = as.vector(unlist(mcols(subsetByOverlaps(dmr_data_culled, dmr_gene_entry)))),
                           expr = as.vector(unlist(expr_data_culled[unlist(mcols(dmr_gene_entry[, id_col_name])), ])))
      ci_bootstrapper(cor_df = cor_df, iterations = iterations, cor_type = cor_type)
    })
    stopCluster(cl)
    gc()
  }
  
  correlations <- t(correlations)
  correlations_df <- cbind(data.frame(overlapping_dmr_gene_anno), correlations)
  correlations_gr <- makeGRangesFromDataFrame(correlations_df, keep.extra.columns = T)
  correlations_gr$nCpGs <- countOverlaps(correlations_gr, meth_data)
  return(correlations_gr)
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
