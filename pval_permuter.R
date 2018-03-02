#' pval_permuter 
#'
#' This function was setup to perform an expression quantitative trait methylation analysis for methylation and expression data obtained from the same samples.
#' @param eqtms A GRanges object whose coordinates correspond to the coordinates of the DMR and whose metadata contains a correlation coefficient as calculated using the eqtm function.
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

pval_permuter <- function(eqtm, meth_data, expr_data, iterations = 1000, alternative = c("two.sided", "greater", "less"), id_col_name, cor_type, ncores = 1, seed = NULL){
  if(is.null(eqtm)) stop("eqtm cannot be found")
  if(class(eqtm) != "GRanges") stop("eqtm must be a GRanges object")
  if(is.null(meth_data)) stop("meth_data cannot be found")
  if(!class(meth_data) %in% c("BSseq", "GRanges")) stop("meth_data must be a GRanges or BSseq object")
  if(is.null(colnames(meth_data))) stop("meth_data must have sample names as column names")
  if(is.null(expr_data)) stop("expr_data cannot be found")
  if(is.null(id_col_name) | !id_col_name %in% colnames(mcols(eqtm))) stop("Cannot find the expression ID column in eqtm")
  if(is.null(cor_type)) stop("cor_type cannot be found")
  if(!is.numeric(ncores)) stop("cores must be a numeric")
  if(ncores > 1){
    #Parallelization under UNIX environments works better if forked (better memory management). However, forking is not an option under windows
    if(Sys.info()['sysname'] != "Windows"){
      cl_type <- "FORK"
    } else{
      cl_type <- "PSOCK"
    }
  }
  
  alternative <- match.arg(alternative)
  
  overlapping_samples <- intersect(colnames(meth_data), colnames(expr_data))
  if(length(overlapping_samples) == 0 ) stop("No overlapping samples found for meth_data and expr_data")
  
  ## Generate null distribution
  permutation_df <- unique(data.frame(chr = as.character(seqnames(eqtm)),
                                      CpGs = eqtm$nCpGs))
  permutation_df <- permutation_df[order(permutation_df$chr, permutation_df$CpGs), ]
  chromosome_list <- with(permutation_df, split(permutation_df, chr))
  
  #Split per chromosome for computational purposes
  permuted_cor <- lapply(names(chromosome_list), function(chromosome){
    cat(paste0("Starting on chromosome ", chromosome, ".\n"))
    chr_data <- meth_data[which(seqnames(meth_data) == chromosome), ]
    dmr_length <- chromosome_list[[chromosome]]
    
    #Find a random set of $nCpG indices $iterations times
    if(ncores == 1){
      if(!is.null(seed)) set.seed(seed)
      permuted_cor <- lapply(dmr_length$CpGs, function(CpG){
        eff_nrow <- nrow(chr_data)-as.numeric(CpG)
        indices <- sample(x = 1:eff_nrow, size = iterations, replace = T)
        
        nulldist <- t(sapply(X = indices, simplify = T, FUN = function(index){
          colMeans(getMeth(chr_data[index:(index+CpG-1), overlapping_samples], type = "smooth"), na.rm = T)
        }))
        return(nulldist)
      })
    } else{
      cl <- makeCluster(spec = ncores, type = cl_type)
      if(!is.null(seed)) clusterSetRNGStream(cl = cl, iseed = seed)
      clusterExport(cl, c("chr_data", "iterations", "overlapping_samples"), envir = environment())
      clusterEvalQ(cl, library(bsseq))
      
      permuted_cor <- parLapply(cl = cl, X = dmr_length$CpGs, fun = function(CpG){
        eff_nrow <- nrow(chr_data)-as.numeric(CpG)
        indices <- sample(x = 1:eff_nrow, size = iterations, replace = T)
        
        nulldist <- t(sapply(X = indices, simplify = T, FUN = function(index){
          colMeans(getMeth(chr_data[index:(index+CpG-1), overlapping_samples], type = "smooth"), na.rm = T)
        }))
        return(nulldist)
      })
      stopCluster(cl)
      gc()
    }
    names(permuted_cor) <- paste0(dmr_length$chr, "_", dmr_length$CpGs)
    return(permuted_cor)
  })
  permuted_cor <- unlist(permuted_cor, recursive = F)
  permuted_cor_chr <- gsub("(.+)_.+", "\\1", names(permuted_cor))
  permuted_cor_nCpGs <- gsub(".+_(.+)", "\\1", names(permuted_cor))
  
  ## Calculate probabilities
  cat("Calculating the probabilities under the null (aka \"pvalues\")\n")
  if(ncores == 1){
    pvals <- lapply(eqtm, function(cor_entry){
      perm_dmrs <- unlist(permuted_cor[[which(permuted_cor_chr == seqnames(cor_entry) & permuted_cor_nCpGs == cor_entry$nCpGs)]][, overlapping_samples])
      transcription <- unlist(expr_data[unlist(mcols(cor_entry)[id_col_name]), overlapping_samples])
      null_corcoefdist <- apply(X = perm_dmrs, MARGIN = 1, FUN = function(perm_dmr){
        cor(perm_dmr, transcription, method = cor_type)
      })
      
      if(alternative == "two.tailed"){
        pval <- mean(abs(cor_entry$corcoef) > abs(null_corcoefdist), na.rm = T)
      } else if(alternative == "greater"){
        pval <- mean(cor_entry$corcoef > null_corcoefdist, na.rm = T)
      } else if(alternative == "less"){
        pval <- mean(cor_entry$corcoef < null_corcoefdist, na.rm = T)
      }
      return(pval)	
    })
  } else{
    cl <- makeCluster(spec = ncores, type = cl_type)
    if(!is.null(seed)) clusterSetRNGStream(cl = cl, iseed = seed)
    clusterExport(cl, c("permuted_cor", "permuted_cor_chr", "permuted_cor_nCpGs", "overlapping_samples", "expr_data", "id_col_name"), envir = environment())
    
    pvals <- parLapply(cl = cl, X = eqtm, fun = function(cor_entry){
      
      perm_dmrs <- unlist(permuted_cor[[which(permuted_cor_chr == seqnames(cor_entry) & permuted_cor_nCpGs == cor_entry$nCpGs)]][, overlapping_samples])
      transcription <- unlist(expr_data[unlist(mcols(cor_entry)[id_col_name]), overlapping_samples])
      
      null_corcoefdist <- apply(X = perm_dmrs, MARGIN = 1, FUN = function(perm_dmr){
        cor(perm_dmr, transcription, method = cor_type)
      })
      
      if(alternative == "two.tailed"){
        pval <- mean(abs(cor_entry$corcoef) < abs(null_corcoefdist), na.rm = T)
      } else if(alternative == "greater"){
        pval <- mean(cor_entry$corcoef < null_corcoefdist, na.rm = T)
      } else if(alternative == "less"){
        pval <- mean(cor_entry$corcoef > null_corcoefdist, na.rm = T)
      }
      return(pval)	
    })
  }
  stopCluster(cl)
  gc()
  eqtm$pvals <- unlist(pvals)
  
  #Sorting according to p-value and corcoef
  eqtm <- eqtm[with(eqtm, order(pvals, rev(corcoef), decreasing = F)), ]
  
  return(eqtm)
}