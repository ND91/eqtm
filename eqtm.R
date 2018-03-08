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
#' @import SummarizedExperiment
#' @import boot
#' @import doParallel
#' @import foreach
#' @import ggplot2
#' @examples 

eqtm <- function(dmrs_gr, gene_col, meth_data, expr_data, meth_anno_gr, cor_method = c("pearson", "kendall", "spearman"), alternative = c("two.sided", "greater", "less"), N = 1000, iseed = NULL, ncores = 1, verbose = T){
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

  cor_method <- match.arg(cor_method)
  alternative <- match.arg(alternative)
  iseed <- as.numeric(iseed)
  ncores <- as.numeric(ncores)
  
  #Subsetting the datasets according to the overlapping genes
  dmr_genes <- unique(data.frame(chr = as.character(seqnames(dmrs_gr)),
                                 start = start(dmrs_gr),
                                 end = end(dmrs_gr),
                                 geneid = as.character(mcols(dmrs_gr)[, gene_col])))
  ol_genes <- intersect(dmr_genes$geneid, rownames(expr_data))
  
  dmr_genes_ol <- makeGRangesFromDataFrame(dmr_genes[dmr_genes$geneid %in% ol_genes,], keep.extra.columns = T)
  dmr_genes_ol$index <- names(dmr_genes_ol)
  expr_ol <- expr_data[rownames(expr_data) %in% ol_genes, ol_samples]
  meth_anno_ol <- meth_anno_gr[unique(subjectHits(findOverlaps(dmr_genes_ol, meth_anno_gr))),]
  meth_bg <- meth_data[, ol_samples]
  
  #Methylation aggregation
  if(verbose) cat(paste0(Sys.time(), " Aggregating the DMRs.\n"))
  dmrs_meth <- dmr_aggregator(dmrs_gr = dmr_genes_ol, meth_data = meth_bg[names(meth_anno_ol),], meth_anno = meth_anno_ol, expr_data = expr_ol)
  
  #Correlation and bootstrapping
  if(verbose) cat(paste0(Sys.time(), " Correlating the DMRs with gene expression.\n"))
  dmrs_cor <- dmr_correlator(dmrs_se = dmrs_meth, cor_method = cor_method, N = N, ncores = ncores, iseed = iseed, verbose = T)
  
  #Permutation
  if(verbose) cat(paste0(Sys.time(), " Calculating the p-values for the correlation coefficient.\n"))
  dmrs_data <- dmr_permutations(dmrs_se = dmrs_cor, meth_bg = meth_bg, meth_anno_gr = meth_anno_ol, cor_method = cor_method, N = N, ncores = ncores, alternative = alternative, iseed = iseed, verbose = T)
  
  return(dmrs_data)
}

dmr_aggregator <- function(dmrs_gr, meth_data, meth_anno, expr_data){
  #TODO: Add a way of calculating the occupancy of the methylation rather than the mean.
  meth_fo <- data.frame(findOverlaps(dmrs_gr, meth_anno))
  meth_fo_split <- split(meth_fo, meth_fo$queryHits)
  dmr_mean <- lapply(meth_fo_split, FUN = function(dmr){
    colMeans(meth_data[dmr$subjectHits,], na.rm = T)
  })
  dmr_mean <- do.call(rbind, dmr_mean)
  names(dmrs_gr) <- rownames(dmr_mean)
  expr_rn <- expr_data[as.character(dmrs_gr$geneid), ]
  rownames(expr_rn) <- rownames(dmr_mean)
  
  agg_dmr_se <- SummarizedExperiment(assays = list(meth_summarized = dmr_mean,
                                                   expr = expr_rn),
                                     rowRanges = dmrs_gr)
  rownames(agg_dmr_se) <- rowData(agg_dmr_se)$geneid
  
  return(agg_dmr_se)
}

dmr_correlator <- function(dmrs_se, cor_method, N, ncores, iseed, verbose){
  expr_entry <- assays(dmrs_se)$expr
  dmr_entry <- assays(dmrs_se)$meth_summarized
  
  cor_boot <- function(data, k, cor_method){
    cor(data[k,1], data[k,2], method = cor_method)
  } 
  
  bootstrapper <- function(i){
    boot_obj <- boot::boot(data = cbind(expr_entry[i,], dmr_entry[i,]), 
                           statistic = cor_boot, 
                           cor_method = cor_method,
                           R = N)
    CI95 <- boot::boot.ci(boot.out = boot_obj, 
                          type = "bca")
    return(c(cor_coef = boot_obj$t0, CI95_lower = CI95$bca[4], CI95_upper = CI95$bca[5]))
  }
  
  if(verbose) cat(paste0(Sys.time(), "\tPerforming ", N, " bootstraps for the ", cor_method, " correlation coefficient on ", ncores, " core(s).\n"))
  if(ncores == 1){
    if(!is.null(iseed)) set.seed(iseed)
    
    cor_vals <- sapply(X = seq.int(nrow(expr_entry)), FUN = bootstrapper)
  } else{
    cl <- makeCluster(ncores)
    if(!is.null(iseed)) clusterSetRNGStream(cl = cl, iseed = iseed)
    
    #BUG: Sometimes clusterExport cannot find "expr_entry" in the environment
    clusterExport(cl = cl, varlist = list("expr_entry", "dmr_entry", "cor_boot", "cor_method", "N"), envir = environment())
    
    cor_vals <- parSapply(cl = cl, X = seq.int(nrow(expr_entry)), FUN = bootstrapper)
    stopCluster(cl)
  }
  if(verbose) cat(paste0(Sys.time(), "\tFinished bootstrapping.\n"))
  
  #Format the results
  mcols(dmrs_se) <- data.frame(mcols(dmrs_se), t(cor_vals))
  
  return(dmrs_se)
}

dmr_permutations <- function(dmrs_se, meth_bg, meth_anno_gr, cor_method, N, ncores, iseed, verbose, alternative){
  #Find the number of CpGs
  ncpgs <- as.numeric(table(queryHits(findOverlaps(dmrs_se, meth_anno_gr))))
  if(length(ncpgs) != nrow(dmrs_se)) stop("Some DMRs cannot be found in the provided annotation file.")
  mcols(dmrs_se)$nCpGs <- ncpgs
  
  #Permutation analyses
  dmrs_split <- split(rowRanges(dmrs_se), seqnames(dmrs_se))
  meth_bg_split <- split(meth_bg, seqnames(meth_anno_gr))
  
  ncpgs_chrom <- lapply(dmrs_split, function(chrom){
    sort(as.vector(unique(chrom$nCpGs)))
  })
  
  if(verbose) cat(paste0(Sys.time(), "\tGenerating ", length(unlist(ncpgs_chrom)), " sets of ", N, " background DMRs.\n"))
  ran_dist <- lapply(names(ncpgs_chrom), function(chrom, ncpgs_chrom, meth_bg_split, N){
    bg_cpgs_chrom <- meth_bg_split[[chrom]]
    ncpgs <- ncpgs_chrom[[chrom]]
    
    null_chrom <- lapply(ncpgs, function(ncpg, bg_cpgs_chrom, N){
      eff_bg <- nrow(bg_cpgs_chrom)-ncpg
      
      #The -1 is to offset the total such that the actual correlation can still be included. A p-value cannot be 0.
      start_nulls <- sample(x = eff_bg, size = N-1)
      
      null_agg <- t(sapply(start_nulls, function(start_null, ncpg, bg_cpgs_chrom){
        end_null <- start_null+ncpg-1
        null_vals <- bg_cpgs_chrom[start_null:end_null,]
        
        #TODO: add other aggregating kernels here
        null_agg <- colMeans(null_vals)
        return(null_agg)
      }, ncpg = ncpg, bg_cpgs_chrom = bg_cpgs_chrom))
      return(null_agg)
    }, bg_cpgs_chrom = bg_cpgs_chrom, N = N)
    names(null_chrom) <- ncpgs
    return(null_chrom)
  }, ncpgs_chrom = ncpgs_chrom, meth_bg_split = meth_bg_split, N = N)
  names(ran_dist) <- names(ncpgs_chrom)
  
  if(verbose) cat(paste0(Sys.time(), "\tGenerating the null distribution and calculating the p-values.\n"))
  
  nd_generator <- function(iterator, dmrs_chrom_df, randmrs_chrom_df, expr_data, cor_method, alternative){
    #Helper function to correlate the gene expression with the randomly generated DMRs methylation and calculate the p-value
    ncpgs <- dmrs_chrom_df$nCpGs[iterator]
    geneid <- dmrs_chrom_df$geneid[iterator]
    permuted <- randmrs_chrom_df[[as.character(ncpgs)]]
    
    null_dist <- apply(X = permuted, MARGIN = 1, FUN = function(meth_data, expr_data, cor_method){
      cor(x = expr_data, y = meth_data, method = cor_method)
    }, expr_data = expr_data[geneid,], cor_method = cor_method)
    
    #P-values cannot be 0
    null_dist <- c(null_dist, dmrs_chrom_df$cor_coef[iterator])
    
    pval <- switch(alternative, 
                   two.sided = mean(abs(null_dist) >= abs(dmrs_chrom_df$cor_coef[iterator]), na.rm = T), 
                   greater = mean(null_dist >= dmrs_chrom_df$cor_coef[iterator], na.rm = T),
                   less = mean(null_dist <= dmrs_chrom_df$cor_coef[iterator], na.rm = T))
    
    return(pval)
  }
  
  dmrs_chrom <- foreach(i = 1:length(dmrs_split), .combine = "rbind") %do% {
    chrom <- names(dmrs_split)[i]
    
    dmrs_chrom_df <- data.frame(dmrs_split[[chrom]], stringsAsFactors = F)
    randmrs_chrom_df <- ran_dist[[as.character(chrom)]]
    expr_data <- assays(dmrs_se)$expr
    
    if(ncores == 1){
      if(!is.null(iseed)) set.seed(iseed)
      pvals_chrom <- foreach(j = 1:nrow(dmrs_chrom_df), .combine = "c") %do% {
        nd_generator(iterator = j, dmrs_chrom_df = dmrs_chrom_df, expr_data = expr_data, randmrs_chrom_df = randmrs_chrom_df, cor_method = cor_method)
      }
    } else{
      cl <- makeCluster(ncores)
      if(!is.null(iseed)) clusterSetRNGStream(cl = cl, iseed = iseed)
      clusterExport(cl = cl, varlist = list("nd_generator", "dmrs_chrom_df", "expr_data", "randmrs_chrom_df", "cor_method", "alternative"), envir = environment())
      registerDoParallel(cl)
      
      pvals_chrom <- foreach(j = 1:nrow(dmrs_chrom_df), .combine = "c") %dopar% {
        nd_generator(iterator = j, dmrs_chrom_df = dmrs_chrom_df, expr_data = expr_data, randmrs_chrom_df = randmrs_chrom_df, cor_method = cor_method, alternative = alternative)
      }
      stopCluster(cl)
    }
    dmrs_chrom_df$pvals <- pvals_chrom
    return(dmrs_chrom_df)
  }
  if(verbose) cat(paste0(Sys.time(), "\tFinished calculating the p-values.\n"))
  
  mcols(dmrs_se)$pval <- dmrs_chrom[match(mcols(dmrs_se)$index, dmrs_chrom$index),]$pval
  
  return(dmrs_se)
}

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