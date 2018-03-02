setClass(Class = "eqtm", 
         contains = "RangedSummarizedExperiment",
         representation(trans = "function", parameters = "list"))

setValidity("eqtm", function(object){
  msg <- validMsg(NULL, .checkAssayNames(object, c("Cov", "M")))
  msg <- validMsg(msg, .checkAssayClasses(object,
                                          c("Cov", "M", "coef", "se.coef")))
  if(class(rowRanges(object)) != "GRanges")
    msg <- validMsg(msg, sprintf("object of class '%s' needs to have a 'GRanges' in slot 'rowRanges'", class(object)))
  ## benchmarking shows that min(assay()) < 0 is faster than any(assay() < 0) if it is false
  if(is.null(colnames(object)))
    msg <- validMsg(msg, "colnames (aka sampleNames) need to be set")
  M <- assay(object, "M", withDimnames = FALSE)
  Cov <- assay(object, "Cov", withDimnames = FALSE)
  if (!.isHDF5ArrayBacked(M) && !.isHDF5ArrayBacked(Cov)) {
    ## TODO: This check is super expensive if M or Cov is a HDF5Matrix, so
    ##       we skip it for the time being
    if(min(assay(object, "M", withDimnames = FALSE)) < 0)
      msg <- validMsg(msg, "the 'M' assay has negative entries")
    if(min(assay(object, "Cov", withDimnames = FALSE)) < 0)
      msg <- validMsg(msg, "the 'Cov' assay has negative entries")
    if(max(assay(object, "M", withDimnames = FALSE) -
           assay(object, "Cov", withDimnames = FALSE)) > 0.5)
      msg <- validMsg(msg, "the 'M' assay has at least one entry bigger than the 'Cov' assay")
  }
  if(!is.null(rownames(M)) || !is.null(rownames(Cov)) ||
     ("coef" %in% assayNames(object) && !is.null(rownames(assay(object, "coef")))) ||
     ("se.coef" %in% assayNames(object) && !is.null(rownames(assay(object, "se.coef")))))
    warning("unnecessary rownames in object")
  if (is.null(msg)) TRUE else msg
})