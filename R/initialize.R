setMethod("initialize", "mcprofile", function(.Object, CM, estimate, vest, vestunsc, dvest, model, SRDP, fsplines, control, method, ...){
  if (any(class(CM) == "contrMat")) class(CM) <- "matrix" 
  .Object@CM <- CM
  .Object@estimate <- estimate
  .Object@vest <- vest
  .Object@vestunsc <- vestunsc
  .Object@dvest <- dvest
  .Object@model <- model
  .Object@SRDP <- SRDP
  .Object@fsplines <- fsplines
  .Object@control <- control
  .Object@method <- method
  return(.Object)
})

setMethod("initialize", "mcprofileRatio", function(.Object, CMn, CMd, estimate, model, SRDP, fsplines, control, method, ...){
  if (any(class(CMn) == "contrMat")) class(CMn) <- "matrix"
  if (any(class(CMd) == "contrMat")) class(CMd) <- "matrix"
  .Object@CMn <- CMn
  .Object@CMd <- CMd
  .Object@estimate <- estimate
  .Object@model <- model
  .Object@SRDP <- SRDP
  .Object@fsplines <- fsplines
  .Object@control <- control
  .Object@method <- method
  return(.Object)
})


