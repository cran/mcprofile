setMethod("show", "mcprofile", function(object){
  cat("\n   Multiple Contrast Profiles\n\n")
  print(data.frame(Estimate=round(object@estimate,2), Stderr=round(sqrt(object@vest),2)))
  cat("\n")
})

setMethod("show", "mcprofileRatio", function(object){
  cat("\n   Multiple Contrast Profiles\n\n")
  print(data.frame(Estimate=round(object@estimate,2)))
  cat("\n")
})


setMethod("show", "mcpconfint", function(object){
  cat("\n   mcprofile - Confidence Intervals\n\n")
  cat("Adjustment:  ", object@adjust, "\n")
  cat("Alternative: ", object@alternative, "\n")
  cat("Quantile:    ", round(object@quantile, 2), "\n\n")
  dat <- data.frame(round(object@estimate,2), round(object@confint,2))
  if (object@alternative == "two.sided") names(dat) <- c("estimate", "lower", "upper")
  if (object@alternative == "less") names(dat) <- c("estimate", "upper")
  if (object@alternative == "greater") names(dat) <- c("estimate", "lower")
  print(dat)
  cat("\n")
})

setMethod("show", "mcptest", function(object){
  cat("\n   mcprofile - Multiple Testing\n\n")
  cat("Adjustment:  ", object@adjust, "\n")
  cat("Alternative: ", object@alternative, "\n\n")
  dat <- data.frame(round(object@estimate,2), round(object@stat,2), round(object@pvalues,3))
  names(dat) <- c("estimate", "statistic", "adj.pvalue")
  print(dat)
  cat("\n")
})

