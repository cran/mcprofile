setClass("mcprofile", representation(CM = "matrix",
                                     estimate = "numeric",
                                     vest = "numeric",
                                     vestunsc = "numeric",
                                     dvest = "numeric",
                                     model = "list",
                                     SRDP = "list",
                                     fsplines = "list",
                                     bsplines = "list",
                                     control = "list",
                                     method = "character"))


setClass(Class="mcpconfint", representation=representation(confint = "data.frame",
                               quantile = "numeric",
                               estimate = "numeric",
                               alternative="character",
                               adjust="character"))


setClass(Class="mcptest", representation=representation(stat = "numeric",
                            pvalues = "numeric",
                            margin = "numeric",
                            estimate = "numeric",
                            alternative="character",
                            adjust="character"))




setClass("mcprofileRatio", representation(CMn = "matrix",
                                          CMd = "matrix",
                                          estimate = "numeric",
                                          model = "list",
                                          SRDP = "list",
                                          fsplines = "list",
                                          bsplines = "list",
                                          control = "list",
                                          method = "character"))   


