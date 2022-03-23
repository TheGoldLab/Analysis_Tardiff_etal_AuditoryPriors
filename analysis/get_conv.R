get_conv <- function(lmlist) {

  conv_warning <- lapply(lmlist,FUN=function(x) x@optinfo$conv$lme4$messages)
  conv_warning2 <- lapply(lmlist,FUN=function(x) x@optinfo$conv$warnings)
  conv_sing <- sapply(lmlist,FUN=isSingular)
  conv_all <- !sapply(conv_warning,is.null) | !sapply(conv_warning2,is.null) | conv_sing

}