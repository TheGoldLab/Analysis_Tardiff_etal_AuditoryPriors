sum_to_df <- function(lmlist,id='model',.confint=T,.pcorr=F,id2num=F,type="sum",anova_type="III") {

  if (is.null(names(lmlist))) stop('sum_to_df only works with named lists!')
  
  if (type=="sum") {
    lmlist.sum <- plyr::ldply(lmlist,
          function(x) data.table::as.data.table(summary(x)$coefficients,keep.rownames="param"),
          .progress="text")#,.id = id)

    #clean up
    if (is(lmlist[[1]],"merMod")) {    
      orig_names <- c("Estimate","Std. Error","df","t value","Pr(>|t|)")
      new_names <- c("B","SE","df","t","p")
    } else if (is(lmlist[[1]],c("lm","lmrob"))) {
      orig_names <- c("Estimate","Std. Error","t value","Pr(>|t|)")
      new_names <- c("B","SE","t","p")
    } else {
        warning('Unrecognized lm class. Column names not changed.')
    }
    names(lmlist.sum)[names(lmlist.sum) %in% orig_names] <- new_names
    
  } else if (type=="anova") {
    if (all(sapply(lmlist,function(x) is(x,"merMod")))) {
      lmlist.sum <- plyr::ldply(lmlist,
            function(x) data.table::as.data.table(anova(x,type=anova_type),keep.rownames="param"),
            .progress="text")
      
    } else if (all(sapply(lmlist,function(x) is(x,"lm")))) {
      lmlist.sum <- plyr::ldply(lmlist,
          function(x) data.table::as.data.table(car::Anova(x,type=anova_type),keep.rownames="param"),
          .progress="text")
    }
    else {
        stop("type=anova currently implemnted only for LMs/LMMs!")
    }
    
    #clean up
    orig_names <- c("F value","Pr(>F)")
    new_names <- c("F","p")
    names(lmlist.sum)[names(lmlist.sum) %in% orig_names] <- new_names
    
  } else {
    stop('Incorrect type option (sum,anova).')
  }
  
  #names(lmlist.sum)[names(lmlist.sum)=='.id'] = id

  
  if (.confint & type=="sum") {
    if (all(sapply(lmlist,function(x) is(x,"merMod")))) {
      lmlist.confint <- plyr::llply(lmlist,confint,method="Wald",.progress="text")
      lmlist.confint <- plyr::ldply(lmlist.confint,data.table::as.data.table,keep.rownames="param",
                              .progress="text") #,.id = id) #convert to DF
      #names(lmlist.confint)[names(lmlist.confint)=='.id'] = id
      ###names(lmlist.confint)[3:4] <- c("CI.lower","CI.upper")
      lmlist.confint <- lmlist.confint[!grepl(".sig*",lmlist.confint$param),] #discard ranef intervals
    } else {
      lmlist.confint <- plyr::llply(lmlist,confint,.progress="text")
      lmlist.confint <- plyr::ldply(lmlist.confint,data.table::as.data.table,keep.rownames="param",
                                    .progress="text") #,.id = id) #convert to DF
    }
    names(lmlist.confint)[3:4] <- c("CI.lower","CI.upper")
    #print(head(lmlist.confint))
    #print(head(lmlist.sum))
    #merge in CIs
    lmlist.sum <- merge(lmlist.sum,lmlist.confint,by=c(".id","param"),sort=F)
  }
  
  #clean up
  names(lmlist.sum)[names(lmlist.sum)=='.id'] = id
  if (id2num) {
    lmlist.sum[id] = as.numeric(lmlist.sum[[id]])
  }
  
  if (.pcorr) {
  	lmlist.sum <- plyr::ddply(lmlist.sum,.(param),mutate,
                                p_corr=p.adjust(p,method="holm"))
  	#this ridiculous .data[[id]] is needed for character vector since I didn't tidy this function
  	#https://dplyr.tidyverse.org/articles/programming.html
  	lmlist.sum <- dplyr::arrange(lmlist.sum,.data[[id]],param) 
  }
  
  return(lmlist.sum)

}