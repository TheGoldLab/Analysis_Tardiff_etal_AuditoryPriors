emm_to_df <- function(emlist,id='model',.confint=T,id2num=F) {

if (is.null(names(emlist))) stop('emm_to_df only works with named lists!')  
    
isList = F

if (class(emlist[[1]]) == 'emmGrid') {
  emlist.sum <- plyr::ldply(emlist,
                            function(x) data.table::as.data.table(x,keep.rownames=F),
                            .progress="text")#,.id = id)
} else if (class(emlist[[1]]) == 'list') {
  emlist.sum <- plyr::ldply(emlist,
                            function(x) dplyr::bind_rows(
                                plyr::llply(x,
                                   function(y) data.table::as.data.table(y,keep.rownames=F))
                                ),
                            .progress="text")
  isList = T
} else {
  error('invalid data type!')
}

#clean up
orig_names <- c("estimate","t.ratio","p.value")
new_names <- c("B","t","p")
names(emlist.sum)[names(emlist.sum) %in% orig_names] <- new_names

names(emlist.sum)[names(emlist.sum)=='.id'] = id

#figure out if we have emmeans or contrasts
if ("contrast" %in% names(emlist.sum)) {
  emtype ="contrast"
} else if ("emmean" %in% names(emlist.sum)){
  emtype = "emmean"
} else {
  warning("Cannot determine emmeans type (emmean or contrast)--assuming is third column name")
  emtype = names(emlist.sum)[3]
}

if (.confint) {
  if (isList) {
    emlist.confint <- plyr::llply(emlist,function(x) plyr::llply(x,confint,.progress="text"))
    
    emlist.confint <- plyr::ldply(emlist.confint,
                                 function(x) dplyr::bind_rows(
                                   plyr::llply(x,
                                       function(y) data.table::as.data.table(y,keep.rownames=F))
                                  ),
                               .progress="text")#,.id=id)
  } else {
    emlist.confint <- plyr::llply(emlist,confint,.progress="text")
    emlist.confint <- plyr::ldply(emlist.confint,data.table::as.data.table,keep.rownames=F,
                                  .progress="text")#,.id = id) #convert to DF
  }
  
  orig_names.confint <- c(".id","lower.CL","upper.CL")
  new_names.confint <- c(id,"CI.lower","CI.upper")
  names(emlist.confint)[names(emlist.confint) %in% orig_names.confint] <- new_names.confint

  #merge in CIs
  emlist.sum <- merge(emlist.sum,emlist.confint[,c(emtype,new_names.confint)],
                       by=c(id,emtype),sort=F)
}

#final cleanup
names(emlist.sum)[names(emlist.sum)==emtype] <- "param"
if (id2num) {
  emlist.sum[id] = as.numeric(emlist.sum[[id]])
}

return(emlist.sum)
}