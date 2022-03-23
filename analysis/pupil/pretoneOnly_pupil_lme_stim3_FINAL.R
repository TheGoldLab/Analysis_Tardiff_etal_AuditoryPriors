#clear memory
rm(list=ls())


### this runs regressions for pretoneOnly stim-aligned pupil data


## LOADING data/libraries ##

#load libraries
library(lme4)
library(lmerTest)
library(car)
library(plyr)
library(dplyr)
library(ggplot2)
#library(afex)
library(emmeans)
emm_options(lmerTest.limit = Inf, lmer.df = "satterthwaite")


switch(Sys.info()[['sysname']],
       Windows = setwd(file.path(
         Sys.getenv('USERPROFILE'),'Dropbox/Goldlab/AuditoryPriors/data processing/pupil')),
       Darwin = setwd('~/Dropbox/Goldlab/AuditoryPriors/data processing/pupil')
)

#path to data files
switch(Sys.info()[['sysname']],
       Windows = DATA_OUT_PATH <- (paste0(
         Sys.getenv('USERPROFILE'),'/OneDrive/Goldlab/AuditoryPriors/cached data/')),
       Darwin = DATA_OUT_PATH <- '~/OneDrive/Goldlab/AuditoryPriors/cached data/'
)
pd_file= paste0(DATA_OUT_PATH,'pupil_data_pretoneOnly_ds50_forR_05-May-2021.csv')
baseline_file = paste0(DATA_OUT_PATH,'pupil_data_pretoneOnly_bl_forR_05-May-2021.csv')

#source my functions if any
source('../sum_to_df.R')
source('../emm_to_df.R')
source('../get_conv.R')


#set temporal parameters [this is easier than pupilnet cause I kept it in actual ms]
stim_period_ms.c <- c(60, 1360)
stim_period_ms.i <- c(60,860)

#this doesn't matter if only running on correct trials
conType = "congruentS" #"congruent"

#load data
pd_df<-read.table(pd_file,sep=',', header=TRUE, stringsAsFactors=FALSE,na.strings = c('NaN'))
bl_df<-read.table(baseline_file,sep=',', header=TRUE, stringsAsFactors=FALSE,na.strings = c('NaN'))

head(pd_df)
head(bl_df)

#cleanup
bl_df$GroupCount <- NULL

#set up variables/factors
my_simple2<-contr.treatment(2,base=2) - matrix(rep(1/2,2))

pd_df$congruent.f <- factor(pd_df[,conType],levels=c(3,2,1,0),
                            labels=c("con","con-incon","incon-con","incon"))
contrasts(pd_df$congruent.f) <- contr.sum(4)
print(conType)
contrasts(pd_df$congruent.f)
unique(pd_df[,c(conType,'congruent.f')])

pd_df$isH.fs <- factor(pd_df$isH,levels=c(1,0),
                      labels=c("high","low"))
contrasts(pd_df$isH.fs) <- my_simple2
contrasts(pd_df$isH.fs)
pd_df$isH.f <- pd_df$isH.fs
contrasts(pd_df$isH.f) <- contr.sum(2)
contrasts(pd_df$isH.f)

#merge in baseline data
pd_df <- left_join(pd_df,bl_df[,c('dataID','trialN','pupilBL2')],by=c('dataID','trialN'))

#get rid of missing before scaling vars (probably already done w/ new preproc)
pd_df <- pd_df[complete.cases(pd_df),]

#get a handle on descriptives to inform scaling
scalevars_stats <- summarise(pd_df,
                             m_bl=mean(pupilBL2,na.rm=T),
                             sd_bl=sd(pupilBL2,na.rm=T),
                             m_x=mean(posXCbl,na.rm=T),
                             sd_x=sd(posXCbl,na.rm=T),
                             m_y=mean(posYCbl,na.rm=T),
                             sd_y=sd(posYCbl,na.rm=T),
                             m_asnr=mean(aSNR,na.rm=T),
                             sd_asnr=sd(aSNR,na.rm=T),
)

#set up variables: centering/scaling
pd_df$zaSNR <- scale(pd_df$aSNR)
pd_df$blz <- scale(pd_df$pupilBL2,center=T,scale=F) #not scaling since already scaled
pd_df$posX <- scale(pd_df$posXCbl)
pd_df$posY <- scale(pd_df$posYCbl)
pd_df$zptlen <- scale(pd_df$pretoneLength)

#set up stimOnT variable (onset of the test tone)
pd_df$trial_time_stimOnT <- pd_df$trial_time_stimOff + 300


### CORRECT TRIALS ####

pd_stimdf.c <- subset(pd_df,success==1 & 
                           trial_time_stimOnT >= stim_period_ms.c[1] & 
                           trial_time_stimOnT <= stim_period_ms.c[2])

if (conType=="congruent") {
  pd_stimdf.cs <- subset(pd_stimdf.c,congruent %in% c(0,3) & aSNR %in% c(.05,.5))
} else {
  pd_stimdf.cs <- subset(pd_stimdf.c,congruentS %in% c(0,3) & aSNR %in% c(.05,.5))
}
unique(pd_stimdf.cs$aSNR)
unique(pd_stimdf.cs[,conType])

pd_stimdf.cs$congruent.f <- factor(pd_stimdf.cs[,conType],levels=c(3,0),
                            labels=c("con","incon"))
contrasts(pd_stimdf.cs$congruent.f) <- contr.sum(2)
contrasts(pd_stimdf.cs$congruent.f)
unique(pd_stimdf.cs[,c(conType,'congruent.f')])

pd_stimdf.cs$aSNR.f <- factor(pd_stimdf.cs$aSNR,levels=c(.5,.05),
                                 labels=c("high SNR","low SNR"))
contrasts(pd_stimdf.cs$aSNR.f) <- contr.sum(2)
contrasts(pd_stimdf.cs$aSNR.f)
unique(pd_stimdf.cs[,c('aSNR','aSNR.f')])


stim_times.c <- sort(unique(pd_stimdf.c$trial_time_stimOnT))

### REGRESSIONS ###

## FIRST: fit no-corr model, reduced dataset to see what we've got

pd_stimdf.cs[,c('congruent.f1','aSNR.f1','congruent_aSNR.f1')] <- 
  model.matrix(~1+pd_stimdf.cs$congruent.f*pd_stimdf.cs$aSNR.f,pd_stimdf.cs)[,2:4]
pd_stimdf.cs$isH.f1 <- 
  model.matrix(~1+pd_stimdf.cs$isH.f,pd_stimdf.cs)[,2]
unique(pd_stimdf.cs[,c("congruent.f","aSNR.f","congruent.f1","aSNR.f1","congruent_aSNR.f1")])
unique(pd_stimdf.cs[,c("isH.f","isH.f1")])

##ZC
if (F) {
  stim.lm.zc.cs <- llply(stim_times.c,
                     function(x) {
                       print(x)
                       this_data <- subset(pd_stimdf.cs,trial_time_stimOnT==x)
                       this_stim_lm <- lmer(pupilCblz2~
                                        congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                          isH.f + blz + posX + posY + zptlen + 
                                        (1 + congruent.f1 + aSNR.f1 + congruent_aSNR.f1 + 
                                           isH.f1 + blz + posX + posY + zptlen||dataID),
                                      data=this_data,
                                      control=lmerControl(optimizer="bobyqa",
                                                          optCtrl=list(maxfun=2e5)))
                     },.progress="text")
  names(stim.lm.zc.cs) <- stim_times.c
  save(stim.lm.zc.cs,file=paste0(DATA_OUT_PATH,'ptO_stim_lm_zc_cs_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_stim_lm_zc_cs_2021-06-09.rda'))
}

conv_all.zc.cs <- get_conv(stim.lm.zc.cs)

#do PCA of rfx and extract loadings/variance explained for each time point
#could just extract variance components since zero correlation but easier to use existing code
rePCA_all <- lapply(stim.lm.zc.cs,rePCA)
pca_rot_bad <- lapply(rePCA_all[conv_all.zc.cs],FUN=function(x) x$dataID$rotation)
#pca_rot_bad <- array(unlist(pca_rot_bad),dim=c(8,8,length(pca_rot_bad)))
pca_sd_bad <- sapply(rePCA_all[conv_all.zc.cs],FUN=function(x) x$dataID$sdev)
pca_sd_bad_all <- rowMeans(pca_sd_bad)


#which rfx components most often load on the least-explanatory PCs? 
#(this is a little janky since averaging)
pca_rot_bad_maxs <- sapply(pca_rot_bad,function(x) apply(x,2,function(y) which.max(abs(y))))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

pca_rot_bad_modes <- apply(pca_rot_bad_maxs,1,Mode)

rfxnames <- names(ranef(stim.lm.zc.cs[[1]])$dataID)
worst_rfx <- rfxnames[pca_rot_bad_modes[pca_sd_bad_all < .1]]


#based on above, first try eliminating aSNR and refit
if (F) {
  stim.lm.zc.cs.2 <- llply(stim_times.c,
                    function(x) {
                      if (conv_all.zc.cs[as.character(x)]) {
                        print(x)
                        this_data <- subset(pd_stimdf.cs,trial_time_stimOnT==x)
                        this_stim_lm <- lmer(pupilCblz2~
                                                 congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                 isH.f + blz + posX + posY + zptlen + 
                                                 (1 + congruent.f1 + congruent_aSNR.f1 + 
                                                    isH.f1 + blz + posX + posY + zptlen||dataID),
                                               data=this_data,
                                               control=lmerControl(optimizer="bobyqa",
                                                                   optCtrl=list(maxfun=2e5)))
                      }
                    },.progress="text")
  names(stim.lm.zc.cs.2) <- stim_times.c
  stim.lm.zc.cs.2 <- stim.lm.zc.cs.2[conv_all.zc.cs]
  save(stim.lm.zc.cs.2,file=paste0(DATA_OUT_PATH,'ptO_stim_lm_zc_cs2_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_stim_lm_zc_cs2_2021-06-09.rda'))
}

conv_all.zc.cs.2 <- get_conv(stim.lm.zc.cs.2)

#based on spot check, remove isH
if (F) {
  stim.lm.zc.cs.3 <- llply(as.integer(names(conv_all.zc.cs.2)),
                              function(x) {
                                if (conv_all.zc.cs.2[as.character(x)]) {
                                  print(x)
                                  this_data <- subset(pd_stimdf.cs,trial_time_stimOnT==x)
                                  this_stim_lm <- lmer(pupilCblz2~
                                                     congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                     isH.f + blz + posX + posY + zptlen + 
                                                     (1 + congruent.f1 + congruent_aSNR.f1 + 
                                                        blz + posX + posY + zptlen||dataID),
                                                     data=this_data,
                                                     control=lmerControl(optimizer="bobyqa",
                                                                         optCtrl=list(maxfun=2e5)))
                                }
                              },.progress="text")
  names(stim.lm.zc.cs.3) <- names(conv_all.zc.cs.2)
  stim.lm.zc.cs.3 <- stim.lm.zc.cs.3[conv_all.zc.cs.2]
  save(stim.lm.zc.cs.3,file=paste0(DATA_OUT_PATH,'ptO_stim_lm_zc_cs3_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_stim_lm_zc_cs3_2021-06-09.rda'))
}

conv_all.zc.cs.3 <- get_conv(stim.lm.zc.cs.3)

#now remove congruent.f (based on spot check)
if (F) {
  stim.lm.zc.cs.4 <- llply(as.integer(names(conv_all.zc.cs.3)),
                              function(x) {
                                if (conv_all.zc.cs.3[as.character(x)]) {
                                  print(x)
                                  this_data <- subset(pd_stimdf.cs,trial_time_stimOnT==x)
                                  this_stim_lm <- lmer(pupilCblz2~
                                                         congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                         isH.f + blz + posX + posY + zptlen + 
                                                         (1 + congruent_aSNR.f1 + 
                                                            blz + posX + posY + zptlen||dataID),
                                                       data=this_data,
                                                       control=lmerControl(optimizer="bobyqa",
                                                                           optCtrl=list(maxfun=2e5)))
                                }
                              },.progress="text")
  names(stim.lm.zc.cs.4) <- names(conv_all.zc.cs.3)
  stim.lm.zc.cs.4 <- stim.lm.zc.cs.4[conv_all.zc.cs.3]
  save(stim.lm.zc.cs.4,file=paste0(DATA_OUT_PATH,'ptO_stim_lm_zc_cs4_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_stim_lm_zc_cs4_2021-06-09.rda'))
}

conv_all.zc.cs.4 <- get_conv(stim.lm.zc.cs.4)


#spot check suggests need to drop interaction for remaining issues
if (F) {
  stim.lm.zc.cs.5 <- llply(as.integer(names(conv_all.zc.cs.4)),
                              function(x) {
                                if (conv_all.zc.cs.4[as.character(x)]) {
                                  print(x)
                                  this_data <- subset(pd_stimdf.cs,trial_time_stimOnT==x)
                                  this_stim_lm <- lmer(pupilCblz2~
                                                         congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                         isH.f + blz + posX + posY + zptlen + 
                                                         (1 + 
                                                            blz + posX + posY + zptlen||dataID),
                                                       data=this_data,
                                                       control=lmerControl(optimizer="bobyqa",
                                                                           optCtrl=list(maxfun=2e5)))
                                }
                              },.progress="text")
  names(stim.lm.zc.cs.5) <- names(conv_all.zc.cs.4)
  stim.lm.zc.cs.5 <- stim.lm.zc.cs.5[conv_all.zc.cs.4]
  save(stim.lm.zc.cs.5,file=paste0(DATA_OUT_PATH,'ptO_stim_lm_zc_cs5_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_stim_lm_zc_cs5_2021-06-09.rda'))
}

conv_all.zc.cs.5 <- get_conv(stim.lm.zc.cs.5)

any(conv_all.zc.cs.5)


## REINTRO CORR IN MODELS THAT FIT
#No converge
# if (T) {
#   stim.lm.cs <- llply(stim_times.c,
#                       function(x) {
#                         if (!conv_all.zc.cs[as.character(x)]) {
#                           print(x)
#                           this_data <- subset(pd_stimdf.cs,trial_time_stimOnT==x)
#                           this_stim_lm <- lmer(pupilCblz2~
#                                                  congruent.f + aSNR.f + congruent.f:aSNR.f +
#                                                  isH.f + blz + posX + posY + zptlen +
#                                                  (1 + congruent.f1 + aSNR.f1 + congruent_aSNR.f1 +
#                                                     isH.f1 + blz + posX + posY + zptlen|dataID),
#                                                data=this_data,
#                                                control=lmerControl(optimizer="bobyqa",
#                                                                    optCtrl=list(maxfun=2e5)))
#                         }
#                       },.progress="text")
#   names(stim.lm.cs) <- stim_times.c
#   stim.lm.cs <- stim.lm.cs[!conv_all.zc.cs]
#   save(stim.lm.cs,file=paste0(DATA_OUT_PATH,'ptO_stim_lm_cs_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_stim_lm_cs_2021-06-03.rda'))
# }
# 
# conv_all.cs <- get_conv(stim.lm.cs)

#noconv
# if (T) {
#   stim.lm.cs.2 <- llply(as.integer(names(conv_all.zc.cs.2)),
#                         function(x) {
#                           if (!conv_all.zc.cs.2[as.character(x)]) {
#                             print(x)
#                             this_data <- subset(pd_stimdf.cs,trial_time_stimOnT==x)
#                             this_stim_lm <- lmer(pupilCblz2~
#                                                    congruent.f + aSNR.f + congruent.f:aSNR.f +
#                                                    isH.f + blz + posX + posY + zptlen +
#                                                    (1 + congruent.f1 + congruent_aSNR.f1 +
#                                                       isH.f1 + blz + posX + posY + zptlen|dataID),
#                                                  data=this_data,
#                                                  control=lmerControl(optimizer="bobyqa",
#                                                                      optCtrl=list(maxfun=2e5)))
#                           }
#                         },.progress="text")
#   names(stim.lm.cs.2) <- names(conv_all.zc.cs.2)
#   stim.lm.cs.2 <- stim.lm.cs.2[!conv_all.zc.cs.2]
#   save(stim.lm.cs.2,file=paste0(DATA_OUT_PATH,'ptO_stim_lm_cs2_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_stim_lm_cs2_2021-06-03.rda'))
# }
# 
# conv_all.cs.2 <- get_conv(stim.lm.cs.2)

#noconv
# if (T) {
#   stim.lm.cs.3 <- llply(as.integer(names(conv_all.zc.cs.3)),
#                         function(x) {
#                           if (!conv_all.zc.cs.3[as.character(x)]) {
#                             print(x)
#                             this_data <- subset(pd_stimdf.cs,trial_time_stimOnT==x)
#                             this_stim_lm <- lmer(pupilCblz2~
#                                                    congruent.f + aSNR.f + congruent.f:aSNR.f +
#                                                    isH.f + blz + posX + posY + zptlen +
#                                                    (1 + congruent.f1 + congruent_aSNR.f1 +
#                                                       blz + posX + posY + zptlen|dataID),
#                                                  data=this_data,
#                                                  control=lmerControl(optimizer="bobyqa",
#                                                                      optCtrl=list(maxfun=2e5)))
#                           }
#                         },.progress="text")
#   names(stim.lm.cs.3) <- names(conv_all.zc.cs.3)
#   stim.lm.cs.3 <- stim.lm.cs.3[!conv_all.zc.cs.3]
#   save(stim.lm.cs.3,file=paste0(DATA_OUT_PATH,'ptO_stim_lm_cs3_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_stim_lm_cs3_2021-06-03.rda'))
# }
# 
# conv_all.cs.3 <- get_conv(stim.lm.cs.3)

#noconv
# if (T) {
#   stim.lm.cs.4 <- llply(as.integer(names(conv_all.zc.cs.4)),
#                         function(x) {
#                           if (!conv_all.zc.cs.4[as.character(x)]) {
#                             print(x)
#                             this_data <- subset(pd_stimdf.cs,trial_time_stimOnT==x)
#                             this_stim_lm <- lmer(pupilCblz2~
#                                                    congruent.f + aSNR.f + congruent.f:aSNR.f +
#                                                    isH.f + blz + posX + posY + zptlen +
#                                                    (1 + congruent_aSNR.f1 +
#                                                       blz + posX + posY + zptlen|dataID),
#                                                  data=this_data,
#                                                  control=lmerControl(optimizer="bobyqa",
#                                                                      optCtrl=list(maxfun=2e5)))
#                           }
#                         },.progress="text")
#   names(stim.lm.cs.4) <- names(conv_all.zc.cs.4)
#   stim.lm.cs.4 <- stim.lm.cs.4[!conv_all.zc.cs.4]
#   save(stim.lm.cs.4,file=paste0(DATA_OUT_PATH,'ptO_stim_lm_cs4_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_stim_lm_cs4_2021-06-03.rda'))
# }
# 
# conv_all.cs.4 <- get_conv(stim.lm.cs.4)

#noconv
# if (T) {
#   stim.lm.cs.5 <- llply(as.integer(names(conv_all.zc.cs.5)),
#                         function(x) {
#                           if (!conv_all.zc.cs.5[as.character(x)]) {
#                             print(x)
#                             this_data <- subset(pd_stimdf.cs,trial_time_stimOnT==x)
#                             this_stim_lm <- lmer(pupilCblz2~
#                                                    congruent.f + aSNR.f + congruent.f:aSNR.f +
#                                                    isH.f + blz + posX + posY + zptlen +
#                                                    (1 +
#                                                       blz + posX + posY + zptlen|dataID),
#                                                  data=this_data,
#                                                  control=lmerControl(optimizer="bobyqa",
#                                                                      optCtrl=list(maxfun=2e5)))
#                           }
#                         },.progress="text")
#   names(stim.lm.cs.5) <- names(conv_all.zc.cs.5)
#   stim.lm.cs.5 <- stim.lm.cs.5[!conv_all.zc.cs.5]
#   save(stim.lm.cs.5,file=paste0(DATA_OUT_PATH,'ptO_stim_lm_cs5_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_stim_lm_cs5_2021-06-09.rda'))
# }
# 
# conv_all.cs.5 <- get_conv(stim.lm.cs.5)

#now the painful merge
stim.lm.cs.final <- stim.lm.zc.cs

if (all(names(stim.lm.zc.cs)[conv_all.zc.cs]==names(stim.lm.zc.cs.2)) & 
    all(names(stim.lm.zc.cs.2)[conv_all.zc.cs.2]==names(stim.lm.zc.cs.3)) &
    all(names(stim.lm.zc.cs.3)[conv_all.zc.cs.3]==names(stim.lm.zc.cs.4)) &
    all(names(stim.lm.zc.cs.4)[conv_all.zc.cs.4]==names(stim.lm.zc.cs.5))) {
  
  stim.lm.cs.final[conv_all.zc.cs] <- stim.lm.zc.cs.2
  stim.lm.cs.final[names(stim.lm.zc.cs.2)[conv_all.zc.cs.2]] <- stim.lm.zc.cs.3
  stim.lm.cs.final[names(stim.lm.zc.cs.3)[conv_all.zc.cs.3]] <- stim.lm.zc.cs.4
  stim.lm.cs.final[names(stim.lm.zc.cs.4)[conv_all.zc.cs.4]] <- stim.lm.zc.cs.5
} else {
  error('Data mismatch between full and reduced models! Cannot combine.')
}

#final check (not a conclusive check but it's something)
if (!all(sapply(stim.lm.cs.final,function(x) nrow(model.matrix(x))) == 
         sapply(stim.lm.zc.cs,function(x) nrow(model.matrix(x)))) | 
    any(get_conv(stim.lm.cs.final))) {
  stop('Data mismatch between full and reduced models!')
}


#extract coefficients/contrasts
stim.lm.cs.final.sum <- sum_to_df(stim.lm.cs.final,id='time',.pcorr=F,.confint=T,id2num=T)

stim.emm.cs.final <- llply(stim.lm.cs.final,
                              function(x) {
                                emm <- emmeans(x,specs=
                                                 list('congruent.f*aSNR.f'=~congruent.f*aSNR.f,
                                                      'congruent.f'=~congruent.f,
                                                      'aSNR.f'=~aSNR.f,
                                                      'isH.f'=~isH.f))
                                emm.pairs <- c(
                                  list('congruent.f*aSNR.f'= contrast(emm[[1]],method="revpairwise",adjust="none")[c(1,6)],
                                       'congruent.f'= contrast(emm[[2]],method="revpairwise",adjust="none")),
                                  pairs(emm,adjust="none",which=c(3,4)))
                                
                                
                              },.progress="text")

stim.emm.cs.final.sum <- emm_to_df(stim.emm.cs.final,id="time",.confint=T,id2num=T)

#OUTPUT
write.csv(stim.lm.cs.final.sum,paste0("stim_ptO_lm_cs_final_",Sys.Date(),".csv"),row.names=F)
write.csv(stim.emm.cs.final.sum,paste0("stim_ptO_emm_cs_final_",Sys.Date(),".csv"),row.names=F)