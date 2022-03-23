#clear memory
rm(list=ls())

#Runs stimulus-aligned regressions for the precue (rule-based) condition pupil data.


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
pd_file= paste0(DATA_OUT_PATH,'pupil_data_ds50_forR_09-Feb-2021.csv')
baseline_file = paste0(DATA_OUT_PATH,'pupil_data_bl_forR_09-Feb-2021.csv')

#source my functions if any
source('../sum_to_df.R')
source('../emm_to_df.R')

#set temporal parameters [this is easier than pupilnet cause I kept it in actual ms]
stim_period_ms.c <- c(60, 1360)
stim_period_ms.i <- c(60,860)
conType = "congruentS" #"congruent"

#load data
pd_df<-read.table(pd_file,sep=',', header=TRUE, stringsAsFactors=FALSE,na.strings = c('NaN'))
bl_df<-read.table(baseline_file,sep=',', header=TRUE, stringsAsFactors=FALSE,na.strings = c('NaN'))

head(pd_df)
head(bl_df)

#cleanup
bl_df$GroupCount <- NULL
pd_df[,c('time','pupilC','posXC','posYC')] <- NULL

#set up variables/factors
my_simple2<-contr.treatment(2,base=2) - matrix(rep(1/2,2))

pd_df$congruent.f <- factor(pd_df[,conType],levels=c(1,0,-1),
                            labels=c("congruent","incongruent","no prior"))
contrasts(pd_df$congruent.f) <- contr.sum(3)
print(conType)
contrasts(pd_df$congruent.f)

pd_df$isH.fs <- factor(pd_df$isH,levels=c(1,0),
                      labels=c("high","low"))
contrasts(pd_df$isH.fs) <- my_simple2
contrasts(pd_df$isH.fs)
pd_df$isH.f <- pd_df$isH.fs
contrasts(pd_df$isH.f) <- contr.sum(2)
contrasts(pd_df$isH.f)

pd_df$aSNR <- abs(pd_df$SNR)
#pd_df$aSNR.L <- round(poly(pd_df$aSNR,degree=1),10)

#merge in baseline data
pd_df <- left_join(pd_df,bl_df[,c('dataID','trialN','pupilBL')],by=c('dataID','trialN'))

#get rid of missing before scaling vars (probably already done w/ new preproc)
pd_df <- pd_df[complete.cases(pd_df),]

#get a handle on descriptives to inform scaling
scalevars_stats <- summarise(pd_df,
                             m_bl=mean(pupilBL,na.rm=T),
                             sd_bl=sd(pupilBL,na.rm=T),
                             m_x=mean(posXCbl,na.rm=T),
                             sd_x=sd(posXCbl,na.rm=T),
                             m_y=mean(posYCbl,na.rm=T),
                             sd_y=sd(posYCbl,na.rm=T),
                             m_asnr=mean(aSNR,na.rm=T),
                             sd_asnr=sd(aSNR,na.rm=T),
)

#set up variables: centering/scaling
pd_df$zaSNR <- scale(pd_df$aSNR)
pd_df$blz <- scale(pd_df$pupilBL,center=T,scale=F) #not scaling since already scaled
pd_df$posX <- scale(pd_df$posXCbl)
pd_df$posY <- scale(pd_df$posYCbl)


### STIM-ALIGNED REGRESSION ###

## CORRECT TRIALS ##
pd_stimdf.c <- subset(pd_df,success==1 & 
                   trial_time_stimOn >= stim_period_ms.c[1] & 
                   trial_time_stimOn <= stim_period_ms.c[2])


#fit full model (will likely be singular everywhere)
stim_times.c <- sort(unique(pd_stimdf.c$trial_time_stimOn))

if (F) {
  stim.lm <- llply(stim_times.c,
                   function(x) {
                     print(x)
                     this_data <- subset(pd_stimdf.c,trial_time_stimOn==x)
                     this_stim_lm <- lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY + 
                                          (1 + congruent.f + isH.f + zaSNR + blz + posX + posY|dataID),
                                        data=this_data,
                                        control=lmerControl(optimizer="bobyqa",
                                                            optCtrl=list(maxfun=2e5)))
                   },.progress="text")
  names(stim.lm) <- stim_times.c
  save(stim.lm,file=paste0(DATA_OUT_PATH,'stim_lm_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'stim_lm_2021-02-12.rda'))
}

conv_warning <- lapply(stim.lm,FUN=function(x) x@optinfo$conv$lme4$messages)
conv_warning2 <- lapply(stim.lm,FUN=function(x) x@optinfo$conv$warnings)
conv_sing <- sapply(stim.lm,FUN=isSingular)
conv_all <- !sapply(conv_warning,is.null) | !sapply(conv_warning2,is.null) | conv_sing

#do PCA of rfx and extract loadings/variance explained for each time point
rePCA_all <- lapply(stim.lm,rePCA)
pca_rot_bad <- lapply(rePCA_all[conv_all],FUN=function(x) x$dataID$rotation)
#pca_rot_bad <- array(unlist(pca_rot_bad),dim=c(8,8,length(pca_rot_bad)))
pca_sd_bad <- sapply(rePCA_all[conv_all],FUN=function(x) x$dataID$sdev)
pca_sd_bad_all <- rowMeans(pca_sd_bad)


#which rfx components most often load on the least-explanatory PCs? 
#(this is a little janky since averaging)
pca_rot_bad_maxs <- sapply(pca_rot_bad,function(x) apply(x,2,function(y) which.max(abs(y))))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

pca_rot_bad_modes <- apply(pca_rot_bad_maxs,1,Mode)

rfxnames <- names(ranef(stim.lm[[1]])$dataID)
worst_rfx <- rfxnames[pca_rot_bad_modes[pca_sd_bad_all < .1]]

#first let's try zero-corr model
pd_stimdf.c[,c('congruent.f1','congruent.f2')] <- 
  model.matrix(~1+pd_stimdf.c$congruent.f,pd_stimdf.c)[,2:3]
pd_stimdf.c$isH.f1 <- 
  model.matrix(~1+pd_stimdf.c$isH.f,pd_stimdf.c)[,2]
if (F) {
  stim.lm.zc <- llply(stim_times.c,
     function(x) {
       if (conv_all[as.character(x)]) {
         print(x)
         this_data <- subset(pd_stimdf.c,trial_time_stimOn==x)
         this_stim_lm <- 
           lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY + 
                            (1 + congruent.f1 + congruent.f2 + isH.f1 + zaSNR + blz + posX + posY||dataID), 
                          data=this_data, 
                          control=lmerControl(optimizer="bobyqa", 
                                              optCtrl=list(maxfun=2e5)))
       }
                     },.progress="text")
  names(stim.lm.zc) <- stim_times.c
  stim.lm.zc <- stim.lm.zc[conv_all]
  save(stim.lm.zc,file=paste0(DATA_OUT_PATH,'stim_lm_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'stim_lm_zc_2021-02-13.rda'))
}

conv_warning.zc <- lapply(stim.lm.zc,FUN=function(x) x@optinfo$conv$lme4$messages)
conv_warning2.zc <- lapply(stim.lm.zc,FUN=function(x) x@optinfo$conv$warnings)
conv_sing.zc <- sapply(stim.lm.zc,FUN=isSingular)
conv_all.zc <- !sapply(conv_warning.zc,is.null) | !sapply(conv_warning2.zc,is.null) | conv_sing.zc


#based on rePCA, eliminate congruent.f2, reintroduce corr
if (F) {
  stim.lm.2 <- llply(as.integer(names(conv_all.zc)), #stim_times.c[conv_all],
       function(x) {
         if (conv_all.zc[as.character(x)]) {
           print(x)
           this_data <- subset(pd_stimdf.c,trial_time_stimOn==x)
           this_stim_lm <-
             lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                    (1 + congruent.f1 + isH.f + zaSNR + blz + posX + posY|dataID),
                  data=this_data,
                  control=lmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))
         }
       },.progress="text")
  names(stim.lm.2) <- names(conv_all.zc) #stim_times.c[conv_all]
  stim.lm.2 <- stim.lm.2[conv_all.zc]
  save(stim.lm.2,file=paste0(DATA_OUT_PATH,'stim_lm_2_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'stim_lm_2_2021-02-13.rda'))
}

conv_warning.2 <- lapply(stim.lm.2,FUN=function(x) x@optinfo$conv$lme4$messages)
conv_warning2.2 <- lapply(stim.lm.2,FUN=function(x) x@optinfo$conv$warnings)
conv_sing.2 <- sapply(stim.lm.2,FUN=isSingular)
conv_all.2 <- !sapply(conv_warning.2,is.null) | !sapply(conv_warning2.2,is.null) | conv_sing.2


#re-suppress corr for failed above
if (F) {
  stim.lm.2.zc <- llply(as.integer(names(conv_all.2)),
        function(x) {
          if (conv_all.2[as.character(x)]) {
            print(x)
            this_data <- subset(pd_stimdf.c,trial_time_stimOn==x)
            this_stim_lm <-
              lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                     (1 + congruent.f1 + isH.f1 + zaSNR + blz + posX + posY||dataID),
                   data=this_data,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e5)))
          }
        },.progress="text")
  names(stim.lm.2.zc) <- names(conv_all.2)
  stim.lm.2.zc <- stim.lm.2.zc[conv_all.2]
  save(stim.lm.2.zc,file=paste0(DATA_OUT_PATH,'stim_lm_2_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'stim_lm_2_zc_2021-02-14.rda'))
}

conv_warning.2.zc <- lapply(stim.lm.2.zc,FUN=function(x) x@optinfo$conv$lme4$messages)
conv_warning2.2.zc <- lapply(stim.lm.2.zc,FUN=function(x) x@optinfo$conv$warnings)
conv_sing.2.zc <- sapply(stim.lm.2.zc,FUN=isSingular)
conv_all.2.zc <- !sapply(conv_warning.2.zc,is.null) | !sapply(conv_warning2.2.zc,is.null) | conv_sing.2.zc

#reintro-corr, remove isH slope
if (F) {
  stim.lm.3 <- llply(as.integer(names(conv_all.2.zc)),
                     function(x) {
                       if (conv_all.2.zc[as.character(x)]) {
                         print(x)
                         this_data <- subset(pd_stimdf.c,trial_time_stimOn==x)
                         this_stim_lm <-
                           lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                                  (1 + congruent.f1 + zaSNR + blz + posX + posY|dataID),
                                data=this_data,
                                control=lmerControl(optimizer="bobyqa",
                                                    optCtrl=list(maxfun=2e5)))
                       }
                     },.progress="text")
  names(stim.lm.3) <- names(conv_all.2.zc)
  stim.lm.3 <- stim.lm.3[conv_all.2.zc]
  save(stim.lm.3,file=paste0(DATA_OUT_PATH,'stim_lm_3_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'stim_lm_3_2021-02-14.rda'))
}

conv_warning.3 <- lapply(stim.lm.3,FUN=function(x) x@optinfo$conv$lme4$messages)
conv_warning2.3 <- lapply(stim.lm.3,FUN=function(x) x@optinfo$conv$warnings)
conv_sing.3 <- sapply(stim.lm.3,FUN=isSingular)
conv_all.3 <- !sapply(conv_warning.3,is.null) | !sapply(conv_warning2.3,is.null) | conv_sing.3

#remove corr
if (F) {
  stim.lm.3.zc <- llply(as.integer(names(conv_all.3)),
                        function(x) {
                          if (conv_all.3[as.character(x)]) {
                            print(x)
                            this_data <- subset(pd_stimdf.c,trial_time_stimOn==x)
                            this_stim_lm <-
                              lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                                     (1 + congruent.f1 + zaSNR + blz + posX + posY||dataID),
                                   data=this_data,
                                   control=lmerControl(optimizer="bobyqa",
                                                       optCtrl=list(maxfun=2e5)))
                          }
                        },.progress="text")
  names(stim.lm.3.zc) <- names(conv_all.3)
  stim.lm.3.zc <- stim.lm.3.zc[conv_all.3]
  save(stim.lm.3.zc,file=paste0(DATA_OUT_PATH,'stim_lm_3_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'stim_lm_3_zc_2021-02-14.rda'))
}

conv_warning.3.zc <- lapply(stim.lm.3.zc,FUN=function(x) x@optinfo$conv$lme4$messages)
conv_warning2.3.zc <- lapply(stim.lm.3.zc,FUN=function(x) x@optinfo$conv$warnings)
conv_sing.3.zc <- sapply(stim.lm.3.zc,FUN=isSingular)
conv_all.3.zc <- !sapply(conv_warning.3.zc,is.null) | !sapply(conv_warning2.3.zc,is.null) | conv_sing.3.zc

#reintro corr, remove zaSNR rfx based on spot check
if (F) {
  stim.lm.4 <- llply(as.integer(names(conv_all.3.zc)),
                     function(x) {
                       if (conv_all.3.zc[as.character(x)]) {
                         print(x)
                         this_data <- subset(pd_stimdf.c,trial_time_stimOn==x)
                         this_stim_lm <-
                           lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                                  (1 + congruent.f1 + blz + posX + posY|dataID),
                                data=this_data,
                                control=lmerControl(optimizer="bobyqa",
                                                    optCtrl=list(maxfun=2e5)))
                       }
                     },.progress="text")
  names(stim.lm.4) <- names(conv_all.3.zc)
  stim.lm.4 <- stim.lm.4[conv_all.3.zc]
  save(stim.lm.4,file=paste0(DATA_OUT_PATH,'stim_lm_4_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'stim_lm_4_2021-02-14.rda'))
}

conv_warning.4 <- lapply(stim.lm.4,FUN=function(x) x@optinfo$conv$lme4$messages)
conv_warning2.4 <- lapply(stim.lm.4,FUN=function(x) x@optinfo$conv$warnings)
conv_sing.4 <- sapply(stim.lm.4,FUN=isSingular)
conv_all.4 <- !sapply(conv_warning.4,is.null) | !sapply(conv_warning2.4,is.null) | conv_sing.4

#remove corr
if (T) {
  stim.lm.4.zc <- llply(as.integer(names(conv_all.4)),
                        function(x) {
                          if (conv_all.4[as.character(x)]) {
                            print(x)
                            this_data <- subset(pd_stimdf.c,trial_time_stimOn==x)
                            this_stim_lm <-
                              lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                                     (1 + congruent.f1 + blz + posX + posY||dataID),
                                   data=this_data,
                                   control=lmerControl(optimizer="bobyqa",
                                                       optCtrl=list(maxfun=2e5)))
                          }
                        },.progress="text")
  names(stim.lm.4.zc) <- names(conv_all.4)
  stim.lm.4.zc <- stim.lm.4.zc[conv_all.4]
  save(stim.lm.4.zc,file=paste0(DATA_OUT_PATH,'stim_lm_4_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'stim_lm_4_zc_2021-02-14.rda'))
}

conv_warning.4.zc <- lapply(stim.lm.4.zc,FUN=function(x) x@optinfo$conv$lme4$messages)
conv_warning2.4.zc <- lapply(stim.lm.4.zc,FUN=function(x) x@optinfo$conv$warnings)
conv_sing.4.zc <- sapply(stim.lm.4.zc,FUN=isSingular)
conv_all.4.zc <- !sapply(conv_warning.4.zc,is.null) | !sapply(conv_warning2.4.zc,is.null) | conv_sing.4.zc

any(conv_all.4.zc)


#combine models
stim.lm.final <- stim.lm

#being paranoid about order here, but looks OK.
if (all(names(stim.lm[conv_all])==names(stim.lm.zc)) & 
    all(names(stim.lm.zc[conv_all.zc])==names(stim.lm.2)) & 
    all(names(stim.lm.2[conv_all.2])==names(stim.lm.2.zc)) & 
    all(names(stim.lm.2.zc[conv_all.2.zc])==names(stim.lm.3)) &
    all(names(stim.lm.3[conv_all.3])==names(stim.lm.3.zc)) &
    all(names(stim.lm.3.zc[conv_all.3.zc])==names(stim.lm.4)) &
    all(names(stim.lm.4[conv_all.4])==names(stim.lm.4.zc))
    ) {
  stim.lm.final[conv_all] <- stim.lm.zc
  stim.lm.final[names(stim.lm.zc[conv_all.zc])] <- stim.lm.2
  stim.lm.final[names(stim.lm.2[conv_all.2])] <- stim.lm.2.zc
  stim.lm.final[names(stim.lm.2.zc[conv_all.2.zc])] <- stim.lm.3
  stim.lm.final[names(stim.lm.3[conv_all.3])] <- stim.lm.3.zc
  stim.lm.final[names(stim.lm.3.zc[conv_all.3.zc])] <- stim.lm.4
  stim.lm.final[names(stim.lm.4[conv_all.4])] <- stim.lm.4.zc
} else {
  error('Data mismatch between full and reduced models! Cannot combine.')
}


#final check (not a conclusive check but it's something)
if (!all(sapply(stim.lm.final,function(x) nrow(model.matrix(x))) == 
         sapply(stim.lm,function(x) nrow(model.matrix(x))))) {
  error('Data mismatch between full and reduced models!')
}
  
source('../sum_to_df.R')
#extract summary tables
stim.lm.final.sum <- sum_to_df(stim.lm.final,id='time',.pcorr=T,.confint=T,id2num=T)

#extract ANOVA tables
stim.lm.final.aov <- sum_to_df(stim.lm.final,id='time',.pcorr=T,id2num=T,type="anova",anova_type="II")

stim.emm.final <- llply(stim.lm.final,
                          function(x) {
                            emm <- emmeans(x,specs=list('congruent.f'=~congruent.f,'isH.f'=~isH.f))
                            emm.pairs <- pairs(emm,adjust="none")
                            
                          },.progress="text")
stim.emm.final.sum <- emm_to_df(stim.emm.final,id="time",.confint=T,id2num=T)


#outout regressions for plotting in matlab
write.csv(stim.lm.final.sum,paste0("stim_lm_final_",Sys.Date(),".csv"),row.names=F)
write.csv(stim.emm.final.sum,paste0("stim_emm_final_",Sys.Date(),".csv"),row.names=F)

