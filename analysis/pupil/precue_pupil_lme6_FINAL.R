#clear memory
rm(list=ls())

### this runs regressions for priorOnly choice-aligned pupil data


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
source('../get_conv.R')

#set temporal parameters [this is easier than pupilnet cause I kept it in actual ms]
choice_period_ms.c <- c(-250, 1000)
choice_period_ms.c.wide <- c(-980, 1000)
choice_period_ms.i <- c(-250,500)
conType = "congruent" #"congruent"
widetime = T

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


### RESPONSE REGRESSION ###

## CORRECT TRIALS ##
if (widetime) {
  pd_choicedf.c <- subset(pd_df,success==1 & 
                     trial_time_choice >= choice_period_ms.c.wide[1] & 
                     trial_time_choice <= choice_period_ms.c.wide[2])
} else {
  pd_choicedf.c <- subset(pd_df,success==1 & 
                            trial_time_choice >= choice_period_ms.c[1] & 
                            trial_time_choice <= choice_period_ms.c[2])
}

#get choice times
choice_times.c <- sort(unique(pd_choicedf.c$trial_time_choice))


#first let's try zero-corr model
pd_choicedf.c[,c('congruent.f1','congruent.f2')] <- 
  model.matrix(~1+pd_choicedf.c$congruent.f,pd_choicedf.c)[,2:3]
pd_choicedf.c$isH.f1 <- 
  model.matrix(~1+pd_choicedf.c$isH.f,pd_choicedf.c)[,2]
if (F) {
  choice.lm.zc <- llply(choice_times.c,
     function(x) {
         print(x)
         this_data <- subset(pd_choicedf.c,trial_time_choice==x)
         this_choice_lm <- 
           lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY + 
                            (1 + congruent.f1 + congruent.f2 + isH.f1 + zaSNR + blz + posX + posY||dataID), 
                          data=this_data, 
                          control=lmerControl(optimizer="bobyqa", 
                                              optCtrl=list(maxfun=2e5)))
                     },.progress="text")
  names(choice.lm.zc) <- choice_times.c
  save(choice.lm.zc,file=paste0(DATA_OUT_PATH,'choice_lm_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_zc_2021-05-18.rda'))
}

conv_all.zc <- get_conv(choice.lm.zc)

#do PCA of rfx and extract loadings/variance explained for each time point
rePCA_all <- lapply(choice.lm.zc,rePCA)
pca_rot_bad <- lapply(rePCA_all[conv_all.zc],FUN=function(x) x$dataID$rotation)
pca_sd_bad <- sapply(rePCA_all[conv_all.zc],FUN=function(x) x$dataID$sdev)
pca_sd_bad_all <- rowMeans(pca_sd_bad)


#which rfx components most often load on the least-explanatory PCs? 
#(this is a little janky since averaging)
pca_rot_bad_maxs <- sapply(pca_rot_bad,function(x) apply(x,2,function(y) which.max(abs(y))))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

pca_rot_bad_modes <- apply(pca_rot_bad_maxs,1,Mode)

rfxnames <- names(ranef(choice.lm.zc[[1]])$dataID)
worst_rfx <- rfxnames[pca_rot_bad_modes[pca_sd_bad_all < .1]]

#for any zc that did converge nonsingularly, reintro corr
if (F) {
  choice.lm <- llply(choice_times.c,
                 function(x) {
                   if (!conv_all.zc[as.character(x)]) {
                     print(x)
                     this_data <- subset(pd_choicedf.c,trial_time_choice==x)
                     this_choice_lm <- lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                                          (1 + congruent.f + isH.f + zaSNR + blz + posX + posY|dataID),
                                        data=this_data,
                                        control=lmerControl(optimizer="bobyqa",
                                                            optCtrl=list(maxfun=2e5)))
                   }
                 },.progress="text")
  names(choice.lm) <- choice_times.c
  choice.lm <- choice.lm.zc[!conv_all.zc]
  save(choice.lm,file=paste0(DATA_OUT_PATH,'choice_lm_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_2021-05-18.rda'))
}

conv_all <- get_conv(choice.lm)

#based on rePCA, eliminate congruent.f2, reintroduce corr
if (F) {
  choice.lm.2 <- llply(as.integer(names(conv_all.zc)), 
       function(x) {
         if (conv_all.zc[as.character(x)]) {
           print(x)
           this_data <- subset(pd_choicedf.c,trial_time_choice==x)
           this_choice_lm <-
             lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                    (1 + congruent.f1 + isH.f + zaSNR + blz + posX + posY|dataID),
                  data=this_data,
                  control=lmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))
         }
       },.progress="text")
  names(choice.lm.2) <- names(conv_all.zc) 
  choice.lm.2 <- choice.lm.2[conv_all.zc]
  save(choice.lm.2,file=paste0(DATA_OUT_PATH,'choice_lm_2_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_2_2021-05-19.rda'))
}

conv_all.2 <- get_conv(choice.lm.2)


#re-suppress corr for failed above
if (F) {
  choice.lm.2.zc <- llply(as.integer(names(conv_all.2)),
        function(x) {
          if (conv_all.2[as.character(x)]) {
            print(x)
            this_data <- subset(pd_choicedf.c,trial_time_choice==x)
            this_choice_lm <-
              lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                     (1 + congruent.f1 + isH.f1 + zaSNR + blz + posX + posY||dataID),
                   data=this_data,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e5)))
          }
        },.progress="text")
  names(choice.lm.2.zc) <- names(conv_all.2)
  choice.lm.2.zc <- choice.lm.2.zc[conv_all.2]
  save(choice.lm.2.zc,file=paste0(DATA_OUT_PATH,'choice_lm_2_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_2_zc_2021-05-19.rda'))
}

conv_all.2.zc <- get_conv(choice.lm.2.zc)

#do PCA of rfx and extract loadings/variance explained for each time point
rePCA_all.2 <- lapply(choice.lm.2.zc,rePCA)
pca_rot_bad.2 <- lapply(rePCA_all.2[conv_all.2.zc],FUN=function(x) x$dataID$rotation)
pca_sd_bad.2 <- sapply(rePCA_all.2[conv_all.2.zc],FUN=function(x) x$dataID$sdev)
pca_sd_bad_all.2 <- rowMeans(pca_sd_bad.2)


#which rfx components most often load on the least-explanatory PCs? 
#(this is a little janky since averaging)
pca_rot_bad_maxs.2 <- sapply(pca_rot_bad.2,function(x) apply(x,2,function(y) which.max(abs(y))))

pca_rot_bad_modes.2 <- apply(pca_rot_bad_maxs.2,1,Mode)

rfxnames.2 <- names(ranef(choice.lm.2.zc[[1]])$dataID)
worst_rfx.2 <- rfxnames.2[pca_rot_bad_modes.2[pca_sd_bad_all.2 < .1]]


#NEXT remove isH
if (F) {
  choice.lm.3.zc <- llply(as.integer(names(conv_all.2.zc)),
                          function(x) {
                            if (conv_all.2.zc[as.character(x)]) {
                              print(x)
                              this_data <- subset(pd_choicedf.c,trial_time_choice==x)
                              this_choice_lm <-
                                lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                                       (1 + congruent.f1 + zaSNR + blz + posX + posY||dataID),
                                     data=this_data,
                                     control=lmerControl(optimizer="bobyqa",
                                                         optCtrl=list(maxfun=2e5)))
                            }
                          },.progress="text")
  names(choice.lm.3.zc) <- names(conv_all.2.zc)
  choice.lm.3.zc <- choice.lm.3.zc[conv_all.2.zc]
  save(choice.lm.3.zc,file=paste0(DATA_OUT_PATH,'choice_lm_3_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_3_zc_2021-05-19.rda'))
}

conv_all.3.zc <- get_conv(choice.lm.3.zc)

#remove zaSNR [in spot check reintro corr failed, so not doing so]
if (F) {
  choice.lm.4.zc <- llply(as.integer(names(conv_all.3.zc)),
                          function(x) {
                            if (conv_all.3.zc[as.character(x)]) {
                              print(x)
                              this_data <- subset(pd_choicedf.c,trial_time_choice==x)
                              this_choice_lm <-
                                lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                                       (1 + congruent.f1 + blz + posX + posY||dataID),
                                     data=this_data,
                                     control=lmerControl(optimizer="bobyqa",
                                                         optCtrl=list(maxfun=2e5)))
                            }
                          },.progress="text")
  names(choice.lm.4.zc) <- names(conv_all.3.zc)
  choice.lm.4.zc <- choice.lm.4.zc[conv_all.3.zc]
  save(choice.lm.4.zc,file=paste0(DATA_OUT_PATH,'choice_lm_4_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_4_zc_2021-05-19.rda'))
}

conv_all.4.zc <- get_conv(choice.lm.4.zc)

#remove congruent.f1 [still can't reintro corr--just one time point!]
if (F) {
  choice.lm.5.zc <- llply(as.integer(names(conv_all.4.zc)),
                          function(x) {
                            if (conv_all.4.zc[as.character(x)]) {
                              print(x)
                              this_data <- subset(pd_choicedf.c,trial_time_choice==x)
                              this_choice_lm <-
                                lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                                       (1 + blz + posX + posY||dataID),
                                     data=this_data,
                                     control=lmerControl(optimizer="bobyqa",
                                                         optCtrl=list(maxfun=2e5)))
                            }
                          },.progress="text")
  names(choice.lm.5.zc) <- names(conv_all.4.zc)
  choice.lm.5.zc <- choice.lm.5.zc[conv_all.4.zc]
  save(choice.lm.5.zc,file=paste0(DATA_OUT_PATH,'choice_lm_5_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_5_zc_2021-05-19.rda'))
}

conv_all.5.zc <- get_conv(choice.lm.5.zc)
any(conv_all.5.zc)

#reintro corr for those that fit [NONE NONSING]
# if (T) {
#   choice.lm.3 <- llply(as.integer(names(conv_all.3.zc)),
#                        function(x) {
#                          if (!conv_all.3.zc[as.character(x)]) {
#                            print(x)
#                            this_data <- subset(pd_choicedf.c,trial_time_choice==x)
#                            this_choice_lm <-
#                              lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
#                                     (1 + congruent.f1 + zaSNR + blz + posX + posY|dataID),
#                                   data=this_data,
#                                   control=lmerControl(optimizer="bobyqa",
#                                                       optCtrl=list(maxfun=2e5)))
#                          }
#                        },.progress="text")
#   names(choice.lm.3) <- names(conv_all.3.zc)
#   choice.lm.3 <- choice.lm.3[!conv_all.3.zc]
#   save(choice.lm.3,file=paste0(DATA_OUT_PATH,'choice_lm_3_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'choice_lm_3_2021-05-19.rda'))
# }
# 
# conv_all.3 <- get_conv(choice.lm.3)


#package this up now
choice.lm.final <- choice.lm.zc

#being paranoid about order here, but looks OK.
if (all(names(choice.lm.zc[conv_all.zc])==names(choice.lm.2)) & 
  all(names(choice.lm.zc[!conv_all.zc])==names(choice.lm)) & 
  !any(conv_all) & #need this here to prevent case of adding unconverged/singular choice.lm
  all(names(choice.lm.2[conv_all.2])==names(choice.lm.2.zc)) & 
  all(names(choice.lm.2.zc[conv_all.2.zc])==names(choice.lm.3.zc)) & 
  all(names(choice.lm.3.zc[conv_all.3.zc])==names(choice.lm.4.zc)) & 
  all(names(choice.lm.4.zc[conv_all.4.zc])==names(choice.lm.5.zc))) {
    
  choice.lm.final[conv_all.zc] <- choice.lm.2
  choice.lm.final[!conv_all.zc] <- choice.lm
  choice.lm.final[names(choice.lm.2[conv_all.2])] <- choice.lm.2.zc
  choice.lm.final[names(choice.lm.2.zc[conv_all.2.zc])] <- choice.lm.3.zc
  choice.lm.final[names(choice.lm.3.zc[conv_all.3.zc])] <- choice.lm.4.zc
  choice.lm.final[names(choice.lm.4.zc[conv_all.4.zc])] <- choice.lm.5.zc

  } else {
    error('Data mismatch between full and reduced models! Cannot combine.')
  }


#final check (not a conclusive check but it's something)
if (!all(sapply(choice.lm.final,function(x) nrow(model.matrix(x))) == 
         sapply(choice.lm.zc,function(x) nrow(model.matrix(x)))) | 
    any(get_conv(choice.lm.final))) {
  error('Data mismatch between full and reduced models!')
}
  
#extract summary tables
choice.lm.final.sum <- sum_to_df(choice.lm.final,id='time',.pcorr=T,.confint=T,id2num=T)

#extract ANOVA tables
choice.lm.final.aov <- sum_to_df(choice.lm.final,id='time',.pcorr=T,id2num=T,type="anova",anova_type="II")

choice.emm.final <- llply(choice.lm.final,
                          function(x) {
                            emm <- emmeans(x,specs=list('congruent.f'=~congruent.f,'isH.f'=~isH.f))
                            emm.pairs <- pairs(emm,adjust="none")
                            
                          },.progress="text")
choice.emm.final.sum <- emm_to_df(choice.emm.final,id="time",.confint=T,id2num=T)


#save output for matlab notebook
write.csv(choice.lm.final.sum,paste0("choice_lm_final_",Sys.Date(),".csv"),row.names=F)
write.csv(choice.emm.final.sum,paste0("choice_emm_final_",Sys.Date(),".csv"),row.names=F)




### CORRECT TRIALS congruentxSNR interaction check (extreme conds only) ###

if (conType=="congruent") {
  pd_choicedf.cs <- subset(pd_choicedf.c,congruent %in% c(0,1) & aSNR %in% c(.05,.5))
} else {
  pd_choicedf.cs <- subset(pd_choicedf.c,congruentS %in% c(0,1) & aSNR %in% c(.05,.5))
}
unique(pd_choicedf.cs$aSNR)
unique(pd_choicedf.cs[,conType])

pd_choicedf.cs$congruent.f <- factor(pd_choicedf.cs[,conType],levels=c(1,0),
                                     labels=c("con","incon"))
contrasts(pd_choicedf.cs$congruent.f) <- contr.sum(2)
contrasts(pd_choicedf.cs$congruent.f)
unique(pd_choicedf.cs[,c(conType,'congruent.f')])

pd_choicedf.cs$aSNR.f <- factor(pd_choicedf.cs$aSNR,levels=c(.5,.05),
                                labels=c("high SNR","low SNR"))
contrasts(pd_choicedf.cs$aSNR.f) <- contr.sum(2)
contrasts(pd_choicedf.cs$aSNR.f)
unique(pd_choicedf.cs[,c('aSNR','aSNR.f')])

## FIRST: fit no-corr model, reduced dataset to see what we've got

pd_choicedf.cs[,c('congruent.f1','aSNR.f1','congruent_aSNR.f1')] <- 
  model.matrix(~1+pd_choicedf.cs$congruent.f*pd_choicedf.cs$aSNR.f,pd_choicedf.cs)[,2:4]
pd_choicedf.cs$isH.f1 <- 
  model.matrix(~1+pd_choicedf.cs$isH.f,pd_choicedf.cs)[,2]
unique(pd_choicedf.cs[,c("congruent.f","aSNR.f","congruent.f1","aSNR.f1","congruent_aSNR.f1")])
unique(pd_choicedf.cs[,c("isH.f","isH.f1")])

##ZC
if (F) {
  choice.lm.zc.cs <- llply(choice_times.c,
                           function(x) {
                             print(x)
                             this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                             this_choice_lm <- lmer(pupilCblz~
                                                      congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                      isH.f + blz + posX + posY + 
                                                      (1 + congruent.f1 + aSNR.f1 + congruent_aSNR.f1 + 
                                                         isH.f1 + blz + posX + posY||dataID),
                                                    data=this_data,
                                                    control=lmerControl(optimizer="bobyqa",
                                                                        optCtrl=list(maxfun=2e5)))
                           },.progress="text")
  names(choice.lm.zc.cs) <- choice_times.c
  save(choice.lm.zc.cs,file=paste0(DATA_OUT_PATH,'choice_lm_zc_cs_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_zc_cs_2021-06-09.rda'))
}

conv_all.zc.cs <- get_conv(choice.lm.zc.cs)

#do PCA of rfx and extract loadings/variance explained for each time point
#could just extract variance components since zero correlation but easier to use existing code
rePCA_all <- lapply(choice.lm.zc.cs,rePCA)
pca_rot_bad <- lapply(rePCA_all[conv_all.zc.cs],FUN=function(x) x$dataID$rotation)
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

rfxnames <- names(ranef(choice.lm.zc.cs[[1]])$dataID)
worst_rfx <- rfxnames[pca_rot_bad_modes[pca_sd_bad_all < .1]]

#based on above, first try eliminating congruent.f rfx and refit
if (T) {
  choice.lm.zc.cs.2 <- llply(choice_times.c,
                             function(x) {
                               if (conv_all.zc.cs[as.character(x)]) {
                                 print(x)
                                 this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                                 this_choice_lm <- lmer(pupilCblz~
                                                          congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                          isH.f + blz + posX + posY + 
                                                          (1 + aSNR.f1 + congruent_aSNR.f1 + 
                                                             isH.f1 + blz + posX + posY||dataID),
                                                        data=this_data,
                                                        control=lmerControl(optimizer="bobyqa",
                                                                            optCtrl=list(maxfun=2e5)))
                               }
                             },.progress="text")
  names(choice.lm.zc.cs.2) <- choice_times.c
  choice.lm.zc.cs.2 <- choice.lm.zc.cs.2[conv_all.zc.cs]
  save(choice.lm.zc.cs.2,file=paste0(DATA_OUT_PATH,'choice_lm_zc_cs2_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_zc_cs2_2021-06-10.rda'))
}

conv_all.zc.cs.2 <- get_conv(choice.lm.zc.cs.2)


#do PCA of rfx and extract loadings/variance explained for each time point
#could just extract variance components since zero correlation but easier to use existing code
rePCA_all.2 <- lapply(choice.lm.zc.cs.2,rePCA)
pca_rot_bad.2 <- lapply(rePCA_all.2[conv_all.zc.cs.2],FUN=function(x) x$dataID$rotation)
pca_sd_bad.2 <- sapply(rePCA_all.2[conv_all.zc.cs.2],FUN=function(x) x$dataID$sdev)
pca_sd_bad.2_all <- rowMeans(pca_sd_bad.2)


#which rfx components most often load on the least-explanatory PCs? 
#(this is a little janky since averaging)
pca_rot_bad.2_maxs <- sapply(pca_rot_bad.2,function(x) apply(x,2,function(y) which.max(abs(y))))

pca_rot_bad.2_modes <- apply(pca_rot_bad.2_maxs,1,Mode)

rfxnames.2 <- names(ranef(choice.lm.zc.cs.2[[1]])$dataID)
worst_rfx.2 <- rfxnames[pca_rot_bad.2_modes[pca_sd_bad.2_all < .1]]

#remove aSNR
if (F) {
  choice.lm.zc.cs.3 <- llply(as.integer(names(conv_all.zc.cs.2)),
                             function(x) {
                               if (conv_all.zc.cs.2[as.character(x)]) {
                                 print(x)
                                 this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                                 this_choice_lm <- lmer(pupilCblz~
                                                          congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                          isH.f + blz + posX + posY + 
                                                          (1 + congruent_aSNR.f1 + 
                                                             isH.f1 + blz + posX + posY||dataID),
                                                        data=this_data,
                                                        control=lmerControl(optimizer="bobyqa",
                                                                            optCtrl=list(maxfun=2e5)))
                               }
                             },.progress="text")
  names(choice.lm.zc.cs.3) <- names(conv_all.zc.cs.2)
  choice.lm.zc.cs.3 <- choice.lm.zc.cs.3[conv_all.zc.cs.2]
  save(choice.lm.zc.cs.3,file=paste0(DATA_OUT_PATH,'choice_lm_zc_cs3_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_zc_cs3_2021-06-10.rda'))
}

conv_all.zc.cs.3 <- get_conv(choice.lm.zc.cs.3)

#remove isH
if (T) {
  choice.lm.zc.cs.4 <- llply(as.integer(names(conv_all.zc.cs.3)),
                             function(x) {
                               if (conv_all.zc.cs.3[as.character(x)]) {
                                 print(x)
                                 this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                                 this_choice_lm <- lmer(pupilCblz~
                                                          congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                          isH.f + blz + posX + posY + 
                                                          (1 + congruent_aSNR.f1 + 
                                                             blz + posX + posY||dataID),
                                                        data=this_data,
                                                        control=lmerControl(optimizer="bobyqa",
                                                                            optCtrl=list(maxfun=2e5)))
                               }
                             },.progress="text")
  names(choice.lm.zc.cs.4) <- names(conv_all.zc.cs.3)
  choice.lm.zc.cs.4 <- choice.lm.zc.cs.4[conv_all.zc.cs.3]
  save(choice.lm.zc.cs.4,file=paste0(DATA_OUT_PATH,'choice_lm_zc_cs4_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_zc_cs4_2021-06-10.rda'))
}

conv_all.zc.cs.4 <- get_conv(choice.lm.zc.cs.4)

#remove interaction (and blz for -980)
if (F) {
  choice.lm.zc.cs.5 <- llply(as.integer(names(conv_all.zc.cs.4)),
                             function(x) {
                               if (conv_all.zc.cs.4[as.character(x)]) {
                                 print(x)
                                 this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                                 if (x==-980) {
                                   this_choice_lm <- lmer(pupilCblz~
                                                            congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                            isH.f + blz + posX + posY + 
                                                            (1 + 
                                                               posX + posY||dataID),
                                                          data=this_data,
                                                          control=lmerControl(optimizer="bobyqa",
                                                                              optCtrl=list(maxfun=2e5)))
                                 } else {
                                   this_choice_lm <- lmer(pupilCblz~
                                                            congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                            isH.f + blz + posX + posY + 
                                                            (1 + 
                                                               blz + posX + posY||dataID),
                                                          data=this_data,
                                                          control=lmerControl(optimizer="bobyqa",
                                                                              optCtrl=list(maxfun=2e5)))
                                 }
                               }
                             },.progress="text")
  names(choice.lm.zc.cs.5) <- names(conv_all.zc.cs.4)
  choice.lm.zc.cs.5 <- choice.lm.zc.cs.5[conv_all.zc.cs.4]
  save(choice.lm.zc.cs.5,file=paste0(DATA_OUT_PATH,'choice_lm_zc_cs5_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_zc_cs5_2021-06-10.rda'))
}

conv_all.zc.cs.5 <- get_conv(choice.lm.zc.cs.5)

#Now try to reintro corr

#no conv
# if (F) {
#   choice.lm.cs <- llply(as.integer(names(conv_all.zc.cs)),
#                         function(x) {
#                           if (!conv_all.zc.cs[as.character(x)]) {
#                             print(x)
#                             this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
#                             this_choice_lm <- lmer(pupilCblz~
#                                                      congruent.f + aSNR.f + congruent.f:aSNR.f + 
#                                                      isH.f + blz + posX + posY + 
#                                                      (1 + congruent.f1 + aSNR.f1 + congruent_aSNR.f1 + 
#                                                         isH.f1 + blz + posX + posY|dataID),
#                                                    data=this_data,
#                                                    control=lmerControl(optimizer="bobyqa",
#                                                                        optCtrl=list(maxfun=2e5)))
#                           }
#                         },.progress="text")
#   names(choice.lm.cs) <- names(conv_all.zc.cs)
#   choice.lm.cs <- choice.lm.cs[!conv_all.zc.cs]
#   save(choice.lm.cs,file=paste0(DATA_OUT_PATH,'choice_lm_cs_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'choice_lm_cs_2021-06-10.rda'))
# }
# 
# conv_all.cs <- get_conv(choice.lm.cs)


if (F) {
  choice.lm.cs.2 <- llply(as.integer(names(conv_all.zc.cs.2)),
                           function(x) {
                             if (!conv_all.zc.cs.2[as.character(x)]) {
                               print(x)
                               this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                               this_choice_lm <- lmer(pupilCblz~
                                                        congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                        isH.f + blz + posX + posY + 
                                                        (1 + aSNR.f1 + congruent_aSNR.f1 + 
                                                           isH.f1 + blz + posX + posY|dataID),
                                                      data=this_data,
                                                      control=lmerControl(optimizer="bobyqa",
                                                                          optCtrl=list(maxfun=2e5)))
                             }
                           },.progress="text")
  names(choice.lm.cs.2) <- names(conv_all.zc.cs.2)
  choice.lm.cs.2 <- choice.lm.cs.2[!conv_all.zc.cs.2]
  save(choice.lm.cs.2,file=paste0(DATA_OUT_PATH,'choice_lm_cs2_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_cs2_2021-06-10.rda'))
}

conv_all.cs.2 <- get_conv(choice.lm.cs.2)

if (F) {
  choice.lm.cs.3 <- llply(as.integer(names(conv_all.zc.cs.3)),
                          function(x) {
                            if (!conv_all.zc.cs.3[as.character(x)]) {
                              print(x)
                              this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                              this_choice_lm <- lmer(pupilCblz~
                                                       congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                       isH.f + blz + posX + posY + 
                                                       (1 + congruent_aSNR.f1 + 
                                                          isH.f1 + blz + posX + posY|dataID),
                                                     data=this_data,
                                                     control=lmerControl(optimizer="bobyqa",
                                                                         optCtrl=list(maxfun=2e5)))
                            }
                          },.progress="text")
  names(choice.lm.cs.3) <- names(conv_all.zc.cs.3)
  choice.lm.cs.3 <- choice.lm.cs.3[!conv_all.zc.cs.3]
  save(choice.lm.cs.3,file=paste0(DATA_OUT_PATH,'choice_lm_cs3_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_cs3_2021-06-10.rda'))
}

conv_all.cs.3 <- get_conv(choice.lm.cs.3)

if (F) {
  choice.lm.cs.4 <- llply(as.integer(names(conv_all.zc.cs.4)),
                          function(x) {
                            if (!conv_all.zc.cs.4[as.character(x)]) {
                              print(x)
                              this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                              this_choice_lm <- lmer(pupilCblz~
                                                       congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                       isH.f + blz + posX + posY + 
                                                       (1 + congruent_aSNR.f1 + 
                                                          blz + posX + posY|dataID),
                                                     data=this_data,
                                                     control=lmerControl(optimizer="bobyqa",
                                                                         optCtrl=list(maxfun=2e5)))
                            }
                          },.progress="text")
  names(choice.lm.cs.4) <- names(conv_all.zc.cs.4)
  choice.lm.cs.4 <- choice.lm.cs.4[!conv_all.zc.cs.4]
  save(choice.lm.cs.4,file=paste0(DATA_OUT_PATH,'choice_lm_cs4_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_lm_cs4_2021-06-10.rda'))
}

conv_all.cs.4 <- get_conv(choice.lm.cs.4)

#all sing
# if (F) {
#   choice.lm.cs.5 <- llply(as.integer(names(conv_all.zc.cs.5)),
#                           function(x) {
#                             if (!conv_all.zc.cs.5[as.character(x)]) {
#                               print(x)
#                               this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
#                               if (x==-980) {
#                                 this_choice_lm <- lmer(pupilCblz~
#                                                          congruent.f + aSNR.f + congruent.f:aSNR.f + 
#                                                          isH.f + blz + posX + posY + 
#                                                          (1 + 
#                                                             posX + posY|dataID),
#                                                        data=this_data,
#                                                        control=lmerControl(optimizer="bobyqa",
#                                                                            optCtrl=list(maxfun=2e5)))
#                               } else {
#                                 this_choice_lm <- lmer(pupilCblz~
#                                                          congruent.f + aSNR.f + congruent.f:aSNR.f + 
#                                                          isH.f + blz + posX + posY + 
#                                                          (1 + 
#                                                             blz + posX + posY|dataID),
#                                                        data=this_data,
#                                                        control=lmerControl(optimizer="bobyqa",
#                                                                            optCtrl=list(maxfun=2e5)))
#                               }
#                             }
#                           },.progress="text")
#   names(choice.lm.cs.5) <- names(conv_all.zc.cs.5)
#   choice.lm.cs.5 <- choice.lm.cs.5[!conv_all.zc.cs.5]
#   save(choice.lm.cs.5,file=paste0(DATA_OUT_PATH,'choice_lm_cs5_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'choice_lm_cs5_2021-06-10.rda'))
# }
# 
# conv_all.cs.5 <- get_conv(choice.lm.cs.5)


#now the merge
choice.lm.cs.final <- choice.lm.zc.cs

if (all(names(choice.lm.zc.cs)[conv_all.zc.cs]==names(choice.lm.zc.cs.2)) & 
    all(names(choice.lm.zc.cs.2)[conv_all.zc.cs.2]==names(choice.lm.zc.cs.3)) &
    all(names(choice.lm.zc.cs.3)[conv_all.zc.cs.3]==names(choice.lm.zc.cs.4)) &
    all(names(choice.lm.zc.cs.4)[conv_all.zc.cs.4]==names(choice.lm.zc.cs.5)) &
    all(names(choice.lm.zc.cs.2)[!conv_all.zc.cs.2]==names(choice.lm.cs.2)) &
    all(names(choice.lm.zc.cs.3)[!conv_all.zc.cs.3]==names(choice.lm.cs.3)) &
    all(names(choice.lm.zc.cs.4)[!conv_all.zc.cs.4]==names(choice.lm.cs.4))) {
  
  choice.lm.cs.final[conv_all.zc.cs] <- choice.lm.zc.cs.2
  choice.lm.cs.final[names(choice.lm.zc.cs.2)[conv_all.zc.cs.2]] <- choice.lm.zc.cs.3
  choice.lm.cs.final[names(choice.lm.zc.cs.3)[conv_all.zc.cs.3]] <- choice.lm.zc.cs.4
  choice.lm.cs.final[names(choice.lm.zc.cs.4)[conv_all.zc.cs.4]] <- choice.lm.zc.cs.5
  choice.lm.cs.final[names(choice.lm.cs.2)[!conv_all.cs.2]] <- choice.lm.cs.2[!conv_all.cs.2] #reintro corr case
  choice.lm.cs.final[names(choice.lm.cs.3)[!conv_all.cs.3]] <- choice.lm.cs.3[!conv_all.cs.3] #reintro corr case
  choice.lm.cs.final[names(choice.lm.cs.4)[!conv_all.cs.4]] <- choice.lm.cs.4[!conv_all.cs.4] #reintro corr case
  
} else {
  error('Data mismatch between full and reduced models! Cannot combine.')
}

#final check (not a conclusive check but it's something)
if (!all(sapply(choice.lm.cs.final,function(x) nrow(model.matrix(x))) == 
         sapply(choice.lm.zc.cs,function(x) nrow(model.matrix(x)))) | 
    any(get_conv(choice.lm.cs.final))) {
  stop('Data mismatch between full and reduced models!')
}

#extract coefficients/contrasts
choice.lm.cs.final.sum <- sum_to_df(choice.lm.cs.final,id='time',.pcorr=F,.confint=T,id2num=T)

choice.emm.cs.final <- llply(choice.lm.cs.final,
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

choice.emm.cs.final.sum <- emm_to_df(choice.emm.cs.final,id="time",.confint=T,id2num=T)

#OUTPUT
write.csv(choice.lm.cs.final.sum,paste0("choice_lm_cs_final_",Sys.Date(),".csv"),row.names=F)
write.csv(choice.emm.cs.final.sum,paste0("choice_emm_cs_final_",Sys.Date(),".csv"),row.names=F)


## INCORRECT TRIALS ##


pd_choicedf.i <- subset(pd_df,success==0 & 
                          trial_time_choice >= choice_period_ms.i[1] & 
                          trial_time_choice <= choice_period_ms.i[2])


unique(pd_choicedf.i[,c('prior','isH','choice01','congruent','congruent.f')])


#fit full model (will likely be singular everywhere)
choice_times.i <- sort(unique(pd_choicedf.i$trial_time_choice))

if (F) {
  choice.i.lm <- llply(choice_times.i,
                       function(x) {
                         print(x)
                         this_data <- subset(pd_choicedf.i,trial_time_choice==x)
                         this_choice_lm <- lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY + 
                                                  (1 + congruent.f + isH.f + zaSNR + blz + posX + posY|dataID),
                                                data=this_data,
                                                control=lmerControl(optimizer="bobyqa",
                                                                    optCtrl=list(maxfun=2e5)))
                       },.progress="text")
  names(choice.i.lm) <- choice_times.i
  save(choice.i.lm,file=paste0(DATA_OUT_PATH,'choice_i_lm_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_i_lm_2021-04-29.rda'))
}


conv.i_warning <- lapply(choice.i.lm,FUN=function(x) x@optinfo$conv.i$lme4$messages)
conv.i_warning2 <- lapply(choice.i.lm,FUN=function(x) x@optinfo$conv.i$warnings)
conv.i_sing <- sapply(choice.i.lm,FUN=isSingular)
conv.i_all <- !sapply(conv.i_warning,is.null) | !sapply(conv.i_warning2,is.null) | conv.i_sing

#do PCA of rfx and extract loadings/variance explained for each time point
rePCA_all.i <- lapply(choice.i.lm,rePCA)
pca_rot_bad.i <- lapply(rePCA_all.i[conv.i_all],FUN=function(x) x$dataID$rotation)
pca_sd_bad.i <- sapply(rePCA_all.i[conv.i_all],FUN=function(x) x$dataID$sdev)
pca_sd_bad_all.i <- rowMeans(pca_sd_bad.i)


#which rfx components most often load on the least-explanatory PCs? 
#(this is a little janky since averaging)
pca_rot_bad_maxs.i <- sapply(pca_rot_bad.i,function(x) apply(x,2,function(y) which.max(abs(y))))

pca_rot_bad_modes.i <- apply(pca_rot_bad_maxs.i,1,Mode)

rfxnames.i <- names(ranef(choice.i.lm[[1]])$dataID)
worst_rfx.i <- rfxnames[pca_rot_bad_modes.i[pca_sd_bad_all.i < .1]]

#in incorrect case, isH seems to be biggest offender, so will try dropping. But first try zc
#first let's try zero-corr model
pd_choicedf.i[,c('congruent.f1','congruent.f2')] <- 
  model.matrix(~1+pd_choicedf.i$congruent.f,pd_choicedf.i)[,2:3]
pd_choicedf.i$isH.f1 <- 
  model.matrix(~1+pd_choicedf.i$isH.f,pd_choicedf.i)[,2]
if (F) {
  choice.i.lm.zc <- llply(choice_times.i,
                          function(x) {
                            if (conv.i_all[as.character(x)]) {
                              print(x)
                              this_data <- subset(pd_choicedf.i,trial_time_choice==x)
                              this_choice_lm <- 
                                lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY + 
                                       (1 + congruent.f1 + congruent.f2 + isH.f1 + zaSNR + blz + posX + posY||dataID), 
                                     data=this_data, 
                                     control=lmerControl(optimizer="bobyqa", 
                                                         optCtrl=list(maxfun=2e5)))
                            }
                          },.progress="text")
  names(choice.i.lm.zc) <- choice_times.i
  choice.i.lm.zc <- choice.i.lm.zc[conv.i_all]
  save(choice.i.lm.zc,file=paste0(DATA_OUT_PATH,'choice_i_lm_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_i_lm_zc_2021-04-29.rda'))
}

conv.i_warning.zc <- lapply(choice.i.lm.zc,FUN=function(x) x@optinfo$conv.i$lme4$messages)
conv.i_warning2.zc <- lapply(choice.i.lm.zc,FUN=function(x) x@optinfo$conv.i$warnings)
conv.i_sing.zc <- sapply(choice.i.lm.zc,FUN=isSingular)
conv.i_all.zc <- !sapply(conv.i_warning.zc,is.null) | !sapply(conv.i_warning2.zc,is.null) | conv.i_sing.zc


#based on rePCA, eliminate isH, reintroduce corr
if (F) {
  choice.i.lm.2 <- llply(as.integer(names(conv.i_all.zc)), 
                         function(x) {
                           if (conv.i_all.zc[as.character(x)]) {
                             print(x)
                             this_data <- subset(pd_choicedf.i,trial_time_choice==x)
                             this_choice_lm <-
                               lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                                      (1 + congruent.f + zaSNR + blz + posX + posY|dataID),
                                    data=this_data,
                                    control=lmerControl(optimizer="bobyqa",
                                                        optCtrl=list(maxfun=2e5)))
                           }
                         },.progress="text")
  names(choice.i.lm.2) <- names(conv.i_all.zc)
  choice.i.lm.2 <- choice.i.lm.2[conv.i_all.zc]
  save(choice.i.lm.2,file=paste0(DATA_OUT_PATH,'choice_i_lm_2_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_i_lm_2_2021-04-29.rda'))
}

conv.i_warning.2 <- lapply(choice.i.lm.2,FUN=function(x) x@optinfo$conv.i$lme4$messages)
conv.i_warning2.2 <- lapply(choice.i.lm.2,FUN=function(x) x@optinfo$conv.i$warnings)
conv.i_sing.2 <- sapply(choice.i.lm.2,FUN=isSingular)
conv.i_all.2 <- !sapply(conv.i_warning.2,is.null) | !sapply(conv.i_warning2.2,is.null) | conv.i_sing.2

#re-suppress corr for failed above
if (F) {
  choice.i.lm.2.zc <- llply(as.integer(names(conv.i_all.2)),
                            function(x) {
                              if (conv.i_all.2[as.character(x)]) {
                                print(x)
                                this_data <- subset(pd_choicedf.i,trial_time_choice==x)
                                this_choice_lm <-
                                  lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                                         (1 + congruent.f1 + congruent.f2 + zaSNR + blz + posX + posY||dataID),
                                       data=this_data,
                                       control=lmerControl(optimizer="bobyqa",
                                                           optCtrl=list(maxfun=2e5)))
                              }
                            },.progress="text")
  names(choice.i.lm.2.zc) <- names(conv.i_all.2)
  choice.i.lm.2.zc <- choice.i.lm.2.zc[conv.i_all.2]
  save(choice.i.lm.2.zc,file=paste0(DATA_OUT_PATH,'choice_i_lm_2_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_i_lm_2_zc_2021-04-29.rda'))
}

conv.i_warning.2.zc <- lapply(choice.i.lm.2.zc,FUN=function(x) x@optinfo$conv$lme4$messages)
conv.i_warning2.2.zc <- lapply(choice.i.lm.2.zc,FUN=function(x) x@optinfo$conv$warnings)
conv.i_sing.2.zc <- sapply(choice.i.lm.2.zc,FUN=isSingular)
conv.i_all.2.zc <- !sapply(conv.i_warning.2.zc,is.null) | !sapply(conv.i_warning2.2.zc,is.null) | conv.i_sing.2.zc

#now remove congruent.f2
if (F) {
  choice.i.lm.3.zc <- llply(as.integer(names(conv.i_all.2.zc)),
                            function(x) {
                              if (conv.i_all.2.zc[as.character(x)]) {
                                print(x)
                                this_data <- subset(pd_choicedf.i,trial_time_choice==x)
                                this_choice_lm <-
                                  lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                                         (1 + congruent.f1 + zaSNR + blz + posX + posY||dataID),
                                       data=this_data,
                                       control=lmerControl(optimizer="bobyqa",
                                                           optCtrl=list(maxfun=2e5)))
                              }
                            },.progress="text")
  names(choice.i.lm.3.zc) <- names(conv.i_all.2.zc)
  choice.i.lm.3.zc <- choice.i.lm.3.zc[conv.i_all.2.zc]
  save(choice.i.lm.3.zc,file=paste0(DATA_OUT_PATH,'choice_i_lm_3_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_i_lm_3_zc_2021-04-29.rda'))
}

conv.i_warning.3.zc <- lapply(choice.i.lm.3.zc,FUN=function(x) x@optinfo$conv$lme4$messages)
conv.i_warning2.3.zc <- lapply(choice.i.lm.3.zc,FUN=function(x) x@optinfo$conv$warnings)
conv.i_sing.3.zc <- sapply(choice.i.lm.3.zc,FUN=isSingular)
conv.i_all.3.zc <- !sapply(conv.i_warning.3.zc,is.null) | !sapply(conv.i_warning2.3.zc,is.null) | conv.i_sing.3.zc

#based on spot checking, remove zaSNR, reintro congruent.f2
if (F) {
  choice.i.lm.4.zc <- llply(as.integer(names(conv.i_all.3.zc)),
                            function(x) {
                              if (conv.i_all.3.zc[as.character(x)]) {
                                print(x)
                                this_data <- subset(pd_choicedf.i,trial_time_choice==x)
                                this_choice_lm <-
                                  lmer(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY +
                                         (1 + congruent.f1 + congruent.f2 + blz + posX + posY||dataID),
                                       data=this_data,
                                       control=lmerControl(optimizer="bobyqa",
                                                           optCtrl=list(maxfun=2e5)))
                              }
                            },.progress="text")
  names(choice.i.lm.4.zc) <- names(conv.i_all.3.zc)
  choice.i.lm.4.zc <- choice.i.lm.4.zc[conv.i_all.3.zc]
  save(choice.i.lm.4.zc,file=paste0(DATA_OUT_PATH,'choice_i_lm_4_zc_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'choice_i_lm_4_zc_2021-04-29.rda'))
}

conv.i_warning.4.zc <- lapply(choice.i.lm.4.zc,FUN=function(x) x@optinfo$conv$lme4$messages)
conv.i_warning2.4.zc <- lapply(choice.i.lm.4.zc,FUN=function(x) x@optinfo$conv$warnings)
conv.i_sing.4.zc <- sapply(choice.i.lm.4.zc,FUN=isSingular)
conv.i_all.4.zc <- !sapply(conv.i_warning.4.zc,is.null) | !sapply(conv.i_warning2.4.zc,is.null) | conv.i_sing.4.zc

any(conv.i_all.4.zc)


#Construct final model set
choice.i.lm.final <- choice.i.lm

#being paranoid about order here, but looks OK.
if (all(names(choice.i.lm[conv.i_all])==names(choice.i.lm.zc)) & 
    all(names(choice.i.lm.zc[conv.i_all.zc])==names(choice.i.lm.2)) & 
    all(names(choice.i.lm.2[conv.i_all.2])==names(choice.i.lm.2.zc)) &
    all(names(choice.i.lm.2.zc[conv.i_all.2.zc])==names(choice.i.lm.3.zc)) &
    all(names(choice.i.lm.3.zc[conv.i_all.3.zc])==names(choice.i.lm.4.zc))) {
  choice.i.lm.final[conv.i_all] <- choice.i.lm.zc
  choice.i.lm.final[names(choice.i.lm.zc[conv.i_all.zc])] <- choice.i.lm.2
  choice.i.lm.final[names(choice.i.lm.2[conv.i_all.2])] <- choice.i.lm.2.zc
  choice.i.lm.final[names(choice.i.lm.2.zc[conv.i_all.2.zc])] <- choice.i.lm.3.zc
  choice.i.lm.final[names(choice.i.lm.3.zc[conv.i_all.3.zc])] <- choice.i.lm.4.zc
} else {
  error('Data mismatch between full and reduced models! Cannot combine.')
}

#final check (not a conclusive check but it's something)
if (!all(sapply(choice.i.lm.final,function(x) nrow(model.matrix(x))) == 
         sapply(choice.i.lm,function(x) nrow(model.matrix(x))))) {
  error('Data mismatch between full and reduced models!')
}

#extract summary tables
choice.i.lm.final.sum <- sum_to_df(choice.i.lm.final,id='time',.pcorr=T,.confint=T,id2num=T)

#extract simple effects
choice.i.emm.final <- llply(choice.i.lm.final,
                            function(x) {
                              emm <- emmeans(x,specs=list('congruent.f'=~congruent.f,'isH.f'=~isH.f))
                              emm.pairs <- pairs(emm,adjust="none")
                              
                            },.progress="text")
choice.i.emm.final.sum <- emm_to_df(choice.i.emm.final,id="time",.confint=T,id2num=T)


#save
write.csv(choice.i.lm.final.sum,paste0("choice_i_lm_final_",Sys.Date(),".csv"),row.names=F)
write.csv(choice.i.emm.final.sum,paste0("choice_i_emm_final_",Sys.Date(),".csv"),row.names=F)