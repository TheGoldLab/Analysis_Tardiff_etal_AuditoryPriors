#clear memory
rm(list=ls())

### this runs regressions for pretoneOnly choice-aligned pupil data

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
choice_period_ms.early <- c(-900,0)
choice_period_ms.c.wide <- c(-980, 1000)

widetime = T
widetime.i = F
conType = "congruent" #"congruent"

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


### CORRECT TRIALS ####

if (widetime) {
  pd_choicedf.c <- subset(pd_df,success==1 & 
                             trial_time_choice >= choice_period_ms.c.wide[1] & 
                             trial_time_choice <= choice_period_ms.c.wide[2])
} else {
  pd_choicedf.c <- subset(pd_df,success==1 & 
                            trial_time_choice >= choice_period_ms.early[1] & 
                            trial_time_choice <= choice_period_ms.early[2])
}

if (conType=="congruent") {
  pd_choicedf.cs <- subset(pd_choicedf.c,congruent %in% c(0,3) & aSNR %in% c(.05,.5))
} else {
  pd_choicedf.cs <- subset(pd_choicedf.c,congruentS %in% c(0,3) & aSNR %in% c(.05,.5))
}
unique(pd_choicedf.cs$aSNR)
unique(pd_choicedf.cs[,conType])

pd_choicedf.cs$congruent.f <- factor(pd_choicedf.cs[,conType],levels=c(3,0),
                            labels=c("con","incon"))
contrasts(pd_choicedf.cs$congruent.f) <- contr.sum(2)
contrasts(pd_choicedf.cs$congruent.f)
unique(pd_choicedf.cs[,c(conType,'congruent.f')])

pd_choicedf.cs$aSNR.f <- factor(pd_choicedf.cs$aSNR,levels=c(.5,.05),
                                 labels=c("high SNR","low SNR"))
contrasts(pd_choicedf.cs$aSNR.f) <- contr.sum(2)
contrasts(pd_choicedf.cs$aSNR.f)
unique(pd_choicedf.cs[,c('aSNR','aSNR.f')])


choice_times.c <- sort(unique(pd_choicedf.c$trial_time_choice))

### REGRESSIONS ###

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
                       this_choice_lm <- lmer(pupilCblz2~
                                        congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                          isH.f + blz + posX + posY + zptlen + 
                                        (1 + congruent.f1 + aSNR.f1 + congruent_aSNR.f1 + 
                                           isH.f1 + blz + posX + posY + zptlen||dataID),
                                      data=this_data,
                                      control=lmerControl(optimizer="bobyqa",
                                                          optCtrl=list(maxfun=2e5)))
                     },.progress="text")
  names(choice.lm.zc.cs) <- choice_times.c
  save(choice.lm.zc.cs,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_cs_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_cs_2021-06-09.rda'))
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
if (F) {
  choice.lm.zc.cs.2 <- llply(choice_times.c,
                    function(x) {
                      if (conv_all.zc.cs[as.character(x)]) {
                        print(x)
                        this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                        this_choice_lm <- lmer(pupilCblz2~
                                                 congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                 isH.f + blz + posX + posY + zptlen + 
                                                 (1 + aSNR.f1 + congruent_aSNR.f1 + 
                                                    isH.f1 + blz + posX + posY + zptlen||dataID),
                                               data=this_data,
                                               control=lmerControl(optimizer="bobyqa",
                                                                   optCtrl=list(maxfun=2e5)))
                      }
                    },.progress="text")
  names(choice.lm.zc.cs.2) <- choice_times.c
  choice.lm.zc.cs.2 <- choice.lm.zc.cs.2[conv_all.zc.cs]
  save(choice.lm.zc.cs.2,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_cs2_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_cs2_2021-06-09.rda'))
}

conv_all.zc.cs.2 <- get_conv(choice.lm.zc.cs.2)

#remove isH
if (F) {
  choice.lm.zc.cs.3 <- llply(as.integer(names(conv_all.zc.cs.2)),
                              function(x) {
                                if (conv_all.zc.cs.2[as.character(x)]) {
                                  print(x)
                                  this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                                  this_choice_lm <- lmer(pupilCblz2~
                                                           congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                           isH.f + blz + posX + posY + zptlen + 
                                                           (1 + aSNR.f1 + congruent_aSNR.f1 + 
                                                              blz + posX + posY + zptlen||dataID),
                                                         data=this_data,
                                                         control=lmerControl(optimizer="bobyqa",
                                                                             optCtrl=list(maxfun=2e5)))
                                }
                              },.progress="text")
  names(choice.lm.zc.cs.3) <- names(conv_all.zc.cs.2)
  choice.lm.zc.cs.3 <- choice.lm.zc.cs.3[conv_all.zc.cs.2]
  save(choice.lm.zc.cs.3,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_cs3_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_cs3_2021-06-09.rda'))
}

conv_all.zc.cs.3 <- get_conv(choice.lm.zc.cs.3)

#now remove interaction
if (F) {
  choice.lm.zc.cs.4 <- llply(as.integer(names(conv_all.zc.cs.3)),
                              function(x) {
                                if (conv_all.zc.cs.3[as.character(x)]) {
                                  print(x)
                                  this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                                  this_choice_lm <- lmer(pupilCblz2~
                                                           congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                           isH.f + blz + posX + posY + zptlen + 
                                                           (1 + aSNR.f1 +
                                                              blz + posX + posY + zptlen||dataID),
                                                         data=this_data,
                                                         control=lmerControl(optimizer="bobyqa",
                                                                             optCtrl=list(maxfun=2e5)))
                                }
                              },.progress="text")
  names(choice.lm.zc.cs.4) <- names(conv_all.zc.cs.3)
  choice.lm.zc.cs.4 <- choice.lm.zc.cs.4[conv_all.zc.cs.3]
  save(choice.lm.zc.cs.4,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_cs4_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_cs4_2021-06-09.rda'))
}

conv_all.zc.cs.4 <- get_conv(choice.lm.zc.cs.4)

summary(choice.lm.zc.cs.4$`-980`)$varcor

#spot check suggests in early cases need to drop zptlen and aSNR
if (F) {
  choice.lm.zc.cs.5 <- llply(as.integer(names(conv_all.zc.cs.4)),
                              function(x) {
                                if (conv_all.zc.cs.4[as.character(x)]) {
                                  print(x)
                                  this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                                  this_choice_lm <- lmer(pupilCblz2~
                                                           congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                           isH.f + blz + posX + posY + zptlen + 
                                                           (1 + 
                                                              blz + posX + posY||dataID),
                                                         data=this_data,
                                                         control=lmerControl(optimizer="bobyqa",
                                                                             optCtrl=list(maxfun=2e5)))
                                }
                              },.progress="text")
  names(choice.lm.zc.cs.5) <- names(conv_all.zc.cs.4)
  choice.lm.zc.cs.5 <- choice.lm.zc.cs.5[conv_all.zc.cs.4]
  save(choice.lm.zc.cs.5,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_cs5_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_cs5_2021-06-09.rda'))
}

conv_all.zc.cs.5 <- get_conv(choice.lm.zc.cs.5)

any(conv_all.zc.cs.5)


## REINTRO CORR IN MODELS THAT FIT
# #ALL SIGNGULAR
# if (T) {
#   choice.lm.cs <- llply(choice_times.c,
#                          function(x) {
#                            if (!conv_all.zc.cs[as.character(x)]) {
#                              print(x)
#                              this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
#                              this_choice_lm <- lmer(pupilCblz2~
#                                                       congruent.f + aSNR.f + congruent.f:aSNR.f +
#                                                       isH.f + blz + posX + posY + zptlen +
#                                                       (1 + congruent.f1 + aSNR.f1 + congruent_aSNR.f1 +
#                                                          isH.f1 + blz + posX + posY + zptlen|dataID),
#                                                     data=this_data,
#                                                     control=lmerControl(optimizer="bobyqa",
#                                                                         optCtrl=list(maxfun=2e5)))
#                            }
#                          },.progress="text")
#   names(choice.lm.cs) <- choice_times.c
#   choice.lm.cs <- choice.lm.cs[!conv_all.zc.cs]
#   save(choice.lm.cs,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_cs_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_cs_2021-06-09.rda'))
# }

#all sing
# if (T) {
#   choice.lm.cs.2 <- llply(as.integer(names(conv_all.zc.cs.2)),
#                            function(x) {
#                              if (!conv_all.zc.cs.2[as.character(x)]) {
#                                print(x)
#                                this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
#                                this_choice_lm <- lmer(pupilCblz2~
#                                                         congruent.f + aSNR.f + congruent.f:aSNR.f +
#                                                         isH.f + blz + posX + posY + zptlen +
#                                                         (1 + aSNR.f1 + congruent_aSNR.f1 +
#                                                            isH.f1 + blz + posX + posY + zptlen|dataID),
#                                                       data=this_data,
#                                                       control=lmerControl(optimizer="bobyqa",
#                                                                           optCtrl=list(maxfun=2e5)))
#                              }
#                            },.progress="text")
#   names(choice.lm.cs.2) <- names(conv_all.zc.cs.2)
#   choice.lm.cs.2 <- choice.lm.cs.2[!conv_all.zc.cs.2]
#   save(choice.lm.cs.2,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_cs2_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_cs2_2021-06-09.rda'))
# }

#nada
# if (T) {
#   choice.lm.cs.3 <- llply(as.integer(names(conv_all.zc.cs.3)),
#                            function(x) {
#                              if (!conv_all.zc.cs.3[as.character(x)]) {
#                                print(x)
#                                this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
#                                 this_choice_lm <- lmer(pupilCblz2~
#                                                          congruent.f + aSNR.f + congruent.f:aSNR.f +
#                                                          isH.f + blz + posX + posY + zptlen +
#                                                          (1 + aSNR.f1 + congruent_aSNR.f1 +
#                                                             blz + posX + posY + zptlen|dataID),
#                                                        data=this_data,
#                                                        control=lmerControl(optimizer="bobyqa",
#                                                                            optCtrl=list(maxfun=2e5)))
#                              }
#                            },.progress="text")
#   names(choice.lm.cs.3) <- names(conv_all.zc.cs.3)
#   choice.lm.cs.3 <- choice.lm.cs.3[!conv_all.zc.cs.3]
#   save(choice.lm.cs.3,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_cs3_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_cs3_2021-06-09.rda'))
# }

#conv_all.cs.3 <- get_conv(choice.lm.cs.3)


# if (F) {
#   choice.lm.cs.4 <- llply(as.integer(names(conv_all.zc.cs.4)),
#                            function(x) {
#                              if (!conv_all.zc.cs.4[as.character(x)]) {
#                                print(x)
#                                this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
#                                this_choice_lm <- lmer(pupilCblz2~
#                                                         congruent.f + aSNR.f + congruent.f:aSNR.f +
#                                                         isH.f + blz + posX + posY + zptlen +
#                                                         (1 + aSNR.f1 +
#                                                            blz + posX + posY + zptlen|dataID),
#                                                       data=this_data,
#                                                       control=lmerControl(optimizer="bobyqa",
#                                                                           optCtrl=list(maxfun=2e5)))
#                              }
#                            },.progress="text")
#   names(choice.lm.cs.4) <- names(conv_all.zc.cs.4)
#   choice.lm.cs.4 <- choice.lm.cs.4[!conv_all.zc.cs.4]
#   save(choice.lm.cs.4,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_cs4_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_cs4_2021-06-09.rda'))
# }
# 
# conv_all.cs.4 <- get_conv(choice.lm.cs.4)

if (F) {
  choice.lm.cs.5 <- llply(as.integer(names(conv_all.zc.cs.5)),
                           function(x) {
                             if (!conv_all.zc.cs.5[as.character(x)]) {
                               print(x)
                               this_data <- subset(pd_choicedf.cs,trial_time_choice==x)
                               this_choice_lm <- lmer(pupilCblz2~
                                                        congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                        isH.f + blz + posX + posY + zptlen + 
                                                        (1 + 
                                                           blz + posX + posY|dataID),
                                                      data=this_data,
                                                      control=lmerControl(optimizer="bobyqa",
                                                                          optCtrl=list(maxfun=2e5)))
                             }
                           },.progress="text")
  names(choice.lm.cs.5) <- names(conv_all.zc.cs.5)
  choice.lm.cs.5 <- choice.lm.cs.5[!conv_all.zc.cs.5]
  save(choice.lm.cs.5,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_cs5_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_cs5_2021-06-09.rda'))
}

#OK a few of these converged....
conv_all.cs.5 <- get_conv(choice.lm.cs.5)


#now the merge
choice.lm.cs.final <- choice.lm.zc.cs

if (all(names(choice.lm.zc.cs)[conv_all.zc.cs]==names(choice.lm.zc.cs.2)) & 
    all(names(choice.lm.zc.cs.2)[conv_all.zc.cs.2]==names(choice.lm.zc.cs.3)) &
    all(names(choice.lm.zc.cs.3)[conv_all.zc.cs.3]==names(choice.lm.zc.cs.4)) &
    all(names(choice.lm.zc.cs.4)[conv_all.zc.cs.4]==names(choice.lm.zc.cs.5)) &
    all(names(choice.lm.zc.cs.5)[!conv_all.zc.cs.5]==names(choice.lm.cs.5))) {
  
  choice.lm.cs.final[conv_all.zc.cs] <- choice.lm.zc.cs.2
  choice.lm.cs.final[names(choice.lm.zc.cs.2)[conv_all.zc.cs.2]] <- choice.lm.zc.cs.3
  choice.lm.cs.final[names(choice.lm.zc.cs.3)[conv_all.zc.cs.3]] <- choice.lm.zc.cs.4
  choice.lm.cs.final[names(choice.lm.zc.cs.4)[conv_all.zc.cs.4]] <- choice.lm.zc.cs.5
  choice.lm.cs.final[names(choice.lm.cs.5)[!conv_all.cs.5]] <- choice.lm.cs.5[!conv_all.cs.5] #reintro corr case

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
write.csv(choice.lm.cs.final.sum,paste0("choice_ptO_lm_cs_final_",Sys.Date(),".csv"),row.names=F)
write.csv(choice.emm.cs.final.sum,paste0("choice_ptO_emm_cs_final_",Sys.Date(),".csv"),row.names=F)



### INCORRECT TRIALS ###
if (widetime.i) {
  pd_choicedf.i <- subset(pd_df,success==0 & 
                            trial_time_choice >= choice_period_ms.c.wide[1] & 
                            trial_time_choice <= choice_period_ms.c.wide[2])
} else {
  pd_choicedf.i <- subset(pd_df,success==0 & 
                            trial_time_choice >= choice_period_ms.early[1] & 
                            trial_time_choice <= choice_period_ms.early[2])
}

pd_choicedf.is <- subset(pd_choicedf.i,congruent %in% c(0,3) & aSNR %in% c(.05,.5))
unique(pd_choicedf.is$aSNR)
unique(pd_choicedf.is$congruent)

#arguably some of this stuff should have been done before splitting into correct/incorrect, but OK...
pd_choicedf.is$congruent.f <- factor(pd_choicedf.is[,conType],levels=c(3,0),
                                     labels=c("con","incon"))
contrasts(pd_choicedf.is$congruent.f) <- contr.sum(2)
contrasts(pd_choicedf.is$congruent.f)
unique(pd_choicedf.is[,c(conType,'congruent.f')])

pd_choicedf.is$aSNR.f <- factor(pd_choicedf.is$aSNR,levels=c(.5,.05),
                                labels=c("high SNR","low SNR"))
contrasts(pd_choicedf.is$aSNR.f) <- contr.sum(2)
contrasts(pd_choicedf.is$aSNR.f)
unique(pd_choicedf.is[,c('aSNR','aSNR.f')])


choice_times.i <- sort(unique(pd_choicedf.i$trial_time_choice))


## FIRST: fit no-corr model, reduced dataset to see what we've got

pd_choicedf.is[,c('congruent.f1','aSNR.f1','congruent_aSNR.f1')] <- 
  model.matrix(~1+pd_choicedf.is$congruent.f*pd_choicedf.is$aSNR.f,pd_choicedf.is)[,2:4]
pd_choicedf.is$isH.f1 <- 
  model.matrix(~1+pd_choicedf.is$isH.f,pd_choicedf.is)[,2]
unique(pd_choicedf.is[,c("congruent.f","aSNR.f","congruent.f1","aSNR.f1","congruent_aSNR.f1")])
unique(pd_choicedf.is[,c("isH.f","isH.f1")])

##ZC
if (F) {
  choice.lm.zc.is <- llply(choice_times.i,
                           function(x) {
                             print(x)
                             this_data <- subset(pd_choicedf.is,trial_time_choice==x)
                             this_choice_lm <- lmer(pupilCblz2~
                                                      congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                      isH.f + blz + posX + posY + zptlen + 
                                                      (1 + congruent.f1 + aSNR.f1 + congruent_aSNR.f1 + 
                                                         isH.f1 + blz + posX + posY + zptlen||dataID),
                                                    data=this_data,
                                                    control=lmerControl(optimizer="bobyqa",
                                                                        optCtrl=list(maxfun=2e5)))
                           },.progress="text")
  names(choice.lm.zc.is) <- choice_times.i
  save(choice.lm.zc.is,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_is_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_is_2021-06-09.rda'))
}

conv_all.zc.is <- get_conv(choice.lm.zc.is)

#NO convergence on this puppy.
#do PCA of rfx and extract loadings/variance explained for each time point
#could just extract variance components since zero correlation but easier to use existing code
rePCA_all.i <- lapply(choice.lm.zc.is,rePCA)
pca_rot_bad.i <- lapply(rePCA_all.i[conv_all.zc.is],FUN=function(x) x$dataID$rotation)
pca_sd_bad.i <- sapply(rePCA_all.i[conv_all.zc.is],FUN=function(x) x$dataID$sdev)
pca_sd_bad_all.i <- rowMeans(pca_sd_bad.i)


#which rfx components most often load on the least-explanatory PCs? 
#(this is a little janky since averaging)
pca_rot_bad_maxs.i <- sapply(pca_rot_bad.i,function(x) apply(x,2,function(y) which.max(abs(y))))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

pca_rot_bad_modes.i <- apply(pca_rot_bad_maxs.i,1,Mode)

rfxnames.i <- names(ranef(choice.lm.zc.is[[1]])$dataID)
worst_rfx.i <- rfxnames.i[pca_rot_bad_modes.i[pca_sd_bad_all.i < .1]]


#based on above, remove congruent and refit
if (F) {
  choice.lm.zc.is.2 <- llply(choice_times.i,
                           function(x) {
                             if (conv_all.zc.is) {
                               print(x)
                               this_data <- subset(pd_choicedf.is,trial_time_choice==x)
                               this_choice_lm <- lmer(pupilCblz2~
                                                        congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                        isH.f + blz + posX + posY + zptlen + 
                                                        (1 + aSNR.f1 + congruent_aSNR.f1 + 
                                                           isH.f1 + blz + posX + posY + zptlen||dataID),
                                                      data=this_data,
                                                      control=lmerControl(optimizer="bobyqa",
                                                                          optCtrl=list(maxfun=2e5)))
                             }
                           },.progress="text")
  names(choice.lm.zc.is.2) <- choice_times.i
  choice.lm.zc.is.2 <- choice.lm.zc.is.2[conv_all.zc.is]
  save(choice.lm.zc.is.2,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_is2_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_is2_2021-06-10.rda'))
}

conv_all.zc.is.2 <- get_conv(choice.lm.zc.is.2)

#NOW REMOVE aSNR!
if (F) {
  choice.lm.zc.is.3 <- llply(as.integer(names(choice.lm.zc.is.2)),
                             function(x) {
                               if (conv_all.zc.is.2[as.character(x)]) {
                                 print(x)
                                 this_data <- subset(pd_choicedf.is,trial_time_choice==x)
                                 this_choice_lm <- lmer(pupilCblz2~
                                                          congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                          isH.f + blz + posX + posY + zptlen + 
                                                          (1 + congruent_aSNR.f1 + 
                                                             isH.f1 + blz + posX + posY + zptlen||dataID),
                                                        data=this_data,
                                                        control=lmerControl(optimizer="bobyqa",
                                                                            optCtrl=list(maxfun=2e5)))
                               }
                             },.progress="text")
  names(choice.lm.zc.is.3) <- names(conv_all.zc.is.2)
  choice.lm.zc.is.3 <- choice.lm.zc.is.3[conv_all.zc.is.2]
  save(choice.lm.zc.is.3,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_is3_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_is3_2021-06-10.rda'))
}

conv_all.zc.is.3 <- get_conv(choice.lm.zc.is.3)

rePCA_all.i.2 <- lapply(choice.lm.zc.is.3,rePCA)
pca_rot_bad.i.2 <- lapply(rePCA_all.i[conv_all.zc.is.3],FUN=function(x) x$dataID$rotation)
pca_sd_bad.i.2 <- sapply(rePCA_all.i[conv_all.zc.is.3],FUN=function(x) x$dataID$sdev)
pca_sd_bad_all.i.2 <- rowMeans(pca_sd_bad.i.2)


#which rfx components most often load on the least-explanatory PCs? 
#(this is a little janky since averaging)
pca_rot_bad_maxs.i.2 <- sapply(pca_rot_bad.i.2,function(x) apply(x,2,function(y) which.max(abs(y))))

pca_rot_bad_modes.i.2 <- apply(pca_rot_bad_maxs.i.2,1,Mode)

rfxnames.i.2 <- names(ranef(choice.lm.zc.is.3[[1]])$dataID)
worst_rfx.i.2 <- rfxnames.i[pca_rot_bad_modes.i.2[pca_sd_bad_all.i.2 < .1]]

#now remove isH
if (T) {
  choice.lm.zc.is.4 <- llply(as.integer(names(choice.lm.zc.is.3)),
                             function(x) {
                               if (conv_all.zc.is.3[as.character(x)]) {
                                 print(x)
                                 this_data <- subset(pd_choicedf.is,trial_time_choice==x)
                                 this_choice_lm <- lmer(pupilCblz2~
                                                          congruent.f + aSNR.f + congruent.f:aSNR.f + 
                                                          isH.f + blz + posX + posY + zptlen + 
                                                          (1 + congruent_aSNR.f1 + 
                                                             blz + posX + posY + zptlen||dataID),
                                                        data=this_data,
                                                        control=lmerControl(optimizer="bobyqa",
                                                                            optCtrl=list(maxfun=2e5)))
                               }
                             },.progress="text")
  names(choice.lm.zc.is.4) <- names(conv_all.zc.is.3)
  choice.lm.zc.is.4 <- choice.lm.zc.is.4[conv_all.zc.is.3]
  save(choice.lm.zc.is.4,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_is4_', Sys.Date(), '.rda'))
} else {
  load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_zc_is4_2021-06-10.rda'))
}

conv_all.zc.is.4 <- get_conv(choice.lm.zc.is.4)

any(conv_all.zc.is.4) 

## Can we reintro corr anywhere??

#no conv
# if (F) {
#   choice.lm.is <- llply(as.integer(names(conv_all.zc.is)),
#                           function(x) {
#                             if (!conv_all.zc.is[as.character(x)]) {
#                               print(x)
#                               this_data <- subset(pd_choicedf.is,trial_time_choice==x)
#                               this_choice_lm <- lmer(pupilCblz2~
#                                                        congruent.f + aSNR.f + congruent.f:aSNR.f + 
#                                                        isH.f + blz + posX + posY + zptlen + 
#                                                        (1 + congruent.f1 + aSNR.f1 + congruent_aSNR.f1 + 
#                                                           isH.f1 + blz + posX + posY + zptlen|dataID),
#                                                      data=this_data,
#                                                      control=lmerControl(optimizer="bobyqa",
#                                                                          optCtrl=list(maxfun=2e5)))
#                             }
#                           },.progress="text")
#   names(choice.lm.is) <- names(conv_all.zc.is)
#   choice.lm.is <- choice.lm.is[!conv_all.zc.is]
#   save(choice.lm.is,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_is_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_is_2021-06-10.rda'))
# }
# 
# conv_all.is <- get_conv(choice.lm.is)

#no convergence
# if (F) {
#   choice.lm.is.2 <- llply(as.integer(names(conv_all.zc.is.2)),
#                           function(x) {
#                             if (!conv_all.zc.is.2[as.character(x)]) {
#                               print(x)
#                               this_data <- subset(pd_choicedf.is,trial_time_choice==x)
#                               this_choice_lm <- lmer(pupilCblz2~
#                                                        congruent.f + aSNR.f + congruent.f:aSNR.f + 
#                                                        isH.f + blz + posX + posY + zptlen + 
#                                                        (1 + aSNR.f1 + congruent_aSNR.f1 + 
#                                                           isH.f1 + blz + posX + posY + zptlen|dataID),
#                                                      data=this_data,
#                                                      control=lmerControl(optimizer="bobyqa",
#                                                                          optCtrl=list(maxfun=2e5)))
#                             }
#                           },.progress="text")
#   names(choice.lm.is.2) <- names(conv_all.zc.is.2)
#   choice.lm.is.2 <- choice.lm.is.2[!conv_all.zc.is.2]
#   save(choice.lm.is.2,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_is2_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_is2_2021-06-10.rda'))
# }
# 
# conv_all.is.2 <- get_conv(choice.lm.is.2)

#no conv
# if (T) {
#   choice.lm.is.3 <- llply(as.integer(names(conv_all.zc.is.3)),
#                           function(x) {
#                             if (!conv_all.zc.is.3[as.character(x)]) {
#                               print(x)
#                               this_data <- subset(pd_choicedf.is,trial_time_choice==x)
#                               this_choice_lm <- lmer(pupilCblz2~
#                                                        congruent.f + aSNR.f + congruent.f:aSNR.f + 
#                                                        isH.f + blz + posX + posY + zptlen + 
#                                                        (1 + congruent_aSNR.f1 + 
#                                                           isH.f1 + blz + posX + posY + zptlen|dataID),
#                                                      data=this_data,
#                                                      control=lmerControl(optimizer="bobyqa",
#                                                                          optCtrl=list(maxfun=2e5)))
#                             }
#                           },.progress="text")
#   names(choice.lm.is.3) <- names(conv_all.zc.is.3)
#   choice.lm.is.3 <- choice.lm.is.3[!conv_all.zc.is.3]
#   save(choice.lm.is.3,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_is3_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_is3_2021-06-10.rda'))
# }
# 
# conv_all.is.3 <- get_conv(choice.lm.is.3)

#no conv
# if (F) {
#   choice.lm.is.4 <- llply(as.integer(names(conv_all.zc.is.4)),
#                           function(x) {
#                             if (!conv_all.zc.is.4[as.character(x)]) {
#                               print(x)
#                               this_data <- subset(pd_choicedf.is,trial_time_choice==x)
#                               this_choice_lm <- lmer(pupilCblz2~
#                                                        congruent.f + aSNR.f + congruent.f:aSNR.f + 
#                                                        isH.f + blz + posX + posY + zptlen + 
#                                                        (1 + congruent_aSNR.f1 + 
#                                                           blz + posX + posY + zptlen|dataID),
#                                                      data=this_data,
#                                                      control=lmerControl(optimizer="bobyqa",
#                                                                          optCtrl=list(maxfun=2e5)))
#                             }
#                           },.progress="text")
#   names(choice.lm.is.4) <- names(conv_all.zc.is.4)
#   choice.lm.is.4 <- choice.lm.is.4[!conv_all.zc.is.4]
#   save(choice.lm.is.4,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_is4_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_is4_2021-06-03.rda'))
# }
# 
# conv_all.is.4 <- get_conv(choice.lm.is.4)

#no conv
# if (F) {
#   choice.lm.is.5 <- llply(as.integer(names(conv_all.zc.is.5)),
#                           function(x) {
#                             if (!conv_all.zc.is.5[as.character(x)]) {
#                               print(x)
#                               this_data <- subset(pd_choicedf.is,trial_time_choice==x)
#                               this_choice_lm <- lmer(pupilCblz2~
#                                                        congruent.f + aSNR.f + congruent.f:aSNR.f +
#                                                        isH.f + blz + posX + posY + zptlen +
#                                                        (1 +
#                                                           blz + posX + posY + zptlen|dataID),
#                                                      data=this_data,
#                                                      control=lmerControl(optimizer="bobyqa",
#                                                                          optCtrl=list(maxfun=2e5)))
#                             }
#                           },.progress="text")
#   names(choice.lm.is.5) <- names(conv_all.zc.is.5)
#   choice.lm.is.5 <- choice.lm.is.5[!conv_all.zc.is.5]
#   save(choice.lm.is.5,file=paste0(DATA_OUT_PATH,'ptO_choice_lm_is5_', Sys.Date(), '.rda'))
# } else {
#   load(file=paste0(DATA_OUT_PATH,'ptO_choice_lm_is5_2021-05-30.rda'))
# }
# 
# conv_all.is.5 <- get_conv(choice.lm.is.5)


#now the painful merge
choice.lm.is.final <- choice.lm.zc.is

if (all(names(choice.lm.zc.is)[conv_all.zc.is]==names(choice.lm.zc.is.2)) & 
    all(names(choice.lm.zc.is.2)[conv_all.zc.is.2]==names(choice.lm.zc.is.3)) &
    all(names(choice.lm.zc.is.3)[conv_all.zc.is.3]==names(choice.lm.zc.is.4))) {
  
  choice.lm.is.final[conv_all.zc.is] <- choice.lm.zc.is.2
  choice.lm.is.final[names(choice.lm.zc.is.2)[conv_all.zc.is.2]] <- choice.lm.zc.is.3
  choice.lm.is.final[names(choice.lm.zc.is.3)[conv_all.zc.is.3]] <- choice.lm.zc.is.4
} else {
  error('Data mismatch between full and reduced models! Cannot combine.')
}

#final check (not a conclusive check but it's something)
if (!all(sapply(choice.lm.is.final,function(x) nrow(model.matrix(x))) == 
         sapply(choice.lm.zc.is,function(x) nrow(model.matrix(x)))) | 
    any(get_conv(choice.lm.is.final))) {
  stop('Data mismatch between full and reduced models!')
}

#extract coefficients/contrasts
choice.lm.is.final.sum <- sum_to_df(choice.lm.is.final,id='time',.pcorr=F,.confint=T,id2num=T)

choice.emm.is.final <- llply(choice.lm.is.final,
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

choice.emm.is.final.sum <- emm_to_df(choice.emm.is.final,id="time",.confint=T,id2num=T)

#OUTPUT
write.csv(choice.lm.is.final.sum,paste0("choice_ptO_lm_is_final_",Sys.Date(),".csv"),row.names=F)
write.csv(choice.emm.is.final.sum,paste0("choice_ptO_emm_is_final_",Sys.Date(),".csv"),row.names=F)