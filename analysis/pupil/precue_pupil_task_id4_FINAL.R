#clear memory
rm(list=ls())

#This script runs individual-subject regressions that are subsequently used in 
#pupil_precue5_FINAL.ipynb to ask how individual differences in pupil response relate to
#individual differences in bias (Fig 5c).


## LOADING data/libraries ##

#load libraries
library(lme4)
library(lmerTest)
library(car)
library(plyr)
library(dplyr)
library(ggplot2)
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
PSYCHO_PATH = '../data/'
pd_file= paste0(DATA_OUT_PATH,'pupil_data_ds50_forR_09-Feb-2021.csv')
baseline_file = paste0(DATA_OUT_PATH,'pupil_data_bl_forR_09-Feb-2021.csv')
bias_file=paste0(PSYCHO_PATH,'fits_lapse_priorOnly_noint_01-Mar-2021.csv')
id_file=paste0(DATA_OUT_PATH,'pupil_data_ds50_id2subj_forR_09-Feb-2021.csv')

#source my functions if any
source('../sum_to_df.R')
source('../emm_to_df.R')
source('../get_conv.R')

#set temporal parameters 
choice_period_ms.c <- c(-250, 1000) 


#load data
pd_df<-read.table(pd_file,sep=',', header=TRUE, stringsAsFactors=FALSE,na.strings = c('NaN'))
bl_df<-read.table(baseline_file,sep=',', header=TRUE, stringsAsFactors=FALSE,na.strings = c('NaN'))
bias_df <- read.table(bias_file,sep=',', header=TRUE, stringsAsFactors=FALSE,na.strings = c('NaN'))
id_df <- read.table(id_file,sep=',', header=TRUE, stringsAsFactors=FALSE,na.strings = c('NaN'))

bias_df <- left_join(bias_df,id_df,by='subject')
#remove subjects w/o data IDs (i.e not in pupil data)
bias_df <- subset(bias_df,!is.na(bias_df$dataID))
rm(id_df)

head(pd_df)
head(bl_df)
head(bias_df)

#cleanup
bl_df$GroupCount <- NULL
pd_df[,c('time','pupilC','posXC','posYC')] <- NULL

#set up variables/factors
my_simple2<-contr.treatment(2,base=2) - matrix(rep(1/2,2))
my_simple3 <- contr.treatment(3,base=3) - matrix(rep(1/3,6),ncol=2)

pd_df$congruent.f <- factor(pd_df$congruent,levels=c(1,0,-1),
                            labels=c("congruent","incongruent","no prior"))
contrasts(pd_df$congruent.f) <- contr.sum(3)
contrasts(pd_df$congruent.f)

pd_df$congruent.fs <- factor(pd_df$congruent,levels=c(0,-1,1),
                             labels=c("incongruent","no prior","congruent"))
contrasts(pd_df$congruent.fs) <- my_simple3
contrasts(pd_df$congruent.fs)

pd_df$isH.fs <- factor(pd_df$isH,levels=c(1,0),
                       labels=c("high","low"))
contrasts(pd_df$isH.fs) <- my_simple2
contrasts(pd_df$isH.fs)
pd_df$isH.f <- pd_df$isH.fs
contrasts(pd_df$isH.f) <- contr.sum(2)
contrasts(pd_df$isH.f)

pd_df$success.f <- factor(pd_df$success,levels=c(1,0),labels=c("correct","incorrect"))
contrasts(pd_df$success.f) <- contr.sum(2)
contrasts(pd_df$success.f)

pd_df$aSNR <- abs(pd_df$SNR)


### RESPONSE REGRESSION ###

## CORRECT TRIALS ##
pd_df.c <- subset(pd_df,success==1 &
                    trial_time_choice >= choice_period_ms.c[1] & 
                    trial_time_choice <= choice_period_ms.c[2])
#gut check
min(pd_df.c$trial_time_choice)
max(pd_df.c$trial_time_choice)
min(pd_df.c$trial_time_stimOn)
unique(pd_df.c$success)

#get rid of missing before averaging/scaling (shouldn't be any but just in case)
pd_df.c <- pd_df.c[complete.cases(pd_df.c),]

#merge in baseline data
pd_df.c <- left_join(pd_df.c,bl_df[,c('dataID','trialN','pupilBL')],by=c('dataID','trialN'))

#set up variables: centering/scaling
pd_df.c$zaSNR <- scale(pd_df.c$aSNR)
pd_df.c$blz <- scale(pd_df.c$pupilBL,center=T,scale=F) #not scaling since already scaled
pd_df.c$posX <- scale(pd_df.c$posXCbl)
pd_df.c$posY <- scale(pd_df.c$posYCbl)

#set up variables for zero-correlation lme
pd_df.c[,c('congruent.f1','congruent.f2')] <- 
  model.matrix(~1+pd_df.c$congruent.f,pd_df.c)[,2:3]
pd_df.c$isH.f1 <- 
  model.matrix(~1+pd_df.c$isH.f,pd_df.c)[,2]

#gut check
unique(pd_df.c[,c('congruent','congruent.f','congruent.f1','congruent.f2')])
unique(pd_df.c[,c('isH.f','isH.f1')])


#looping variables
choice_times.c <- sort(unique(pd_df.c$trial_time_choice))
IDs <- unique(pd_df.c$dataID)

#subset bias by pupil subjects and create bias variables
bias_df <- bias_df[bias_df$dataID %in% IDs,]
bias_df$bias <- bias_df$bias_high - bias_df$bias_low
bias_df$zbias <- scale(bias_df$bias)


### RUN: single-subject LM at every time point implementing same regression as in mixed-effects
# choice model ###

#Run individual subject lms
if (F) {
  choice.lm_subj <- llply(choice_times.c,
                        function(t) {
                          this_time.lm<- llply(IDs,
                             function(x) {
                               this_data <- subset(pd_df.c,dataID==x & trial_time_choice==t)
                               #print(head(this_data))
                               tryCatch( 
                                 this_lm <- lm(pupilCblz~congruent.f + isH.f + zaSNR + blz + posX + posY,
                                               data=this_data),
                                 error=function(c) print(sprintf("%s time: %d, dataID: %d",c,t,x)),
                                 warning=function(c) print(sprintf("%s, time: %d, dataID: %d",c,t,x))
                               )
                             })
                          names(this_time.lm) <- IDs
                          #drop empty
                          valid_lm <- sapply(this_time.lm,function(x) is(x,"lm"))
                          return(this_time.lm[valid_lm])
                        },.progress="text")
  names(choice.lm_subj) <- choice_times.c
  save(choice.lm_subj,file=paste0(DATA_OUT_PATH,'choice_lm_subj_', Sys.Date(), '.rda'),compress = T)
} else {
  load(paste0(DATA_OUT_PATH,'choice_lm_subj_2021-03-31.rda'))
}

#extract contrasts
if (F) {
  choice.lm_subj.emm <- llply(choice.lm_subj,function(t) {
    llply(t,
          function(x) {
            emm <- emmeans(x,specs=~congruent.f)
            emm.contr <- contrast(emm,
                                  list("incongruent - congruent"=c(-1,1,0),
                                       "no prior - congruent"=c(-1,0,1),
                                       "no prior - incongruent"=c(0,-1,1)))  
            #"incongruent_no prior - congruent"=c(-1,1/2,1/2)))
          })  
  },.progress="text")
  save(choice.lm_subj.emm,file=paste0(DATA_OUT_PATH,'choice_lm_subj_emm_', Sys.Date(), '.rda'),
       compress = T)
} else {
  load(paste0(DATA_OUT_PATH,'choice_lm_subj_emm_2021-03-31.rda'))
}

choice.lm_subj.emm.sum <- ldply(choice.lm_subj.emm,function(x) emm_to_df(x,id="dataID",
                                                                         id2num = T,.confint=F))
choice.lm_subj.emm.sum <- rename(choice.lm_subj.emm.sum,time=.id)
choice.lm_subj.emm.sum$time <- as.numeric(choice.lm_subj.emm.sum$time)

#package contrast/bias data for output
conbiasdf <- left_join(choice.lm_subj.emm.sum[,c("time","dataID","param","B")],
                       bias_df[,c("dataID","zbias")])


#ALTERNATIVE estimates of conditions rather than condition contrasts (using simple lm) [DEPRECATED]

if (F) {
  choice.lm_subj.emm.1 <- llply(choice.lm_subj,function(t) {
    llply(t,
          function(x) {
            emm <- emmeans(x,specs=~congruent.f)
          })  
  },.progress="text")
  save(choice.lm_subj.emm.1,file=paste0(DATA_OUT_PATH,'choice_lm_subj_emm1_', Sys.Date(), '.rda'),
       compress = T)
} else {
  load(paste0(DATA_OUT_PATH,'choice_lm_subj_emm1_2021-04-01.rda'))
}

choice.lm_subj.emm.1.sum <- ldply(choice.lm_subj.emm.1,function(x) emm_to_df(x,id="dataID",
                                                                         id2num = T,.confint=F))
choice.lm_subj.emm.1.sum <- rename(choice.lm_subj.emm.1.sum,time=.id,B=param,param=congruent.f) #got names a bit wrong
choice.lm_subj.emm.1.sum$time <- as.numeric(choice.lm_subj.emm.1.sum$time)


#package for output
conbiasdf.1 <- left_join(choice.lm_subj.emm.1.sum[,c("time","dataID","param","B")],
                       bias_df[,c("dataID","zbias")])


#Export individual subject running regressions for looking at spearman corr w/ behavioral bias
#stitch together individual congruent cond estimates and congruent contrast estimates and save
conbiasdf_full <- rbind(conbiasdf,conbiasdf.1) 
conbiasdf_full <- arrange(conbiasdf_full,time,dataID,param)
write.csv(conbiasdf_full,file=paste0(DATA_OUT_PATH,'conbias_full_',Sys.Date(),'.csv'),row.names=F)
