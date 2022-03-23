# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 13:23:40 2020

helper functions for AuditoryPriors data

@author: ntard
"""
import numpy as np
import pandas as pd
import gddmwrapper as gdw


#base preprocessing for all data
def base_preproc(data):
    #there shouldn't be any other missing data in other columns,
    #can't do a blanket dropna b/c of pretone columns...
    data.dropna(subset=['rt','correct'],inplace=True)
    data.loc[:,'rt'] = data.rt/1000. #rt to seconds
    data.loc[:,'SNR'] = data.SNR.abs() #pyDDM does not use stimulus coding

    return data

#preproc for blocks w/ pretones or pts and priors (only used for reduced models cause doesn't recode pretones)
#we remove no-pretone trials from pretone blocks b/c they are confusing
def pretoneOnly_prior_preproc(data):
    data = base_preproc(data)
    data = data[~((data.pretoneLength==0) & 
              (data.optionName.isin(
                  ('pretoneOnly','pretone_pL','pretone_pH'))))].copy()
    
    return data


#pretone only for modeling individual pretones
#note: must specify the columns containing the pretones as  slice object
#this allows control over how many pretones to pass
#note that for pt_range, numbering is time-reversed 
#[i.e, (closest in time to test tone, furthest in time to test tone)]
def pretoneOnly_ind_preproc(data,pt_max=14,bias_col='pretoneSeq'):
    data = pretoneOnly_prior_preproc(data)
    

    #covert Hs and Ls to '2s' and '1s' in reverse order, ignoring terminal 'X'
    pt_cats = data[bias_col].apply(
        lambda x: ['2' if s[0]=='H' else '1' 
                   for s in x[-2::-1]][:pt_max])
    pt_cats = [int(''.join(x)) if x else 0 for x in pt_cats] #join into ints,no pretones=0
    
    
    #add single pretone column
    data['bias'] = pt_cats
    
    return data


#loader for precue
def load_data_priorOnly(data_file,**kwargs):
    
    return gdw.load_data(data_file,rt='RT',correct='success',
                                conds=['SNR','prior','isH'],
                                preproc=base_preproc,
                                **kwargs)

#loader for reduced pretone models
def load_data_pretoneOnly_prior(data_file,**kwargs):
    
    return gdw.load_data(data_file,rt='RT',correct='success',
                                conds=['SNR','prior','isH'],
                                preproc=pretoneOnly_prior_preproc,
                                **kwargs)


#loader for pretone w/ individual pretones
def load_data_pretoneOnly_ind(data_file,pt_max,**kwargs):
    
    preproc = lambda d: pretoneOnly_ind_preproc(d,pt_max)
    
    return gdw.load_data(data_file,rt='RT',correct='success',
                                conds=['SNR','bias','isH'],
                                preproc=preproc,
                                **kwargs)

#loader for pretone w/ individual pretones + prior
def load_data_all(data_file,pt_max,**kwargs):
    
    preproc = lambda d: pretoneOnly_ind_preproc(d,pt_max)
    
    return gdw.load_data(data_file,rt='RT',correct='success',
                                conds=['SNR','bias','prior','isH'],
                                preproc=preproc,
                                **kwargs)

def format_predicted(soldf):
     '''
     reformat code produced by get_predicted for use in AuditoryPriors
     ''' 
     
     #try to figure out whether this is a precue or pretone model
     bias_cols = np.array(['prior','bias'])
     bias_here = np.isin(bias_cols,soldf.columns)
     
     if 'isH' not in soldf.columns:
         #make no bias model (m0) equivilant to other models
         #no bias means no asymmetry, so just mirror
         soldf_low = soldf.copy()
         soldf_low['SNR']*=-1
         soldf['isH']=1
         soldf_low['isH']=0
         soldf = pd.concat([soldf,soldf_low])
         
         #now make copies for the prior (only do this if bias not present)
         #might want to put this elsewhere if run a no bias model on pretone...
         if not bias_here.any():
             soldf_ps = list()
             for p in [-2,0,2]:
                 this_p = soldf.copy()
                 this_p['prior'] = p
                 soldf_ps.append(this_p)
             soldf = pd.concat(soldf_ps)
     else:
         #reintroduce -SNR for low tones
         soldf.loc[soldf.isH==0,'SNR']*=-1
     
     #sort df
     soldf.sort_values([*bias_cols[bias_here].tolist(),'SNR'],
                       inplace=True,ignore_index=True)

     
     #compute chooseH from cross of mean_corr/error and isH (i.e.,
     #for isH==0, choose_H is mean_err)
     #in principle mean_err = 1-mean_corr, but slightly off due to
     #predicted undecided trials.
     soldf['mean_chooseH'] = soldf['mean_corr']
     soldf.loc[soldf.isH==0,'mean_chooseH'] = soldf.loc[soldf.isH==0,'mean_err']
         
     return soldf
