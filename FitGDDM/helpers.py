# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 13:23:40 2020

helper functions for AuditoryPriors data

@author: ntard
"""
import numpy as np
import pandas as pd
import gddmwrapper as gdw
import copy
from ddm import Sample
from ddm.fitresult import FitResult
from ddm.models import LossLikelihood


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

#load data for precue w/ uniform drift variability
def load_data_priorOnly_dv(data_file,**kwargs):
    
    sample,subj =  gdw.load_data(data_file,rt='RT',correct='success',
                                conds=['SNR','prior','isH'],
                                preproc=base_preproc,
                                **kwargs)            
    if isinstance(sample,dict):
        sample = {k:prepare_sample_for_variable_drift(s) for k,s in sample.items()}
    else:
        sample = prepare_sample_for_variable_drift(sample)
    return sample,subj

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
     
#this is for preparing for fitting uniform drift rate variability
#FROM: https://pyddm.readthedocs.io/en/stable/cookbook/driftnoise.html#drift-uniform
def prepare_sample_for_variable_drift(sample, resolution=11):
    new_samples = []
    for i in range(0, resolution):
        corr = sample.corr.copy()
        err = sample.err.copy()
        undecided = sample.undecided
        conditions = copy.deepcopy(sample.conditions)
        conditions['driftnum'] = (np.asarray([i]*len(corr)),
                                  np.asarray([i]*len(err)),
                                  np.asarray([i]*undecided))
        new_samples.append(Sample(corr, err, undecided, **conditions))
    new_sample = new_samples.pop()
    for s in new_samples:
        new_sample += s
    return new_sample
     
#this function corrects the driftvar log-likelihood, which is incorrect in 
#model output due to the sample replication necessary to run this model
def driftvar_ll(model,sample):
    print('Correcting driftvar ll: This may take a minute...')
    #solve model
    sol = gdw.solve_some_conditions(model,sample)
    
    #get all driftnums (the discretization of drift variability)
    driftnums = frozenset(('driftnum',x) for x in sample.condition_values('driftnum'))
    #dict to translate from full sample conditions to reduced conditions w/o driftnum
    full_to_reduc = {c:c-driftnums for c in sol.keys()}
    
    fullconds = np.array(list(sol.keys())) #list of all conditions for filtering
    #dict to translate from reduced conditiions to full conditions
    reduc_to_full = {u:fullconds[[u.issubset(f) for f in fullconds]] 
              for u in frozenset(full_to_reduc.values())}
    
    #integrate out driftnum (aka average over driftnum pdfs to achieve reduced pdfs)
    pdf_corr_ave = {k:np.mean([sol[s].pdf_corr() for s in v],axis=0) 
                for k,v in reduc_to_full.items()}
    pdf_err_ave = {k:np.mean([sol[s].pdf_err() for s in v],axis=0) 
                    for k,v in reduc_to_full.items()}
    
    #let's do a bunch of sanity checks
    pdf_check = np.array([np.trapz(c,dx=model.dt)+np.trapz(e,dx=model.dt) 
                 for c,e in zip(pdf_corr_ave.values(),pdf_err_ave.values())])
    pdf_undec_ave = {k:np.mean([sol[s].prob_undecided() for s in v],axis=0) 
                    for k,v in reduc_to_full.items()} #add in undec prob so pdfs sum to ~1
    pdf_check += list(pdf_undec_ave.values())

    assert np.all((1-pdf_check) < 0.001), "pdfs do not sum to 1!" 
    #print("max pdf deviation from 1: %s" % np.max(np.abs(1-pdf_check)))
    assert np.all([s.undec==None for s in sol.values()]), "Handling of undecided trials not implemented!"
    
    corr_ave_check = {k:np.mean([sol[s].prob_correct() for s in v])
                    for k,v in reduc_to_full.items()}
    corr_ave_check2 = {k:np.trapz(v,dx=model.dt) for k,v in pdf_corr_ave.items()}
    assert np.all([(c1-c2)<0.001] for c1,c2 in zip(corr_ave_check.values(),corr_ave_check2.values())), "Mean prob corrects do not match"
    
    err_ave_check = {k:np.mean([sol[s].prob_error() for s in v])
                    for k,v in reduc_to_full.items()}
    err_ave_check2 = {k:np.trapz(v,dx=model.dt) for k,v in pdf_err_ave.items()}
    assert np.all([(c1-c2)<0.001] for c1,c2 in zip(err_ave_check.values(),err_ave_check2.values())), "Mean prob errors do not match"
       
    #get losslikelihood for only one replicate of sample (i.e. original sample)
    samp0 = sample.subset(driftnum=0)
    loss = LossLikelihood(samp0,dt=model.dt,T_dur=model.T_dur)
    
    #now compute log likelihood
    loglikelihood = 0
    for k in samp0.condition_combinations():
        #convert dict returned by condition_combinations to frozenset for use as dict key
        k = frozenset(k.items())
        
        #pdfs are indexed by empirical histograms
        #correct pdf is chosen by indexing averaged pdfs for each full condiiton using the
        #full_to_reduc dict
        loglikelihood += np.sum(np.log(pdf_corr_ave[full_to_reduc[k]][loss.hist_indexes[k][0]]))
        loglikelihood += np.sum(np.log(pdf_err_ave[full_to_reduc[k]][loss.hist_indexes[k][1]]))
    
    #create new fitresult object with corrected ll/sample size
    res = FitResult(method=model.get_fit_result().method,
                fitting_method=model.get_fit_result().fitting_method, 
                loss=model.get_fit_result().loss, value=-loglikelihood,
                nparams=model.get_fit_result().properties['nparams'], 
                samplesize=len(samp0),
                mess=model.get_fit_result().properties['mess'])
        
    
    return res #-loglikelihood
