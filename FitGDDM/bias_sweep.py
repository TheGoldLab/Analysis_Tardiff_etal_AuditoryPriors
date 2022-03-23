# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 12:24:21 2021

Used to solve individual subject models at different values of the bias parameters in order
to construct an optimality landscape

@author: ntard
"""

import helpers
import fitinfo

import gddmwrapper as gdw
import pandas as pd
import numpy as np
import os.path
import copy
import argparse

from ddm import set_N_cpus
from datetime import date

def main():
    FIT_DIR = './fits'
    SWEEP_DIR = os.path.join(FIT_DIR,'sweeps')
    
    #parse command line arguments
    parser = argparse.ArgumentParser(description='Parameter sweep for AuditoryPriors optimality landscapes.')
    parser.add_argument('--subject',action='store',dest='subnum',type=int,required=False)
    parser.add_argument('--session',action='store',dest='sess',required=True)
    parser.add_argument('--parallel',action='store',dest='par',type=int)

    args = parser.parse_args()           
    
    if args.par:
        set_N_cpus(args.par)
    
    if args.sess=='priorOnly':
        model_info = fitinfo.precue_models['m3lbn']
    elif args.sess=='pretoneOnly':
        model_info = fitinfo.pretone_models['m14']
    else:
        raise NotImplementedError('Only priorOnly and pretoneOnly blocks supported!')
    
    models,samples=gdw.load_subjects(model_info.data,FIT_DIR,model_info.model_name,
                                     data_loader=model_info.data_loader)
    
    #sanity checks
    if samples.keys()!=models.keys():
        raise ValueError('Samples and modesl do not match!')
    
    #mnum_check = np.array([isinstance(m,list) for m in models.values()])
    mnum_check = np.array([len(m) if isinstance(m,list) else -1 
                           for m in models.values()])
    mnum_bad = mnum_check!=-1
    if mnum_bad.any():
        raise ValueError('Multiple models or missing models for some subjects!')
    
    #this should maybe go in loading code but complicated by potential to return
    #multiple models so doing here for now...
    msubj_check = np.array([v.subject==k for k,v in models.items()])
    if not msubj_check.all():
        raise ValueError('model.subject and model subject key mismatch present!')
    
    
    #set up variables
    #get some info from the data
    data = pd.read_csv(model_info.data)
    subj = np.sort(data.subject.unique())
    SNR = np.sort(data.SNR.unique())
    
    priors = dict()
    if args.sess=='priorOnly':
        priors['P'] = [1/6,.5,5/6]
        priors['Label'] = ['Low','No','High']
        priors['Cond'] = np.sort(data.prior.unique())
        v_vars = ['v_' + x for x in priors['Label']]
        z_vars = ['z_' + x for x in priors['Label']]
        biasvar = 'prior'
        cond_filter = None
    else:
        priors['P'] = [.5,.5]
        priors['Label'] = ['LL','HH']
        priors['Cond'] = [11,22]
        v_vars = 'v_Bias'
        z_vars = 'z_Bias'
        biasvar = 'bias'
        cond_filter = {biasvar:priors['Cond']}

    
    n=55 #set granularity of sweep
    zs = np.linspace(0, 1, n) 
    mes = np.linspace(-4, 4, n)

    
    #convert zs to logit scale for pretoneOnly, trimming ends. A bit janky frankly
    logit = lambda x: np.log(x/(1-x))
    if args.sess=='pretoneOnly':
        zs = logit(zs[1:-1])
    

    #construct proportions
    prop = []
    for prior,pl in zip(priors['P'],priors['Cond']):
            
        proportions = np.ones(SNR.shape)
        Lleft = SNR<0
        Lright = SNR>=0
        proportions[Lleft] = (1-prior)/np.sum(Lleft)
        proportions[Lright] = prior/np.sum(Lright)
        assert(proportions.sum()==1.0)
        prop.append(pd.DataFrame({biasvar:pl, 'SNR':SNR,
                                  'proportions':proportions}))
    prop = pd.concat(prop)
    
    if args.subnum is not None:
        runlist = subj[[args.subnum]]
    else: #default to running all
        runlist = subj
        
    
    #loop!
    for s in runlist:
        print(f'Running prediction sweep for subject {s}')
        this_mod = copy.deepcopy(models[s])

        paramvals = np.asarray(this_mod.get_model_parameters())
        paramnames = np.asarray(this_mod.get_model_parameter_names())
        
        preds = []
        #pcs = {p:np.zeros([n,n]) for p in priors['Cond']}
        pcs = []
        for d1,zz in enumerate(zs):
            print(f'z={zz}')
            for d2,mm in enumerate(mes):
                #print(f'z={zz},v={mm}')
                this_pval = paramvals.copy()
                this_pval[np.isin(paramnames,v_vars)] = mm
                this_pval[np.isin(paramnames,z_vars)] = zz
                this_mod.set_model_parameters(this_pval)
                pred=helpers.format_predicted(
                    gdw.get_predicted(this_mod,samples[s],cond_filter=cond_filter))

                pred[biasvar]=pred[biasvar].astype('int')
                
                pred['v'] = mm
                pred['z'] = zz
                
                preds.append(pred)
                
                #compute pct correct point on landscape
                pred = pred.merge(prop,on=[biasvar,'SNR'])
                pcts = pred.groupby([biasvar]).apply(
                    lambda subf: np.dot(subf['mean_corr'],subf['proportions']))
                pcts = pcts.reset_index()
                pcts.columns = (biasvar,'mean_corr')
                pcts['subject'] = s
                pcts['v'] = mm
                pcts['z'] = zz
                pcs.append(pcts)

                    
        #write results
        out_suffix = f'_{args.sess}_{s}_z{zs[0]}_{zs[-1]}_v{mes[0]}_{mes[-1]}_n{n}_{date.today()}.csv'

        preds = pd.concat(preds)
        preds.to_csv(os.path.join(SWEEP_DIR,'predsweep'+out_suffix),
                     index=False)
        
        pcs = pd.concat(pcs)
        pcs.to_csv(os.path.join(SWEEP_DIR,'pcs'+out_suffix),
             index=False)
