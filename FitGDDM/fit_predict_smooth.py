# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 11:04:11 2020

@author: ntard
"""

import helpers
import fitinfo

#import ddm.plot
import gddmwrapper as gdw
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os.path
import warnings

#from ddm.functions import solve_partial_conditions,solve_all_conditions,display_model
from ddm import set_N_cpus
from datetime import date

#import paranoid as pns
#pns.settings.Settings.set(enabled=False)

def main():
    set_N_cpus(6)
    
    model_info = fitinfo.pretone_models['m14'] #fitinfo.precue_models['m3lbn'] #fitinfo.all_models['m16pt5']  #fitinfo.all_models['m16pt5'] #fitinfo.precue_models['m3lb'] #fitinfo.all_models['m16'] #fitinfo.pretone_models['m11lb-14']
    output = True
    undec=True
    output_undec=False
    forced=True
    test=False
    test_dir = False
    FIT_DIR = './fits_test' if test_dir else './fits'

    
    models,samples=gdw.load_subjects(model_info.data,FIT_DIR,model_info.model_name,
                                     data_loader=model_info.data_loader)
    #print(models['AdL'].has_analytical_solution())
    #sanity checks
    if samples.keys()!=models.keys():
        raise ValueError('Samples and modesl do not match!')
    
    #mnum_check = np.array([isinstance(m,list) for m in models.values()])
    mnum_check = np.array([len(m) if isinstance(m,list) else -1 
                           for m in models.values()])
    mnum_bad = mnum_check!=-1
    if mnum_bad.any():
        #subj = np.array(list(models.keys()))
        #print(subj[mnum_check > 1])
        if test:
            warnings.warn('Testing mode: removing missing/multiple models!')
            to_del=np.asarray(list(models.keys()))[mnum_bad]
            for k in to_del:
                models.pop(k,None)
                samples.pop(k,None)
        else:
            raise ValueError('Multiple models or missing models for some subjects!')
    
    #this should maybe go in loading code but complicated by potential to return
    #multiple models so doing here for now...
    msubj_check = np.array([v.subject==k for k,v in models.items()])
    if not msubj_check.all():
        raise ValueError('model.subject and model subject key mismatch present!')
    
    
   
    #extract solutions
    #sol = gdw.get_predicted(models['AdL'],samples['AdL'])
    #I hope after above checks this is safe...
    cond_filter = None# {'SNR': [0.05,0.5]}
    cond_replace = None #{'SNR':list(np.hstack((0.025,np.linspace(0.05,0.5,10))))}
    cond_augment = {'SNR':list(np.hstack((0.025,np.linspace(0.05,0.5,10))))}
    soldf = [helpers.format_predicted(gdw.get_predicted(mod,samp,
                                                        undec=undec,forced=forced,
                                                        cond_filter=cond_filter,
                                                        cond_replace=cond_replace,
                                                        cond_augment=cond_augment)) 
            for mod,samp in zip(models.values(),samples.values())]
    
    soldf = pd.concat(soldf)
    
    #fix annoying int32 issue with bias values
    if model_info.model_name in ['m12-1_pretoneOnly','m11lb-14_pretoneOnly']:
        with open('pt_int_dict.pickle',"rb") as f:
            import pickle
            pt_int_dict=pickle.load(f)
        soldf['bias'].replace(to_replace=pt_int_dict,inplace=True)
    
    if undec:
        soldf.mean_undec.hist()
        plt.show()
        
    if model_info.model_name=='m14-1_pretone_pL':
        soldf['prior'] = -2
    elif model_info.model_name=='m14-1_pretone_pH':
        soldf['prior'] = 2
    
    if output:
        sol_outfile = '{model_name}_predicted_smooth_{date}.csv'.format(
            model_name=model_info.model_name,  
            date=date.today())
        if not output_undec:
            soldf.drop('mean_undec',axis=1,inplace=True)
        soldf.to_csv(os.path.join(FIT_DIR,sol_outfile),index=False)
    

  
    
if __name__=='__main__':
    main()