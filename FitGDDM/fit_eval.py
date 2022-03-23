# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 11:04:11 2020

@author: ntard


This script is used for processing fits and outputing csv files containing fit parameters,
predicted performance, and fit stats for all subjects.

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

from ddm import set_N_cpus
from datetime import date


def main():
    set_N_cpus(6)
    
    model_info = fitinfo.all_models['m16pt5'] #set model params specified in fitinfo here.
    output = True
    undec=True
    output_undec=False
    forced=True
    test=False
    test_dir = False
    FIT_DIR = './fits_test' if test_dir else './fits'

    
    models,samples=gdw.load_subjects(model_info.data,FIT_DIR,model_info.model_name,
                                     data_loader=model_info.data_loader)

    #sanity checks
    if samples.keys()!=models.keys():
        raise ValueError('Samples and modesl do not match!')
    
    mnum_check = np.array([len(m) if isinstance(m,list) else -1 
                           for m in models.values()])
    mnum_bad = mnum_check!=-1
    if mnum_bad.any():
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
    
    
    #are we hitting the minval-maxval boundary for any params?
    params = [gdw.get_params(mod,diagnostics=True) for mod in models.values()]
    params = pd.concat(params)
    

    if output:
        param_outfile = '{model_name}_params_{date}.csv'.format(
            model_name=model_info.model_name,  
            date=date.today())
        params.to_csv(os.path.join(FIT_DIR,param_outfile),index=True)
    
    bound_issues = params[params.hit_boundary==True]
    if len(bound_issues)>0:
        print('BOUNDARY ISSUES DETECTED!')
        print(bound_issues)
        
    #visualize params
    params.hist(by='param')
    plt.show()
    
    #extract fit indicies and save
    fit_stats = pd.DataFrame.from_records(
        gdw.get_fit_stats(m) for m in models.values())
    fit_stats.sort_values('subject',inplace=True)
    
    if output:
        fitstats_outfile = '{model_name}_lle_{date}.csv'.format(
            model_name=model_info.model_name,  
            date=date.today())
        fit_stats.to_csv(os.path.join(FIT_DIR,fitstats_outfile),index=False)
    
    #extract solutions
    soldf = [helpers.format_predicted(gdw.get_predicted(mod,samp,
                                                        undec=undec,forced=forced)) 
            for mod,samp in zip(models.values(),samples.values())]
    
    soldf = pd.concat(soldf)
    
    
    if undec:
        soldf.mean_undec.hist()
        plt.show()
        
    
    if output:
        sol_outfile = '{model_name}_predicted_{date}.csv'.format(
            model_name=model_info.model_name,  
            date=date.today())
        if not output_undec:
            soldf.drop('mean_undec',axis=1,inplace=True)
        soldf.to_csv(os.path.join(FIT_DIR,sol_outfile),index=False)
    

  
    
if __name__=='__main__':
    main()