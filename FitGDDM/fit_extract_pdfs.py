# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 11:04:11 2020

@author: ntard
"""

import helpers
import fitinfo

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
    
    model_info = fitinfo.pretone_models['m0intlb'] #fitinfo.precue_models['m3lbn'] #fitinfo.all_models['m16pt5']  #fitinfo.all_models['m16pt5'] #fitinfo.precue_models['m3lb'] #fitinfo.all_models['m16'] #fitinfo.pretone_models['m11lb-14']
    output = True
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
    pdfs = [gdw.get_pdfs(mod,samp,method=model_info.method) 
            for mod,samp in zip(models.values(),samples.values())]
    
    pdf_df = pd.concat(pdfs)
    
     
    if output:
        pdf_outfile = '{model_name}_pdfs_{date}.csv'.format(
            model_name=model_info.model_name,  
            date=date.today())
        pdf_df.to_csv(os.path.join(FIT_DIR,pdf_outfile),index=False)
 
    
if __name__=='__main__':
    main()