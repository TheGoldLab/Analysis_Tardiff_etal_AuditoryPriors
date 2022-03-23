#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:47:26 2020
    
@author: ntardiff


Command-line program for fitting GDDM models to AuditoryPriors data.

"""

import model_defs
import helpers

import argparse
import copy
import math
import gddmwrapper as gdw
from ddm.functions import display_model
from ddm import set_N_cpus


def main():
    #set constants
    DATA_FILE_DEFAULT = './data/priorOnly_28-Jan-2020.csv'
    FITS_DIR = './fits'
    FITS_TEST_DIR = './fits_test'
    PT_IND_MODELS = ['m14','m17']
    PT_MAX = 14
        
    #parse command line arguments
    parser = argparse.ArgumentParser(description='Fit PyDDM models.')
    parser.add_argument('--model_id',action='store',dest='model_id',required=True)
    parser.add_argument('--subject',action='store',dest='subnum',type=int,required=True)
    parser.add_argument('--iteration',action='store',dest='it',type=int)
    parser.add_argument('--session',action='store',dest='sess')
    parser.add_argument('--data',action='store',dest='data')
    parser.add_argument('--parallel',action='store',dest='par',type=int)
    parser.add_argument('--test',action='store_true')

    args = parser.parse_args()           
    
    if args.data:
        DATA_FILE = args.data
    else:
        DATA_FILE = DATA_FILE_DEFAULT
    
    if args.par:
        set_N_cpus(args.par)
        
    if args.test:
        OUT_DIR = FITS_TEST_DIR
        load_args = {'verbose': True}
    else:
        OUT_DIR = FITS_DIR
        load_args = dict()
    
    #get model
    model = copy.deepcopy(getattr(model_defs,args.model_id))

    #load data
    if args.sess=='priorOnly':
        data_loader = helpers.load_data_priorOnly
    elif args.sess=='pretoneOnly':
        #for pretones, loader depends on model
        if model.name in ['m0int','m0intlb']:
            data_loader= helpers.load_data_pretoneOnly_prior
        elif model.name in PT_IND_MODELS: 
            data_loader = lambda d,**kwargs: \
                helpers.load_data_pretoneOnly_ind(
                    d,PT_MAX,**kwargs)
        else:
            raise NotImplementedError('Undefined pretone model!')
    elif args.sess in ['all','pretone_pLH','pretone5']:
        data_loader = lambda d,**kwargs: \
            helpers.load_data_all(
                d,PT_MAX,**kwargs)
    else:
        data_loader = helpers.load_data_priorOnly #legacy behavior
    
    sample,subj = data_loader(DATA_FILE,subjnum=args.subnum,
                                              **load_args)
            

    #run model
    print('Fitting model %s for subject %s iteration %s session %s...' % 
          (model.name,subj,args.it,args.sess))
    print('Data file: %s\n' % DATA_FILE)
    if args.sess=='all' or model.name in PT_IND_MODELS: 
        print("PT_MAX: %d\n" % PT_MAX)
    
    if args.test:
        import time
        stime = time.time()
    
    #keep popsize ~40 
    popsize = math.ceil(40./len(model.get_model_parameter_names()))
    model_fit = gdw.run_model(sample,model,subj,args.it,args.sess,OUT_DIR,
                              fitparams={'popsize':popsize,'strategy':'rand1bin',
                                         'recombination':0.9})
    
    if args.test:
        rtime = time.time() - stime
        print('Fitting time: %.2f s (%.2f h)' % (rtime,rtime/60**2))
        
    display_model(model_fit)
    
   
if __name__=='__main__':
    main()

