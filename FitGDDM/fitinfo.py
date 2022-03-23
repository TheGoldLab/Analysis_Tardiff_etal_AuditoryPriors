# -*- coding: utf-8 -*-
"""
Created on Wed May  6 16:02:25 2020

@author: ntard


Parameters used in fit_eval.py for producing model output.

"""
import helpers


class FitInfo:
    
    def __init__(self,sess,data,data_loader,model_name,rundate=None):
        self.sess = sess
        self.data = data
        self.data_loader = data_loader
        self.model_name = model_name
        self.rundate = rundate
        
        
class PriorFitInfo(FitInfo):
    
    def __init__(self,model_name,rundate=None):
        super().__init__(sess='priorOnly',
                         data='./data/priorOnly_28-Jan-2020.csv',
                         data_loader=helpers.load_data_priorOnly,
                         model_name=model_name,
                         rundate=rundate)
        
class PretoneFitInfo(FitInfo):
    
    def __init__(self,data_loader,model_name,rundate=None):
        super().__init__(sess='pretoneOnly',
                         data='./data/pretoneOnly_26-Apr-2020.csv',
                         data_loader=data_loader,
                         model_name=model_name,
                         rundate=rundate)
        
class PT5FitInfo(FitInfo):
    def __init__(self,model_name,rundate=None):
        super().__init__(sess='pretone5',
                         data='./data/pretone5_28-Sep-2020.csv',
                         data_loader=lambda d,**kwargs: \
                                        helpers.load_data_all(
                                            d,14,**kwargs),
                         model_name=model_name,
                         rundate=rundate)
class PTLHFitInfo(FitInfo):
    def __init__(self,model_name,
                 data='./data/pretone_pL_pretone_pH_28-Sep-2020.csv',
                 rundate=None):
        super().__init__(sess='pretone_pLH',
                         data=data,
                         data_loader=lambda d,**kwargs: \
                                        helpers.load_data_all(
                                            d,14,**kwargs),
                         model_name=model_name,
                         rundate=rundate)
        
precue_models = {
    'm0int': PriorFitInfo(model_name='m0int-1_priorOnly'),
    'm0intlb': PriorFitInfo(model_name='m0intlb-1_priorOnly'),
    'm3lbn' : PriorFitInfo(model_name='m3lb-1'), #fit w/ numerical solver (using this)
    'm3lba' : PriorFitInfo(model_name='m3lb-1_priorOnly'), #fit w/ analytical solver
    }
pretone_models = {
    'm0int': PretoneFitInfo(
        data_loader=helpers.load_data_pretoneOnly_prior,
        model_name='m0int-1_pretoneOnly'),
    'm0intlb': PretoneFitInfo(
        data_loader=helpers.load_data_pretoneOnly_prior,
        model_name='m0intlb-1_pretoneOnly'),
    'm14': PretoneFitInfo(
        data_loader = lambda d,**kwargs: \
                helpers.load_data_pretoneOnly_ind(
                    d,14,**kwargs),
        model_name='m14-1_pretoneOnly'),
    'm17': PretoneFitInfo(
        data_loader = lambda d,**kwargs: \
                helpers.load_data_pretoneOnly_ind(
                    d,14,**kwargs),
        model_name='m17-1_pretoneOnly'),
    }

all_models = {
    'm16pt5': PT5FitInfo(model_name='m16pt5-1_pretone5'),
    'm16ptlh': PTLHFitInfo(model_name='m16ptlh-1_pretone_pLH'),
    }