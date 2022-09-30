# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 10:11:44 2020

This file includes custom model components for fitting GDDMs,
based in part on code from: https://pyddm.readthedocs.io/en/latest/cookbook/index.html

@author: ntard
"""

import numpy as np
from ddm.models import Drift,InitialCondition,Overlay,ICPoint
from ddm import Solution


class DriftCohBias(Drift):
    name = "Drift depends linearly on coherence, with additive prior bias"
    required_parameters = ["v_SNR","v_High","v_No","v_Low"] # <-- Parameters we want to include in the model
    #required_parameters = ["v_SNR","v_High","v_Low"]
    required_conditions = ["SNR","prior","isH"] # <-- Task parameters ("conditions"). Should be the same name as in the sample.
    
    # We must always define the get_drift function, which is used to compute the instantaneous value of drift.
    def get_drift(self, conditions, **kwargs):
        v_switch = {
            -2 : self.v_Low,
            0 : self.v_No,
            2 : self.v_High
        }
        v_bias = v_switch[conditions['prior']]
        # We need to change the directoin of the biases to reflect accuracy coding
        # Note that we want v_Low to be negative, as in stimulus coding
        # ie when is high, should push toward the error bound.
        if not conditions['isH']:
            v_bias = -v_bias
        
        return self.v_SNR * conditions['SNR'] + v_bias
 
 
#FROM: https://pyddm.readthedocs.io/en/stable/cookbook/driftnoise.html#drift-uniform
RESOLUTION = 11    
class DriftUniformCohBias(Drift):
    """Drift with trial-to-trial variability.
    
    Note that this is a numerical approximation to trial-to-trial
    drift rate variability and is inefficient.  It also requires
    running the "prepare_sample_for_variable_drift" function above on
    the data.
    """
    name = "Uniformly-distributed drift, Drift depends linearly on coherence, with additive prior bias"
    resolution = RESOLUTION # Number of bins, should be an odd number
    required_parameters = ["v_SNR","v_High","v_No","v_Low","v_Width"] # Mean drift and the width of the uniform distribution
    required_conditions = ["SNR","prior","isH","driftnum"]
    def get_drift(self, conditions, **kwargs):
        v_switch = {
            -2 : self.v_Low,
            0 : self.v_No,
            2 : self.v_High
        }
        v_bias = v_switch[conditions['prior']]
        # We need to change the directoin of the biases to reflect accuracy coding
        # Note that we want v_Low to be negative, as in stimulus coding
        # ie when is high, should push toward the error bound.
        if not conditions['isH']:
            v_bias = -v_bias
        
        #pick drift from uniform
        stepsize = self.v_Width/(self.resolution-1)
        mindrift = self.v_SNR - self.v_Width/2
        v_snr = mindrift + stepsize*conditions['driftnum']
        
        return v_snr * conditions['SNR'] + v_bias
    
class DriftCohInt(Drift):
    name = "Drift depends linearly on coherence, with additive offset"
    required_parameters = ["v_SNR","v_No"] # <-- Parameters we want to include in the model
    required_conditions = ["SNR","isH"] # <-- Task parameters ("conditions"). Should be the same name as in the sample.
    
    # We must always define the get_drift function, which is used to compute the instantaneous value of drift.
    def get_drift(self, conditions, **kwargs):
        
        v_bias = self.v_No
        
        # We need to change the directoin of the biases to reflect accuracy coding
        # Note that we want v_Low to be negative, as in stimulus coding
        # ie when is high, should push toward the error bound.
        if not conditions['isH']:
            v_bias = -v_bias
        
        return self.v_SNR * conditions['SNR'] + v_bias
 
            
class DriftCohExpBias1pExpAdaptPT(Drift):
    name = "Drift depends linearly on coherence, exponentially-weighted pretone bias/adaptation, shared bias weight"
    required_parameters = ["v_SNR","v_No",
                           "v_Bias","tau_Bias",
                           "wa_High","wa_Low","tau_Adapt"] # <-- Parameters we want to include in the model
    #required_parameters = ["v_SNR","v_High","v_Low"]
    required_conditions = ["SNR","bias","isH"] # <-- Task parameters ("conditions"). Should be the same name as in the sample.

    
    # We must always define the get_drift function, which is used to compute the instantaneous value of drift.
    def get_drift(self, conditions, **kwargs):
        
            
        v_intercept = self.v_No
        v_Bias = self.v_Bias
        wa_High = self.wa_High
        wa_Low = self.wa_Low
            
        # We need to change the direction of the biases to reflect accuracy coding
        # Note that we want v_Low to be negative, as in stimulus coding
        # ie when is high, should push toward the error bound.
        if not conditions['isH']:
            v_intercept*=-1
            v_Bias*=-1
            wa_High*=-1
            wa_Low*=-1
            
        #find positions of low and high pretones
        #pt_isH = np.asarray([p=='2' for p in str(conditions['bias'])])
        #pt_high = np.where(pt_isH)
        #pt_low = np.where(~pt_isH)
        
        pt_low,pt_high = pt_loc(conditions['bias'])
       
        #exponential bias function
        #v_bias_exp = pt_exp(self.tau_Bias,v_Bias,pt_high) - \
        #    pt_exp(self.tau_Bias,v_Bias,pt_low) 
        v_bias_exp = pt_bias_shared(self.tau_Bias,v_Bias,pt_low,pt_high)
            
        #exponential adaptation function
        v_adapt_exp = pt_exp(self.tau_Adapt,wa_High,pt_high) - \
            pt_exp(self.tau_Adapt,wa_Low,pt_low) 
        
        return v_intercept + self.v_SNR * conditions['SNR'] + \
            v_bias_exp + v_adapt_exp * conditions['SNR']

class DriftCohExpBias1pExpAdaptPTPrior(Drift):
    name = "Drift depends linearly on coherence, exponentially-weighted pretone bias/adaptation, shared bias weight, and additional prior"
    required_parameters = ["v_SNR","v_No","v_High","v_Low",
                           "v_Bias","tau_Bias",
                           "wa_High","wa_Low","tau_Adapt"] # <-- Parameters we want to include in the model
    #required_parameters = ["v_SNR","v_High","v_Low"]
    required_conditions = ["SNR","bias","isH","prior"] # <-- Task parameters ("conditions"). Should be the same name as in the sample.

    
    # We must always define the get_drift function, which is used to compute the instantaneous value of drift.
    def get_drift(self, conditions, **kwargs):
        
            
        #v_intercept = self.v_No
        v_Bias = self.v_Bias
        wa_High = self.wa_High
        wa_Low = self.wa_Low
        
        v_switch = {
            -2 : self.v_Low,
            0 : self.v_No,
            2 : self.v_High
        }
        v_intercept = v_switch[conditions['prior']]
       
        
        # We need to change the direction of the biases to reflect accuracy coding
        # Note that we want v_Low to be negative, as in stimulus coding
        # ie when is high, should push toward the error bound.
        if not conditions['isH']:
            v_intercept*=-1
            v_Bias*=-1
            wa_High*=-1
            wa_Low*=-1
        
        pt_low,pt_high = pt_loc(conditions['bias'])
       
        #exponential bias function
        v_bias_exp = pt_bias_shared(self.tau_Bias,v_Bias,pt_low,pt_high)
            
        #exponential adaptation function
        v_adapt_exp = pt_exp(self.tau_Adapt,wa_High,pt_high) - \
            pt_exp(self.tau_Adapt,wa_Low,pt_low) 
        
        return v_intercept + self.v_SNR * conditions['SNR'] + \
            v_bias_exp + v_adapt_exp * conditions['SNR']

#helper function to calculate pretone bias            

#convenience function for calculating bias for shared gain and tau between
#high/low pts
def pt_bias_shared(tau,gain,pt_low,pt_high):
    return pt_exp(tau,gain,pt_high) - \
            pt_exp(tau,gain,pt_low) 

def pt_exp(tau,gain,pt_pos):
    #using np.multiply b/c easy for scalars (tau,gain) to get 
    #implicitly converted to Pyton basic types without realizing it
    #even if you have at some point cast to np.array
    #and this causes errors/bad behavior
    if len(pt_pos)==0:
        #this will actually return 0 w/o this check b/c np.sum of empty list
        #returns 0, but this seems safer...
        return np.float64(0.0)
    else:
        pt_weight = np.exp(np.multiply(-tau,pt_pos))
        return gain*np.sum(pt_weight)

def pt_loc(pts):
    #find positions of low and high pretones
    if pts != 0:
        pt_isH = np.asarray([p=='2' for p in str(pts)])
        pt_high = np.where(pt_isH)
        pt_low = np.where(~pt_isH)
    else:
        pt_high = np.asarray([])
        pt_low = np.asarray([])
    
    return pt_low,pt_high


#based on: https://pyddm.readthedocs.io/en/latest/cookbook/initialconditions.html#ic-biased
#NOTE: original fits for legacy reasons used the numerical solver
#change superclass to ICPoint for analytical
#no, pass solving method to priorGDDM.py!
class ICPointPrior(ICPoint):
    name = "A prior-biased starting point."
    required_parameters = ["z_No","z_High","z_Low"]
    required_conditions = ["prior","isH"]
    def get_IC(self, x, dx, conditions):
        z_switch = {
            -2 : self.z_Low,
            0 : self.z_No,
            2 : self.z_High
        }
        z_bias = z_switch[conditions['prior']]
        # We need to change the directoin of the biases to reflect accuracy coding
        # Note that we want z_Low to be negative, as in stimulus coding
        # ie when is high, should push toward the error bound.
        if not conditions['isH']:
            z_bias = 1-z_bias
        shift_i = int((len(x)-1)*z_bias)
        assert shift_i >= 0 and shift_i < len(x), "Invalid initial conditions"
        pdf = np.zeros(len(x))
        pdf[shift_i] = 1. # Initial condition at x=self.x0.
        return pdf

#NOTE: for legacy reasons this uses the numerical solver
#change to ICPoint for analytical    
#no, pass solving method to priorGDDM.py!
class ICPointInt(ICPoint):
    name = "A biased starting point."
    required_parameters = ["z_No"]
    required_conditions = ["isH"]
    def get_IC(self, x, dx, conditions):
        
        z_bias = self.z_No
        # We need to change the directoin of the biases to reflect accuracy coding
        # Note that we want z_Low to be negative, as in stimulus coding
        # ie when is high, should push toward the error bound.
        if not conditions['isH']:
            z_bias = 1-z_bias
        shift_i = int((len(x)-1)*z_bias)
        assert shift_i >= 0 and shift_i < len(x), "Invalid initial conditions"
        pdf = np.zeros(len(x))
        pdf[shift_i] = 1. # Initial condition at x=self.x0.
        return pdf
    

class ICPointPTBias(ICPoint):
    name = "A biased starting point as a function of pt bias"
    required_parameters = ["z_No","z_Bias","tau_Bias"]
    required_conditions = ["bias","isH"]
    def get_IC(self, x, dx, conditions):
        
        
        #get low/high pretones
        pt_low,pt_high = pt_loc(conditions['bias'])
       
        #exponential bias function
        z_bias_exp = pt_bias_shared(self.tau_Bias,self.z_Bias,pt_low,pt_high)
        
        z_bias_lin = self.z_No + z_bias_exp
        z_bias = 1/(1+np.exp(-z_bias_lin))
        
        # We need to change the direction of the biases to reflect accuracy coding
        if not conditions['isH']:
            z_bias = 1-z_bias
        shift_i = int((len(x)-1)*z_bias)
        assert shift_i >= 0 and shift_i < len(x), "Invalid initial conditions"
        pdf = np.zeros(len(x))
        pdf[shift_i] = 1. # Initial condition at x=self.x0.
        return pdf
    
class ICPointPTBiasPrior(ICPoint):
    name = "A biased starting point as a function of pt bias with additional prior"
    required_parameters = ["z_No","z_High","z_Low","z_Bias","tau_Bias"]
    required_conditions = ["bias","prior","isH"]
    def get_IC(self, x, dx, conditions):
        #get bias-influenced intercept
        z_switch = {
            -2 : self.z_Low,
            0 : self.z_No,
            2 : self.z_High
        }
        z_intercept = z_switch[conditions['prior']]
        
        #get low/high pretones
        pt_low,pt_high = pt_loc(conditions['bias'])
       
        #exponential bias function
        z_bias_exp = pt_bias_shared(self.tau_Bias,self.z_Bias,pt_low,pt_high)
        
        z_bias_lin = z_intercept + z_bias_exp
        z_bias = 1/(1+np.exp(-z_bias_lin))
        
        # We need to change the direction of the biases to reflect accuracy coding
        if not conditions['isH']:
            z_bias = 1-z_bias
        shift_i = int((len(x)-1)*z_bias)
        assert shift_i >= 0 and shift_i < len(x), "Invalid initial conditions"
        pdf = np.zeros(len(x))
        pdf[shift_i] = 1. # Initial condition at x=self.x0.
        return pdf

   
#implements a nondecision time w/ a linear (absolute) dependence on 
#pt bias (individual pretones)
#adapted from: https://pyddm.readthedocs.io/en/latest/cookbook/nondecision.html#nd-lr    
class OverlayNonDecisionPTBiasCentered(Overlay):
    name = "Nondecision time as a function of pt bias (centered)"
    required_parameters = ["nondectime_No", "nondectime_Bias","tau_Bias"]
    required_conditions = ["bias"] # continuous measure of bias
    def apply(self, solution):
        # Unpack solution object
        corr = solution.corr
        err = solution.err
        m = solution.model
        cond = solution.conditions
        undec = solution.undec
        
        #compute pretone bias
        pt_low,pt_high = pt_loc(cond['bias'])
        
        #absolute value of exponential bias function
        #we use a weight of 1 before the absolute value and then multiply by
        #a gain term after
        #we then subtract the bias expected for two high pretones
        #as a rough way of mean-centering the absoulte bias
        apt_bias_exp = np.abs(pt_bias_shared(self.tau_Bias,1,pt_low,pt_high)) - \
            pt_exp(self.tau_Bias,1,np.asarray([0,1]))
        
        # Compute non-decision time
        ndtime = self.nondectime_No + self.nondectime_Bias * apt_bias_exp
        shifts = int(ndtime/m.dt) # truncate
        # Shift the distribution
        newcorr = np.zeros(corr.shape, dtype=corr.dtype)
        newerr = np.zeros(err.shape, dtype=err.dtype)
        if shifts > 0:
            newcorr[shifts:] = corr[:-shifts]
            newerr[shifts:] = err[:-shifts]
        elif shifts < 0:
            newcorr[:shifts] = corr[-shifts:]
            newerr[:shifts] = err[-shifts:]
        else:
            newcorr = corr
            newerr = err
        return Solution(newcorr, newerr, m, cond, undec)
