# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 15:37:26 2020

Defines models used to fit the AuditoryPriors dataset

@author: ntard
"""

from gddmwrapper.models import DriftCohInt,DriftCohBias,DriftCohExpBias1pExpAdaptPT,DriftCohExpBias1pExpAdaptPTPrior,DriftUniformCohBias
from gddmwrapper.models import ICPointInt,ICPointPrior,ICPointPTBias,ICPointPTBiasPrior
from gddmwrapper.models import OverlayNonDecisionPTBiasCentered
from ddm import Model,Fittable
from ddm.models import NoiseConstant, BoundConstant, BoundCollapsingLinear
from ddm.models import OverlayChain, OverlayNonDecision, OverlayUniformMixture

#common model parameter ranges
DX_ALT = .005 #####
DT_ALT2 = .005 ####
T_DUR = 3 ####
T_DUR_PT = 2.05 ####
NOISE = 1 ####

DRIFT_MIN = 0 ######
DRIFT_MAX2 = 18 #####
DRIFT_VAR_MIN = 0 #####
DRIFT_VAR_MAX = 5 #####

BOUND_MIN=.1 ######
BOUND_MAX_ALT = 2.5 ######
BOUND_T_MIN = 0 ######
BOUND_T_MAX = BOUND_MAX_ALT ######

ND_MIN = 0 #####
ND_MAX = .5 #####
ND_MAX_PT = .7 #####
U_MIX_MIN = 0 ######
U_MIX_MIN2 = 0.00005 ######
U_MIX_MAX = .1 ######

IC_DEFAULT = 0.5 ######
IC_MIN2 = .1 ######
IC_MAX2 = .9 ######

IC_LOG_MIN3 = -1.8 #####
IC_LOG_MAX3 = 1.8 #####
IC_LOG_DEFAULT = 0 #####

DRIFT_BIAS_MIN2 = -5 #####
DRIFT_BIAS_MAX2 = 5 ######

ND_BIAS_MIN = -ND_MAX_PT ####
ND_BIAS_MAX = ND_MAX_PT ####

DRIFT_WB_MIN2 = DRIFT_BIAS_MIN2 #####
DRIFT_WB_MAX2 = DRIFT_BIAS_MAX2 #####
DRIFT_WA_MIN2 = DRIFT_WB_MIN2*2 #####
DRIFT_WA_MAX2 = DRIFT_WB_MAX2*2 #####

DRIFT_TAU_MIN = 0 ####
DRIFT_TAU_MAX = 20 ####

#Model definitions

#basic DDM,drift varies with SNR


#basic DDMs w/ biases for drift and starting point that do not depend on priors
m0int = Model(name='m0int',
             drift=DriftCohInt(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                               v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2)),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundConstant(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT)),
             IC=ICPointInt(z_No=Fittable(minval=IC_MIN2, maxval=IC_MAX2,default=IC_DEFAULT)),
             overlay=OverlayChain(overlays=[OverlayNonDecision(
                 nondectime=Fittable(minval=ND_MIN, maxval=ND_MAX)),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR)

#m0int for pretoneOnly
m0intpt = Model(name='m0int',
             drift=DriftCohInt(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                               v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2)),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundConstant(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT)),
             IC=ICPointInt(z_No=Fittable(minval=IC_MIN2, maxval=IC_MAX2,default=IC_DEFAULT)),
             overlay=OverlayChain(overlays=[OverlayNonDecision(
                 nondectime=Fittable(minval=ND_MIN, maxval=ND_MAX_PT)),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR_PT)

#m0int + collapsing bound (linear)
m0intlb = Model(name='m0intlb',
             drift=DriftCohInt(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                               v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2)),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundCollapsingLinear(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT),
                                              t=Fittable(minval=BOUND_T_MIN,maxval=BOUND_T_MAX)),
             IC=ICPointInt(z_No=Fittable(minval=IC_MIN2, maxval=IC_MAX2,default=IC_DEFAULT)),
             overlay=OverlayChain(overlays=[OverlayNonDecision(
                 nondectime=Fittable(minval=ND_MIN, maxval=ND_MAX)),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR)

#m0int + collapsing bound (linear) for pretoneOnly
m0intlbpt = Model(name='m0intlb',
             drift=DriftCohInt(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                               v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2)),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundCollapsingLinear(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT),
                                              t=Fittable(minval=BOUND_T_MIN,maxval=BOUND_T_MAX)),
             IC=ICPointInt(z_No=Fittable(minval=IC_MIN2, maxval=IC_MAX2,default=IC_DEFAULT)),
             overlay=OverlayChain(overlays=[OverlayNonDecision(
                 nondectime=Fittable(minval=ND_MIN, maxval=ND_MAX_PT)),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR_PT)
             
#m3 start shift + drift bias, no collapsing bound
m3 = Model(name='m3',
             drift=DriftCohBias(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                                v_High=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2),
                                v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2),      
                                v_Low=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2)),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundConstant(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT)),
             IC=ICPointPrior(z_No=Fittable(minval=IC_MIN2, maxval=IC_MAX2,default=IC_DEFAULT),
                             z_Low=Fittable(minval=IC_MIN2,maxval=IC_MAX2),
                             z_High=Fittable(minval=IC_MIN2,maxval=IC_MAX2)),
             overlay=OverlayChain(overlays=[OverlayNonDecision(
                 nondectime=Fittable(minval=ND_MIN, maxval=ND_MAX)),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR)    

#m3 start shift + drift bias, linearly collapsing bounds
m3lb = Model(name='m3lb',
             drift=DriftCohBias(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                                v_High=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2),
                                v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2),      
                                v_Low=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2)),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundCollapsingLinear(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT),
                                              t=Fittable(minval=BOUND_T_MIN,maxval=BOUND_T_MAX)),
             IC=ICPointPrior(z_No=Fittable(minval=IC_MIN2, maxval=IC_MAX2,default=IC_DEFAULT),
                             z_Low=Fittable(minval=IC_MIN2,maxval=IC_MAX2),
                             z_High=Fittable(minval=IC_MIN2,maxval=IC_MAX2)),
             overlay=OverlayChain(overlays=[OverlayNonDecision(
                 nondectime=Fittable(minval=ND_MIN, maxval=ND_MAX)),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR)

#m3dv start shift + drift bias, uniform drift rate variability
m3dv = Model(name='m3dv',
             drift=DriftUniformCohBias(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                                v_High=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2),
                                v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2),      
                                v_Low=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2),
                                v_Width=Fittable(minval=DRIFT_VAR_MIN, maxval=DRIFT_VAR_MAX)),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundConstant(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT)),
             IC=ICPointPrior(z_No=Fittable(minval=IC_MIN2, maxval=IC_MAX2,default=IC_DEFAULT),
                             z_Low=Fittable(minval=IC_MIN2,maxval=IC_MAX2),
                             z_High=Fittable(minval=IC_MIN2,maxval=IC_MAX2)),
             overlay=OverlayChain(overlays=[OverlayNonDecision(
                 nondectime=Fittable(minval=ND_MIN, maxval=ND_MAX)),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR)

#m1lb start shift, linearly collapsing bounds
m1lb = Model(name='m1lb',
             drift=DriftCohInt(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                               v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2)),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundCollapsingLinear(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT),
                                              t=Fittable(minval=BOUND_T_MIN,maxval=BOUND_T_MAX)),
             IC=ICPointPrior(z_No=Fittable(minval=IC_MIN2, maxval=IC_MAX2,default=IC_DEFAULT),
                             z_Low=Fittable(minval=IC_MIN2,maxval=IC_MAX2),
                             z_High=Fittable(minval=IC_MIN2,maxval=IC_MAX2)),
             overlay=OverlayChain(overlays=[OverlayNonDecision(
                 nondectime=Fittable(minval=ND_MIN, maxval=ND_MAX)),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR)

#m2lb drift bias, linearly collapsing bounds
m2lb = Model(name='m2lb',
             drift=DriftCohBias(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                                v_High=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2),
                                v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2),      
                                v_Low=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2)),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundCollapsingLinear(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT),
                                              t=Fittable(minval=BOUND_T_MIN,maxval=BOUND_T_MAX)),
             IC=ICPointInt(z_No=Fittable(minval=IC_MIN2, maxval=IC_MAX2,default=IC_DEFAULT)),
             overlay=OverlayChain(overlays=[OverlayNonDecision(
                 nondectime=Fittable(minval=ND_MIN, maxval=ND_MAX)),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR)

### PRETONE MODELS ###

#for all pretone models the time constant for the starting points/evidence biases is shared
tau_bias_mx = Fittable(minval=DRIFT_TAU_MIN,maxval=DRIFT_TAU_MAX)

#full pretoneOnly model with bias and adaptation terms
m14 = Model(name='m14',
             drift=DriftCohExpBias1pExpAdaptPT(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                                             v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2), 
                                             v_Bias=Fittable(minval=DRIFT_WB_MIN2, maxval=DRIFT_WB_MAX2),
                                             tau_Bias=tau_bias_mx, 
                                             wa_High=Fittable(minval=DRIFT_WA_MIN2, maxval=DRIFT_WA_MAX2),
                                             wa_Low=Fittable(minval=DRIFT_WA_MIN2, maxval=DRIFT_WA_MAX2),
                                             tau_Adapt=Fittable(minval=DRIFT_TAU_MIN,maxval=DRIFT_TAU_MAX)),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundCollapsingLinear(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT),
                                              t=Fittable(minval=BOUND_T_MIN,maxval=BOUND_T_MAX)),
             IC=ICPointPTBias(z_No=Fittable(minval=IC_LOG_MIN3, maxval=IC_LOG_MAX3,default=IC_LOG_DEFAULT),
                              z_Bias=Fittable(minval=IC_LOG_MIN3,maxval=IC_LOG_MAX3),
                              tau_Bias=tau_bias_mx),
             overlay=OverlayChain(overlays=[OverlayNonDecisionPTBiasCentered(
                 nondectime_No=Fittable(minval=ND_MIN, maxval=ND_MAX_PT),
                 nondectime_Bias=Fittable(minval=ND_BIAS_MIN,maxval=ND_BIAS_MAX),
                 tau_Bias=tau_bias_mx),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR_PT)

#same as m13ndc, but idiosyncratic biases fit for starting point
#i.e. full drift bias+adaptation, but no cue-based starting point bias
m13int = Model(name='m13int',
             drift=DriftCohExpBias1pExpAdaptPT(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                                             v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2), 
                                             v_Bias=Fittable(minval=DRIFT_WB_MIN2, maxval=DRIFT_WB_MAX2),
                                             tau_Bias=tau_bias_mx, 
                                             wa_High=Fittable(minval=DRIFT_WA_MIN2, maxval=DRIFT_WA_MAX2),
                                             wa_Low=Fittable(minval=DRIFT_WA_MIN2, maxval=DRIFT_WA_MAX2),
                                             tau_Adapt=Fittable(minval=DRIFT_TAU_MIN,maxval=DRIFT_TAU_MAX)),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundCollapsingLinear(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT),
                                              t=Fittable(minval=BOUND_T_MIN,maxval=BOUND_T_MAX)),
             IC=ICPointPTBias(z_No=Fittable(minval=IC_LOG_MIN3, maxval=IC_LOG_MAX3,default=IC_LOG_DEFAULT),
                              z_Bias=0,
                              tau_Bias=0),
             overlay=OverlayChain(overlays=[OverlayNonDecisionPTBiasCentered(
                 nondectime_No=Fittable(minval=ND_MIN, maxval=ND_MAX_PT),
                 nondectime_Bias=Fittable(minval=ND_BIAS_MIN,maxval=ND_BIAS_MAX),
                 tau_Bias=tau_bias_mx),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR_PT)

#same as m14, but no cue-based drift bias (just idiosyncratic bias + adaptation)
m19 = Model(name='m19',
             drift=DriftCohExpBias1pExpAdaptPT(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                                             v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2), 
                                             v_Bias=0,
                                             tau_Bias=0, 
                                             wa_High=Fittable(minval=DRIFT_WA_MIN2, maxval=DRIFT_WA_MAX2),
                                             wa_Low=Fittable(minval=DRIFT_WA_MIN2, maxval=DRIFT_WA_MAX2),
                                             tau_Adapt=Fittable(minval=DRIFT_TAU_MIN,maxval=DRIFT_TAU_MAX)),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundCollapsingLinear(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT),
                                              t=Fittable(minval=BOUND_T_MIN,maxval=BOUND_T_MAX)),
             IC=ICPointPTBias(z_No=Fittable(minval=IC_LOG_MIN3, maxval=IC_LOG_MAX3,default=IC_LOG_DEFAULT),
                              z_Bias=Fittable(minval=IC_LOG_MIN3,maxval=IC_LOG_MAX3),
                              tau_Bias=tau_bias_mx),
             overlay=OverlayChain(overlays=[OverlayNonDecisionPTBiasCentered(
                 nondectime_No=Fittable(minval=ND_MIN, maxval=ND_MAX_PT),
                 nondectime_Bias=Fittable(minval=ND_BIAS_MIN,maxval=ND_BIAS_MAX),
                 tau_Bias=tau_bias_mx),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR_PT)

#same as m14, but no adaptation
m17 = Model(name='m17',
             drift=DriftCohExpBias1pExpAdaptPT(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                                             v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2), 
                                             v_Bias=Fittable(minval=DRIFT_WB_MIN2, maxval=DRIFT_WB_MAX2),
                                             tau_Bias=tau_bias_mx, 
                                             wa_High=0,
                                             wa_Low=0,
                                             tau_Adapt=0),
             noise=NoiseConstant(noise=NOISE),
             bound=BoundCollapsingLinear(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT),
                                              t=Fittable(minval=BOUND_T_MIN,maxval=BOUND_T_MAX)),
             IC=ICPointPTBias(z_No=Fittable(minval=IC_LOG_MIN3, maxval=IC_LOG_MAX3,default=IC_LOG_DEFAULT),
                              z_Bias=Fittable(minval=IC_LOG_MIN3,maxval=IC_LOG_MAX3),
                              tau_Bias=tau_bias_mx),
             overlay=OverlayChain(overlays=[OverlayNonDecisionPTBiasCentered(
                 nondectime_No=Fittable(minval=ND_MIN, maxval=ND_MAX_PT),
                 nondectime_Bias=Fittable(minval=ND_BIAS_MIN,maxval=ND_BIAS_MAX),
                 tau_Bias=tau_bias_mx),
                                            OverlayUniformMixture(
                                                umixturecoef=Fittable(
                                                    minval=U_MIX_MIN, 
                                                    maxval=U_MIX_MAX))]),
             dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR_PT)



#same as m14, with prior terms for mixed blocks
m16pt5 = Model(name='m16pt5',
    drift=DriftCohExpBias1pExpAdaptPTPrior(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                                    v_No=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2),
                                    v_Low=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2), 
                                    v_High=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2), 
                                    v_Bias=Fittable(minval=DRIFT_WB_MIN2, maxval=DRIFT_WB_MAX2),
                                    tau_Bias=tau_bias_mx, 
                                    wa_High=Fittable(minval=DRIFT_WA_MIN2, maxval=DRIFT_WA_MAX2),
                                    wa_Low=Fittable(minval=DRIFT_WA_MIN2, maxval=DRIFT_WA_MAX2),
                                    tau_Adapt=Fittable(minval=DRIFT_TAU_MIN,maxval=DRIFT_TAU_MAX)),
    noise=NoiseConstant(noise=NOISE),
    bound=BoundCollapsingLinear(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT),
                                     t=Fittable(minval=BOUND_T_MIN,maxval=BOUND_T_MAX)),
    IC=ICPointPTBiasPrior(z_No=Fittable(minval=IC_LOG_MIN3, maxval=IC_LOG_MAX3,default=IC_LOG_DEFAULT),
                          z_Low=Fittable(minval=IC_LOG_MIN3, maxval=IC_LOG_MAX3),
                          z_High=Fittable(minval=IC_LOG_MIN3, maxval=IC_LOG_MAX3),
                     z_Bias=Fittable(minval=IC_LOG_MIN3,maxval=IC_LOG_MAX3),
                     tau_Bias=tau_bias_mx),
    overlay=OverlayChain(overlays=[OverlayNonDecisionPTBiasCentered(
        nondectime_No=Fittable(minval=ND_MIN, maxval=ND_MAX_PT),
        nondectime_Bias=Fittable(minval=ND_BIAS_MIN,maxval=ND_BIAS_MAX),
        tau_Bias=tau_bias_mx),
                                   OverlayUniformMixture(
                                       umixturecoef=Fittable(
                                           minval=U_MIX_MIN2, 
                                           maxval=U_MIX_MAX))]),
    dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR_PT)

#same as m16pt5 but all neutral prior intercepts fixed at 0 for fitting
#pretone_pLH
m16ptlh = Model(name='m16ptlh',
    drift=DriftCohExpBias1pExpAdaptPTPrior(v_SNR=Fittable(minval=DRIFT_MIN, maxval=DRIFT_MAX2),
                                    v_No=0,
                                    v_Low=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2), 
                                    v_High=Fittable(minval=DRIFT_BIAS_MIN2, maxval=DRIFT_BIAS_MAX2), 
                                    v_Bias=Fittable(minval=DRIFT_WB_MIN2, maxval=DRIFT_WB_MAX2),
                                    tau_Bias=tau_bias_mx, 
                                    wa_High=Fittable(minval=DRIFT_WA_MIN2, maxval=DRIFT_WA_MAX2),
                                    wa_Low=Fittable(minval=DRIFT_WA_MIN2, maxval=DRIFT_WA_MAX2),
                                    tau_Adapt=Fittable(minval=DRIFT_TAU_MIN,maxval=DRIFT_TAU_MAX)),
    noise=NoiseConstant(noise=NOISE),
    bound=BoundCollapsingLinear(B=Fittable(minval=BOUND_MIN, maxval=BOUND_MAX_ALT),
                                     t=Fittable(minval=BOUND_T_MIN,maxval=BOUND_T_MAX)),
    IC=ICPointPTBiasPrior(z_No=IC_LOG_DEFAULT,
                          z_Low=Fittable(minval=IC_LOG_MIN3, maxval=IC_LOG_MAX3),
                          z_High=Fittable(minval=IC_LOG_MIN3, maxval=IC_LOG_MAX3),
                     z_Bias=Fittable(minval=IC_LOG_MIN3,maxval=IC_LOG_MAX3),
                     tau_Bias=tau_bias_mx),
    overlay=OverlayChain(overlays=[OverlayNonDecisionPTBiasCentered(
        nondectime_No=Fittable(minval=ND_MIN, maxval=ND_MAX_PT),
        nondectime_Bias=Fittable(minval=ND_BIAS_MIN,maxval=ND_BIAS_MAX),
        tau_Bias=tau_bias_mx),
                                   OverlayUniformMixture(
                                       umixturecoef=Fittable(
                                           minval=U_MIX_MIN2, 
                                           maxval=U_MIX_MAX))]),
    dx=DX_ALT, dt=DT_ALT2, T_dur=T_DUR_PT)
