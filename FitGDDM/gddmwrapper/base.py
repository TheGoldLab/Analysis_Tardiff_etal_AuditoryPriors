# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 10:08:07 2020

@author: ntardiff
"""
import pandas as pd
import numpy as np
import pickle
import os.path
import warnings
import itertools
import copy
from ddm import Sample,Fitted
from ddm.functions import fit_adjust_model,dependence_hit_boundary,solve_all_conditions
from datetime import date
from glob import glob


def load_data(data_file,rt='rt',correct='correct',conds=None,preproc=None,
              subjnum=None,selex=None,split_subj=True,verbose=False,debug=False):
    """
    Load data from data_file and convert to format preferred for GDDM
    Can optionally specify a query string to be passed to pd.query() to only keep
    a subset of the data (i.e., individual subjects)
    subjnum - numeric value of subject, in alpha order (used for cluster arrayjob)
    Right now column names are hardcoded based on AuditoryPriors data set.
    Note:debugging does not support split_subj atm.
    """
    
    with open(data_file, "r") as f:
        data = pd.read_csv(f)
    
    #note that pyDDM can handle any name you want, just doing this
    #for consistency/ease
    data.rename(columns={rt:"rt",correct:"correct"},inplace=True) #correct names
    
    if preproc is not None:
        data = preproc(data)
       
    if selex:
        data.query(selex,inplace=True)
        
    if subjnum is not None:
        subj = data.subject.unique()
        subj = np.sort(subj)
        this_subj = subj[subjnum]
        data = data[data.subject==this_subj]
    else:
        this_subj = data.subject.unique().squeeze()
        this_subj = this_subj.item() if this_subj.ndim == 0 else this_subj.tolist()
    
    #only keep columns necessary for fitting
    keep_cols = ['correct','rt']
    if conds is not None:
        if not isinstance(conds,list):
            conds = [conds]
        keep_cols = conds + keep_cols
    
    #could this be cleaned up to be less redudant?
    if split_subj and isinstance(this_subj,list):
        sample = dict()
        for su in this_subj:
            this_data = data.loc[data.subject==su,keep_cols]
            if verbose:
                print(this_data.head(20))
            sample[su] = Sample.from_pandas_dataframe(this_data,
                                                      rt_column_name="rt", 
                                                      correct_column_name="correct")
    else:
        data = data.loc[:,keep_cols]
        if verbose:
            print(data.head(20))
        sample = Sample.from_pandas_dataframe(data,
                                              rt_column_name="rt", 
                                              correct_column_name="correct")
    
    if debug:
        return data,this_subj
    else:
        return sample,this_subj


def run_model(sample,model,subj,it=None,sess=None,out_dir='.',**kwargs):
    """
    Fit and save PyDDM model.
    """
    #build filename
    outprefix = os.path.join(out_dir,model.name)
    if it is not None:
        outprefix = outprefix + '-%s' % it
    if sess is not None:
        outprefix = outprefix + '_%s' % sess
    outprefix = outprefix + '_%s_%s' % (subj,date.today())
  
    model_fit = fit_adjust_model(sample=sample, model=model,**kwargs)
    model_fit.subject = subj
    
    save_model(model_fit,outprefix)
    get_params(model_fit,outprefix+'_params.csv')
    
    return model_fit

def save_model(model,file):
    """
    Save PyDDM model to pickle file.
    """
    with open(file,'wb') as f:
        pickle.dump(model,f)
        
def load_model(file):
    """
    Load a PyDDM pickled model.
    """
    with open(file,'rb') as f:
        model=pickle.load(f)
        
    return model

def load_models(files,verbose=False):
    """
    load multiple models.'
    """
    if isinstance(files,str):
        if os.path.isfile(files):
            files = [files]
        else:
            files = glob(files)
    
    if verbose:
        print('loading models:\n%s' % '\n'.join(files))
    
    models = [load_model(f) for f in files if os.path.isfile(f)]
    return models

def load_subjects(data_file,model_dir,model_id,selex=None,data_loader=load_data,
                  simplify=True):
    '''
    Convenience function for loading multiple subject's data as PyDDM
    samples and their associated models.
    '''
    samples,subj = data_loader(data_file,selex=selex)
    
    #paying the price for not returning a list from load_data for single subject
    #probably should change that behavior instead of doing this...
    if not isinstance(subj,list):
        subj = [subj]
    
    model_files = {su:glob(os.path.join(model_dir,model_id+'_'+su+'_*[!_params.csv]')) 
                   for su in subj}
    
    models = {k:load_models(v) for k,v in model_files.items()}
    
    if simplify:
        m_count = np.array([len(m) for m in models.values()])
        if (m_count<=1).all():
            for k in models:
                try:
                    models[k] = models[k][0]
                except IndexError:
                    pass #ignore missing models
        else:
            warnings.warn('Note: Can\'t simplify model lists.')
            
    if (m_count==0).any():
        warnings.warn('Models missing for some subjects!')
        
    if(m_count>1).any():
        warnings.warn('Multiply models returned for some subjects.')
    
    return models,samples
    

def get_params(model,outfile=None,diagnostics=False):
    params = {p : getattr(component,p).real 
              for component in model.dependencies
              for p in component.required_parameters}
    
    params = pd.DataFrame.from_dict(params,orient='index',columns=['value'])
    params.index.name = 'param'
    
    if diagnostics:
        bound_check = {}
        for component in model.dependencies:
            for p in component.required_parameters:
                pv=getattr(component,p)
                if isinstance(pv, Fitted):
                    bound_check[p] = dependence_hit_boundary(pv)
                    
        bound_check = pd.DataFrame.from_dict(bound_check,orient='index',columns=['hit_boundary'])
        params = pd.merge(
            params,bound_check,left_index=True,right_index=True,how='left')
        
    if hasattr(model,'subject'):
        params['subject'] = model.subject
    
    if outfile is not None:
        params.to_csv(outfile)
    
    return params

def get_fit_stats(model):
    fitres = model.fitresult
    assert fitres.loss=='Negative log likelihood', \
        'Only fits using negative log likelihood loss are supported.'
    nparams = fitres.properties['nparams']
    nlle = fitres.value()
    fit_stats = {
    'subject':model.subject,
    'model':model.name,
    'nparams':nparams,
    'nlle':nlle,
    #already negative lle, so + 2*nlle == -2*lle
    'bic':np.log(fitres.properties['samplesize'])*nparams + 2.*nlle,
    'aic':2.*nparams + 2.*nlle
    }
     
    return fit_stats


def get_predicted(model,sample,undec=False,forced=False,**kwargs):
    '''
    gets predicted correct and errror probabilities and mean RTs for a given
    model. Useful for plotting psychometric/chronometric functions.
    Values are based on analytical/numerical solutions to model,
    not sampling.
    
    NOTE: currently only return correct RTs.
    '''
    #solve the model for all conditions
    print('Solving model for all conditions. May take a minute...')
    sols = solve_some_conditions(model,sample,**kwargs)
    
    #extract values for each condition and produce dataframe
    soldf = {k:[] for k in list(sols.values())[0].conditions.keys()}
    soldf.update({'mean_corr':[],'mean_err':[],'mean_RT_corr':[]})
    if undec:
        soldf.update({'mean_undec':[]})
        
    for sol in sols.values():
        for k,v in sol.conditions.items():
            soldf[k].append(v.item())
            
        soldf['mean_RT_corr'].append(sol.mean_decision_time())
        if forced:
            soldf['mean_corr'].append(sol.prob_correct_forced())
            soldf['mean_err'].append(sol.prob_error_forced())
        else:
            soldf['mean_corr'].append(sol.prob_correct())
            soldf['mean_err'].append(sol.prob_error())
        
        if undec:
            soldf['mean_undec'].append(sol.prob_undecided())

    soldf = pd.DataFrame.from_dict(soldf)
    
    if hasattr(model,'subject'):
        soldf['subject']=model.subject
    else:
        warnings.warn('No subject name detected in model.')
    
    return soldf


#adapted from PyDDM code
def condition_combinations(sample, required_conditions=None,cond_filter=None,
                           cond_replace=None,cond_augment=None):
    """Get all values for set conditions and return every combination of them.

    Since PDFs of solved models in general depend on all of the
    conditions, this returns a list of dictionaries.  The keys of
    each dictionary are the names of conditions, and the value is
    a particular value held by at least one element in the sample.
    Each list contains all possible combinations of condition values.

    If `required_conditions` is iterable, only the conditions with
    names found within `required_conditions` will be included.
    
    If cond_filter is included, only specific subconditions 
    (e.g. specific SNRs) will be included in output. Only filters on conditions
    that are specified, so will keep all conditions values that are not
    filtered on. To be clear, this filters on conditions present IN THE SAMPLE.
    Will not give you conditions that don`t exist in the sample 
    (e..g SNRs that weren`t tested)
    
    To specify conditions that weren''t tested, use cond_replace or cond_augment. 
    Cond_replace will generate all combinations of conditions using the new condiitons,
    regardless of what was present in the data.
    Cond_augment will keep only the combinations of other conditions that were present
    in the data, and then combine those with the new conditions provided. Only ONE
    condition can be augmented!!!
   These different methods cannot be combined. 
   BE CAREFUL that the model is computable with the synthetic conditions.
   
   This has become super ugly/hacky but it works on at least limited cases
    """
    
    cs = sample.conditions
    conditions = []
    names = sample.condition_names()
    if required_conditions is not None:
        names = [n for n in names if n in required_conditions]
    for c in names:
        undecided = cs[c][2] if len(cs[c]) == 3 else np.asarray([])
        joined = np.concatenate([cs[c][0], cs[c][1], undecided])
        conditions.append(joined)
        
    if (cond_replace is None) and (cond_augment is None):
        alljoined = list(zip(*conditions))
        combs = list(set(alljoined))
    elif cond_replace is not None:
        assert cond_filter is None, 'Can''t use cond_filter and cond_replace in same call!'
        assert cond_augment is None, 'Can''t use cond_augment and cond_replace is same call!'
        uconds = [np.unique(x) for x in conditions]
        for k,v in cond_replace.items():
            if not isinstance(v,(list,tuple)):
                v = [v]
            uconds[names.index(k)] = v #replace conditions in sample w/ new list   
        combs=list(itertools.product(*uconds))
    else:
        assert cond_filter is None, 'Can''t use cond_filter and cond_augment in same call!'
        assert len(cond_augment)==1, 'Currently only augmenting one condition supported!'
        conds_reduced = copy.copy(conditions)
        #remove augmented conditions (have some structure here in case wanted to allow multi-augment in the future)
        aug_index = []
        aug_names = []
        aug_vals = []
        for k,v in cond_augment.items():
            conds_reduced.pop(names.index(k)) #this is why we currently only allow one argument
            aug_index.append(names.index(k))
            aug_names.append(k)
            aug_vals.append(v)
            
        alljoined = list(zip(*conds_reduced))
        combs0 = list(set(alljoined))
        combs0 = [list(x) for x in combs0] #make all unique combos of the other conds
        #now create all unique combos of the other conds w/ the new cond
        combs0 = list(itertools.product(*[*aug_vals,combs0])) 
        
        #this is why we currently only allow 1 augment--would need to do this another way
        combs=[]
        for x,y in combs0:
            z = copy.copy(y)
            z.insert(aug_index[0],x)
            combs.append(tuple(z))

    #filter out specific condition instances that we want (e.g. specific SNRs)
    if cond_filter:
        for k,v in cond_filter.items():
            if not isinstance(v,(list,tuple)):
                v = [v]
            #might be a little inefficient to filter multiple times but OK
            combs = [x for x in combs if x[names.index(k)] in v]            
    
    if len(combs) == 0: # Generally not needed since iterools.product does this
        return [{}]
    return [dict(zip(names, c)) for c in combs]


#Adapted from PyDDM code
def solve_some_conditions(model, sample, cond_filter=None, cond_replace=None, cond_augment=None, method=None):
    """Solve the model for all conditions relevant to the sample.

    This takes the following parameters:

    - `model` - A Model() object
    - `sample` - A Sample() object which has conditions for each of
      the required conditions in `model`
    - `conditions` - Restrict to specific conditions
    - `method` - A string describing the solver method.  Can be
      "analytical", "numerical", "cn", "implicit", or "explicit".

    For each value of each relevant condition in sample (i.e. those in
    the model's required conditions), this will solve the model for
    that set of parameters.  It returns a dictionary indexed by a
    frozenset of the condition names and values, with the Solution
    object as the value, e.g.:
    
        {frozenset({('reward', 3)}): <Solution object>,
         frozenset({('reward', 1)}): <Solution object>}

    This function will automatically parallelize if set_N_cpus() has
    been called.
    
    NOTE: cond_filter and cond_replace are mutually exclusive. cond_filter filters on
    conditions that exist in the sample! Cond_replace replaces a 
    condition or conditions in the sample with a unique set of values. Use carefully!
    
    """
    from ddm.functions import _parallel_pool,paranoid_settings
    
    conds = condition_combinations(sample,required_conditions=model.required_conditions,
                                              cond_filter=cond_filter,
                                              cond_replace=cond_replace,
                                              cond_augment=cond_augment)
        
    if method is None:
        meth = model.solve
    elif method == "analytical":
        meth = model.solve_analytical
    elif method == "numerical":
        meth = model.solve_numerical
    elif method == "cn":
        meth = model.solve_numerical_cn
    elif method == "implicit":
        meth = model.solve_numerical_implicit
    elif method == "explicit":
        meth = model.solve_numerical_explicit
    else:
        raise ValueError("Invalid method "+method)

    cache = {}
    if _parallel_pool is None: # No parallelization
        for c in conds:
            cache[frozenset(c.items())] = meth(conditions=c)
        return cache
    else: # Parallelize across pool
        if paranoid_settings.get('enabled') is False:
            # The *2 makes sure that this runs on all subprocesses,
            # since you can't broadcast commands to all processes
            _parallel_pool.map(lambda x : paranoid_settings.set(enabled=False), [None]*_parallel_pool.n_cpus*2)
        sols = _parallel_pool.map(meth, conds, chunksize=1)
        for c,s in zip(conds, sols):
            cache[frozenset(c.items())] = s
        return cache
