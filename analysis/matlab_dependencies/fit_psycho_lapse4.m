function [fits_lapse,fit_stats] = fit_psycho_lapse4(data,fit_options) %optionName,max_ptlen)
%% code for fitting logistic function w/ lapse rate for AuditoryPriors

addpath('../logistic')

%get subjects
subj = unique(data.subject);

%% now with fit lapse
fits_lapse = table();
fit_stats = table();
lapse_vars_full = [fit_options.base_params,'lapse'];
stats_vars = {'LL','dev','adev','p'};

for s=1:length(subj)
    this_subj = subj(s);
    
    %pull out subject data %and all necessary variables
    this_datamat = data(strcmp(data.subject,this_subj),:); %...
        %fit_options.dmat_invars);

    %construct design matrix (each two columns are intercept/slope for a
    %condition, and last column is choice
    dat2Fit = fit_options.dmat_func(this_datamat);
    
    % do fit
    [this_fit,~,this_stats] = logist_fit(dat2Fit);
    
    %construct table from params and fit stats
    this_fit = array2table(this_fit','VariableNames',lapse_vars_full);
    this_stats = array2table(this_stats,'VariableNames',stats_vars);
    %this_fit = [this_fit this_stats];
    %this_fit.nparam = length(lapse_vars_full);
    this_stats.nparam = length(lapse_vars_full);
    this_stats.subject = this_subj;
    this_fit.subject = this_subj;
    
    metavars = {'mname','model_id','max_ptlen'};
    for mv = 1:length(metavars)
        if isfield(fit_options,metavars{mv})
            this_meta = fit_options.(metavars{mv});
            if ischar(this_meta)
                this_meta = cellstr(this_meta);
            end
            this_fit.(metavars{mv}) = repmat(this_meta,1);
            this_stats.(metavars{mv}) = repmat(this_meta,1);
        end
    end
    %{
    if isfield(fit_options,'mname')
        this_fit.mname = repmat(fit_options.mname,1);
        this_stats.mname = repmat(fit_options.mname,1);
    end
    if isfield(fit_options,'model_id')
        this_fit.model_id = repmat(fit_options.model_id,1);
        this_stats.model_id = repmat(fit_options.model_id,1);
    end
    if isfield(fit_options,'max_ptlen')
        this_fit.max_ptlen = repmat(fit_options.max_ptlen,1);
        this_stats.max_ptlen = repmat(fit_options.max_ptlen,1);
    end
    %}
    fits_lapse = [fits_lapse;this_fit];
    fit_stats = [fit_stats;this_stats];
    
    %clear this_datamat dat2Fit this_fit this_subj
end


end