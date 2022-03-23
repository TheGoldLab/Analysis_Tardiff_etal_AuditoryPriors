function [lm]=lm_adjust(lm,group_correct)
%function for performing multiple comparisons correction for auditorypriors
%pupil regressions. Correction is performed per variable across time,
%separately for correct and incorrect models. If group_correct is passed,
%these variables are corrected together across variables and time (e.g. for
%pairwise contrasts).

if ~exist('group_correct','var')
    group_correct = {};
end

all_params = unique(lm.param);
lm.p_corr = nan(height(lm),1);

if ismember('success',lm.Properties.VariableNames)
    success = unique(lm.success);
    for s=1:length(success)
        lm.p_corr(ismember(lm.param,group_correct) & ...
            lm.success==success(s)) = ...
            pval_adjust(lm.p(ismember(lm.param,group_correct) & ...
            lm.success==success(s)),'fdr');
        for p=1:length(all_params)
            if ~ismember(all_params{p},group_correct)
                lm.p_corr(strcmp(lm.param,all_params{p}) & ...
                    lm.success==success(s)) = ...
                    pval_adjust(lm.p(strcmp(lm.param,all_params{p}) & ...
                    lm.success==success(s)),'fdr');
            end
        end
    end
else
    lm.p_corr(ismember(lm.param,group_correct)) = ...
            pval_adjust(lm.p(ismember(lm.param,group_correct)),'fdr');
    for p=1:length(all_params)
        if ~ismember(all_params{p},group_correct)
            lm.p_corr(strcmp(lm.param,all_params{p})) = ...
                pval_adjust(lm.p(strcmp(lm.param,all_params{p})),'fdr');
        end
    end
end

end
