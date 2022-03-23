function [pcs,subj_z,subj_v,zs,mes]=sweep2mat(param_file,sweep_dir,bias_type,zrange,merange)
%convenience function for packaging PyDDM parameter sweeps identically to
%the analytical sweeps produced by tradeoff_contour_analysis w/ the persubj
%flag
% input:
%   param_file - recovered params from actual fitting (for convenience in
%       returning the subject's true params
%   sweeps_dir - directory containing all param sweep files (should only
%   contain files for ONE sweep (though sweep can be spread across multiple
%   files)
%  bias_type - whether priorOnly (1) or pretoneOnly (2)
% zrange,merange - restrict output to subset of sweep range

debug = false;

if ~exist('zrange','var')
    zrange = [];
end
if ~exist('merange','var')
    merange = [];
end

mparams = readtable(param_file);
mparams.hit_boundary = [];
mparamsW = unstack(mparams,'value','param');
subj = sort(unique(mparams.subject));

if bias_type==1
    priorCond = [-2;0;2];
    priorLabel = {'Low';'No';'High'};
    biasvar = 'prior';
elseif bias_type==2
    priorCond = [11;22];
    pt_vector.LL = [-1 -1];
    pt_vector.HH = [1 1];
    priorLabel = {'LL';'HH'};
    biasvar = 'bias';
end

pcs = [];
for s=1:length(subj)
    subj_files = dir(fullfile(sweep_dir,['pcs_*' subj{s} '_*']));
    sweeps = [];
    for f=1:length(subj_files)
        sweeps = [sweeps;readtable(fullfile(subj_files(f).folder,subj_files(f).name),...
            'Format','%f%f%s%f%f')];
    end
    try
        sweeps = sortrows(sweeps,{'z','v'},'ascend');
    catch ME
        if isempty(sweeps)
            warning('subject %s sweep is missing!',subj{s});
            continue;
        else
            error(ME.message);
        end
    end
    assert(length(unique(sweeps.subject))==1,'Multiple subjects in single file!');
    if s==1
        %full sweep range
        zs_all = sweeps.z;
        mes_all = sweeps.v;
        
        %unique sweep values
        zs = unique(zs_all);
        mes = unique(mes_all);
        
        if ~isempty(zrange)
            zs = zs(zs >= zrange(1) & zs <= zrange(2));
        end
        if ~isempty(merange)
            mes = mes(mes >= merange(1) & mes <= merange(2));
        end
    else
        assert(isequal(sweeps.z,zs_all) && isequal(sweeps.v,mes_all),...
            'Inconsistent sweep ranges across subject!');
    end
    for p=1:length(priorCond)
        this_mcorr = sweeps.mean_corr(sweeps.(biasvar)==priorCond(p) & ...
            ismember(sweeps.z,zs) & ismember(sweeps.v,mes));
        %need to transpose so zs are rows
        pcs.(priorLabel{p}).(subj{s}) = ...
            reshape(this_mcorr,...
            length(mes),length(zs))';
        
        %debugging
        if debug
            for zz=1:length(zs)
                for mm=1:length(mes)
                    assert(isequal(pcs.(priorLabel{p}).(subj{s})(zz,mm),...
                        sweeps.mean_corr(sweeps.(biasvar)==priorCond(p) & ...
                        sweeps.z==zs(zz) & sweeps.v==mes(mm))),'Bad reshape!')
                end
            end
        end
    end
end

%package subject fits for output
for p=1:length(priorLabel)
    if bias_type==1
        subj_z.(priorLabel{p}) = mparamsW.(['z_' priorLabel{p}]);
        subj_v.(priorLabel{p}) = mparamsW.(['v_' priorLabel{p}]);
    elseif bias_type==2
        %{
        subj_z.(priorLabel{p}) = mparamsW.z_Bias.*...
            pt_exp(mparamsW.tau_Bias,pt_vector.(priorLabel{p}));
        subj_v.(priorLabel{p}) = mparamsW.v_Bias.*...
            pt_exp(mparamsW.tau_Bias,pt_vector.(priorLabel{p}));
        %}
        subj_z.(priorLabel{p}) = mparamsW.z_Bias;
        subj_v.(priorLabel{p}) = mparamsW.v_Bias;
    end
end

end