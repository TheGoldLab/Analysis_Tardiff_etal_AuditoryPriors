%% set options and read in data

optionName = 'pretoneOnly'; %'priorOnly';
model_id = 4;
max_ptlen = 6;
omit_pt = [];

DATA_DIR = './data';
DATA_FILE = '/priorOnly_pretoneOnly_pretone_pL_pretone_pH_pretone5_05-Feb-2021.csv';

%read in data
data = read_data_csv(...
    fullfile(DATA_DIR,DATA_FILE));

%% model setup
prior_labels = {'low','no','high'};
params = {'zbias_','zSNR_'};
switch model_id
    case 1 %full precue model w/ separate intercepts/slopes for each cue
        fit_options.base_params = [strcat(params,prior_labels{1}) strcat(params,prior_labels{2}) ...
            strcat(params,prior_labels{3})];
        fit_options.dmat_func = @make_precue_dmat;
        fit_options.mname = 'priorOnly_full';
        fit_options.model_id = model_id;
        %fit_options.dmat_invars = {'prior','zSNR','choice01'};
    case 2 %precue shared slopes model
        %nomenclature between bias/zbias getting confusing here...bias
        %means adjusted for any zSNR interaction...doesn't mean not zscore
        %scale...
        fit_options.base_params = {'bias_low','bias_no','bias_high','zSNR'};
        fit_options.dmat_func = @make_precue_noint_dmat;
        fit_options.mname = 'priorOnly_noint';
        fit_options.model_id = model_id;
    case 3 %base model, no bias/adapt
        fit_options.base_params = {'bias','zSNR'};
        fit_options.dmat_func = @make_basic_dmat;
        if strcmp(optionName,'priorOnly')
            fit_options.mname = 'priorOnly_noprior';
        elseif strcmp(optionName,'pretoneOnly')
            fit_options.mname = 'pretoneOnly_nopt';
        end
        fit_options.model_id = model_id;
    case 4 % full pretone model
        pt_lens = arrayfun(@num2str,1:max_ptlen,'UniformOutput',false);
        if ~isempty(omit_pt)
            pt_lens(omit_pt) = [];
            fit_options.mname = sprintf('pretoneOnly_full%do%d',max_ptlen,omit_pt);
        else
            fit_options.mname = sprintf('pretoneOnly_full%d',max_ptlen);
        end
        %not sure if this should be zbias or bias for consistent
        %nomenclature...think about it.
        fit_options.base_params = ['zbias','zSNR',strcat('pt',pt_lens),strcat('pt',pt_lens,'_aSNR')];
        fit_options.max_ptlen = max_ptlen;
        fit_options.omit_pt = omit_pt;
        fit_options.dmat_func = @(x) make_pretone_dmat(x,max_ptlen,omit_pt);
        fit_options.model_id = model_id;
        %fit_options.dmat_invars = {'zSNR','aSNR','pretones','choice01'};
    case 5 %pretone w/o adaptation terms
        pt_lens = arrayfun(@num2str,1:max_ptlen,'UniformOutput',false);
        if ~isempty(omit_pt)
            pt_lens(omit_pt) = [];
            fit_options.mname = sprintf('pretoneOnly_noint%do%d',max_ptlen,omit_pt);
        else
            fit_options.mname = sprintf('pretoneOnly_noint%d',max_ptlen);
        end
        %not sure if this should be zbias or bias for consistent
        %nomenclature...think about it.
        fit_options.base_params = ['zbias','zSNR',strcat('pt',pt_lens)];
        fit_options.max_ptlen = max_ptlen;
        fit_options.omit_pt = omit_pt;
        fit_options.dmat_func = @(x) make_pretone_noint_dmat(x,max_ptlen,omit_pt);
        fit_options.model_id = model_id;
    case 6 %precue w/ "adaptation" terms
        params = {'zbias_','aSNR_'};
        fit_options.base_params = [strcat(params,prior_labels{1}) strcat(params,prior_labels{2}) ...
            strcat(params,prior_labels{3}) 'zSNR'];
        fit_options.dmat_func = @make_precue_adapt_dmat;
        fit_options.mname = 'priorOnly_adapt';
        fit_options.model_id = model_id;
end
%% data setup

datapc = data(strcmp(data.optionName,optionName),:);
%remove missing trials (should already be done)
datapc = datapc(~isnan(datapc.choice01),:);

if strcmp(optionName,'pretoneOnly')
    datapc(datapc.pretoneLength==0,:) = [];
    datapc.pretones = cellfun(@(x) fliplr(replace(x,{'X','H','L'},{'','1','0'})), ...
        datapc.pretoneSeq,'UniformOutput',false);
    %convert to numeric and recode as 0.5/-0.5
    datapc.pretones = cellfun(@(x) (arrayfun(@str2double,x)-0.5),datapc.pretones,'UniformOutput',false);
end

if strcmp(optionName,'pretoneOnly') || model_id==6
    datapc.aSNR = zscore(abs(datapc.SNR));
    
    norms.aSNR.mean = mean(abs(datapc.SNR));
    norms.aSNR.sd = std(abs(datapc.SNR));
    norms.aSNR.zero = -norms.aSNR.mean/norms.aSNR.sd;
end

%zscore SNR (remember this is w/o missing trials)
datapc.zSNR = zscore(datapc.SNR);
norms.SNR.mean = mean(datapc.SNR);
norms.SNR.sd = std(datapc.SNR);
norms.SNR.zero = -norms.SNR.mean/norms.SNR.sd; %need SNR 0 for calculating intercepts in presence of interaction w/ SNR

%% do fit 
[fits_lapse,fit_stats] = fit_psycho_lapse4(datapc,fit_options);

%check for boundary issues
fit_params = [fit_options.base_params,'lapse'];
min_param = varfun(@min,fits_lapse,'InputVariables',fit_params);
max_param = varfun(@max,fits_lapse,'InputVariables',fit_params);

if any(min_param{:,end-1} <= -299 | max_param{:,end-1} >= 299 | max_param.max_lapse>=.44)
    warning('parameter boundary hit!');
    pause(2)
end

figure();
for i=1:min(length(fit_params),16)
    subplot(4,4,i);
    histogram(fits_lapse.(fit_params{i}),25)
    title(fit_params{i});
end

near_bound=sum(abs(abs(fits_lapse{:,1:length(fit_options.base_params)})-300) <= 1);
near_bound(end+1) = sum(abs(abs(fits_lapse{:,length(fit_options.base_params)+1})-.45) <= .02);

figure();bar(near_bound);
xticklabels([fit_options.base_params,'lapse']);
xtickangle(45)


%% save output

if strcmp(optionName,'priorOnly')
    %now compute bias for SNR0
    prior_suffix = {'low','no','high'};
    for p=1:length(prior_suffix)
        if model_id==1
        fits_lapse.(['bias_' prior_suffix{p}]) = fits_lapse.(['zbias_' prior_suffix{p}]) + ...
            fits_lapse.(['zSNR_' prior_suffix{p}]).*norms.SNR.zero;
        end
    end
    
    save(fullfile(DATA_DIR,sprintf('fits_lapse_%s_%s.mat',fit_options.mname,date())),...
        'fits_lapse','norms','fit_stats');

elseif strcmp(optionName,'pretoneOnly')
    aSNRdc = unique(datapc.aSNR);
    aSNRdc0 = [norms.aSNR.zero;aSNRdc];
    aSNR0 = [0;unique(abs(datapc.SNR))];
    
    if model_id==4 && isempty(omit_pt)
        pt_total_fx = get_pt_effects_mat(fits_lapse,aSNRdc0,max_ptlen,aSNR0);
        writetable(pt_total_fx,fullfile(DATA_DIR,sprintf('pt_total_fx_%s_%s.csv',...
            fit_options.mname,date())));
        save(fullfile(DATA_DIR,sprintf('fits_lapse_%s_%s.mat',fit_options.mname,date())),...
            'fits_lapse','norms','pt_total_fx','fit_stats');
    else
        save(fullfile(DATA_DIR,sprintf('fits_lapse_%s_%s.mat',fit_options.mname,date())),...
            'fits_lapse','norms','fit_stats');
    end

end
writetable(fits_lapse,fullfile(DATA_DIR,sprintf('fits_lapse_%s_%s.csv',...
    fit_options.mname,date())));
writetable(fit_stats,fullfile(DATA_DIR,sprintf('fit_stats_%s_%s.csv',...
    fit_options.mname,date())));