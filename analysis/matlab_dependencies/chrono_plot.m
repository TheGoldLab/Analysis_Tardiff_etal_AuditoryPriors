function [h_chrono]=chrono_plot(data,condvar,condstyle,meanvar,semvar,success,verbose,subplots,SNRvar)

if ~exist('success','var') || isempty(success)
    success = [0:1];
end

h_chrono = [];
chrono_title = {'error','correct'};

condList = unique(data.(condvar));
if ~exist('condstyle','var') || isempty(condstyle)
    condstyle = cellstr(repmat(' ',size(condList)));
end
if ischar(condstyle)
	condstyle = cellstr(repmat(condstyle,size(condList)));
end


if ~exist('verbose','var')
    verbose = false;
end

if ~exist('meanvar','var')
    meanvar = 'nanmean_nanmean_RT';
end

%legacy code butting up against model display code...
if ~exist('semvar','var') %|| isempty(semvar)
    semvar = 'nansem_nanmean_RT';
end

if ~exist('SNRvar','var')
    SNRvar = 'SNR';
end

%figure();
if ~ismember('success',data.Properties.VariableNames)
    if length(success)==1
        data.success = repmat(success,height(data),1);
        warning('No success column found. Infering from success parameter.');
    else
        error('Data has no success variable.');
    end
end

for s=success
    if length(success) > 1
        if exist('subplots','var') && ~isempty(subplots)
            set(gcf,'CurrentAxes',subplots(2-s));
        else
            subplot(2,1,2-s)
        end
    end
    hold on;
    this_h = [];
    for p=1:length(condList)
        if verbose
            sprintf('prior? %d',condList(p))
        end
        
        if iscell(condList)
            this_data = data(strcmp(data.(condvar),condList{p}) & data.success==s,:);
        else
            this_data = data(data.(condvar)==condList(p) & data.success==s,:);
        end
        
        if isempty(semvar)
            this_h = [this_h;plot(this_data.(SNRvar)(this_data.(SNRvar) <= 0),...
                this_data.(meanvar)(this_data.(SNRvar) <= 0),condstyle{p});        
                plot(this_data.(SNRvar)(this_data.(SNRvar) >= 0),...
                this_data.(meanvar)(this_data.(SNRvar) >= 0),condstyle{p})];
        else
            this_h = [this_h;errorbar(this_data.(SNRvar)(this_data.(SNRvar) <= 0),...
                this_data.(meanvar)(this_data.(SNRvar) <= 0),...
                this_data.(semvar)(this_data.(SNRvar) <= 0),condstyle{p});        
                errorbar(this_data.(SNRvar)(this_data.(SNRvar) >= 0),...
                this_data.(meanvar)(this_data.(SNRvar) >= 0),...
                this_data.(semvar)(this_data.(SNRvar) >= 0),condstyle{p})];
        end
        
    end
    %xlabel('SNR')
    %ylabel('RT')
    xlim([min(data.(SNRvar))-.11, max(data.(SNRvar))+.11])
    title(chrono_title{s+1});
    h_chrono = [this_h,h_chrono];
end