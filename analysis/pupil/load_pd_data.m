function [out,data,metadata,baselines,task_data,ID2subj,validTrials,validTrialsBL] = ...
	load_pd_data(data_file,meta_file,bl_file,valid_file,task_data,taskvars,rmNaN)

if ~exist('rmNaN','var')
    rmNaN = true;
end

load(valid_file);
load(bl_file);
FM = load(meta_file);
F = load(data_file);
if isfield(F,'datads')
    data = F.datads;
    metadata = FM.metadatads;
else
    data = F.data;
    metadata = FM.metadata;
end
clear F FM

%add dataID to task_data for merging later (this will remove subjects we
%don't have pupil for)
ID2subj = table([metadata.subject]',[metadata.dataID]','VariableNames',{'subject','dataID'});
task_data = innerjoin(task_data,ID2subj,'Keys','subject');

%verify data match
to_verify = {'trialStart','trialEnd','stimOn','response';'trialStarts',...
    'trialStops','stimStarts','postResponse'};
[badtt, vtrials] = verifyTimeStamps(metadata,task_data,to_verify);

assert(isempty(badtt),'Pupil-task data mismatch detected!');

%now add eyelink timestamps to task data in case needed
task_data = add_mtimes(metadata,task_data);

%now remove trials w/ no response (or the weird cases of no recorded
%response despite RT
task_data(isnan(task_data.choice01) | isnan(task_data.choiceTimeStamp) | ...
    task_data.postResponse==0,:) = [];

%% data cleanup
%do not want any data before fixation or long periods between fix and stimOn
data = data(data.trial_time_fixOn >= 0 & data.trial_time_stimOn >= -275,:);

%remove trials deemed invalid by QA check
data = data(ismember(data(:,{'dataID','trialN'}),validTrials,'rows'),:);
baselines = baselines(ismember(baselines(:,{'dataID','trialN'}),validTrialsBL,'rows'),:);

%finally remove all nans if requested
if rmNaN
    data(isnan(data.pupilCblz),:) = [];
    %baselines(isnan(baselines.pupilBL),:) = []; %should be handled above
end

%confirm deletions
assert(isequal(validTrials,unique(data(:,{'dataID','trialN'}),'rows')),...
    'pupil data improperly removed!')

assert(isequal(validTrialsBL,baselines(:,{'dataID','trialN'})),...
    'pupil data improperly removed!');


%% merge pupil and task data, including task vars and timing offsets 
%note: this will remove trials w/o response from pupil data!

%now add task variables
data = innerjoin(data,task_data,'LeftKeys',{'dataID','trialN'},'RightKeys',{'dataID','trialID'},...
    'RightVariables',taskvars);

%confirm no trials remain for which time stamps were not verified
vsubj = unique(vtrials.dataID);
for i=1:length(vsubj)
    this_trialN = unique(data.trialN(data.dataID==vsubj(i)));
    for j=1:size(to_verify,2)
        assert(all(ismember(this_trialN,vtrials.trialID(vtrials.dataID==vsubj(i) & ...
            strcmp(vtrials.event,to_verify(1,j))))))
    end
    clear this_trialN
end

%also merge baseline data
baselines = innerjoin(baselines,task_data,'LeftKeys',{'dataID','trialN'},'RightKeys',{'dataID','trialID'},...
    'RightVariables',taskvars);

out = struct('data',data,'metadata',metadata,'baselines',baselines,'task_data',task_data,...
    'ID2subj',ID2subj,'validTrials',validTrials,'validTrialsBL',validTrialsBL);

end