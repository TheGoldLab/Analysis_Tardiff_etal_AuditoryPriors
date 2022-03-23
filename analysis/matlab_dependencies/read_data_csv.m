function [datapt] = read_data_csv(file)
%general loading/light processing of AuditoryPriors data from csv files
%generated from analyze_priorOnly.m. Used in displaying/processing data for
%jupyter notebooks analysis.

datapt = readtable(file);

try
    datapt = datapt(strcmp(datapt.optionName,'priorOnly') | datapt.pretoneLength>0,:);
    datapt.bias = pt2num(datapt.pretoneSeq); %recode pretones for PyDDM format
    %create categories for plotting
    datapt.ptC = cellfun(@(x) fliplr(x(1:end-1)),datapt.pretoneSeq,'UniformOutput',false);
    datapt.ptC = cellfun(@(x) x(1:min(2,length(x))),datapt.ptC,'UniformOutput',false);
    datapt.ptC = categorical(datapt.ptC);
catch ME
    warning(ME.message);
    warning('No pretone data detected. Only reading in file.');
end

end
