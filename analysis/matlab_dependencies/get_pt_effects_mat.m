function [pt_total_fx] = get_pt_effects_mat(ptfits,zaSNR,max_ptlen,aSNR) %,SNR_leg) %,plotfx)
%used to calculate effects of each pretone from Matlab psychometric
%function fitting of all data. See get_pt_effects for deriving from mixed
%models fit on each pretone length.

%{
if ~exist('plotfx','var')
    plotfx=true;
end
%}
if ~exist('aSNR','var')
    aSNR = zaSNR;
end


pt_total_fx = table();

subj = unique(ptfits.subject);
main_effect_cols = arrayfun(@(x) sprintf('pt%d',x),[1:max_ptlen],'UniformOutput',false);
main_effect_cols = ismember(ptfits.Properties.VariableNames,main_effect_cols);
interaction_effect_cols = arrayfun(@(x) sprintf('pt%d_aSNR',x),[1:max_ptlen],'UniformOutput',false);
interaction_effect_cols = ismember(ptfits.Properties.VariableNames,interaction_effect_cols);
    
for i=1:length(subj)
    this_fit = ptfits(strcmp(ptfits.subject,subj{i}),:);
    main_fx = this_fit{:,main_effect_cols};
    interaction_fx= this_fit{:,interaction_effect_cols} .* zaSNR;
    total_fx = interaction_fx + main_fx;
    subject = repmat(subj(i),size(total_fx,1),1);
    %pretone_pos = (1:size(total_fx,2))';
    this_table = table(subject,aSNR,total_fx);
    pt_total_fx = [pt_total_fx;this_table];
end

%{
if plotfx
    for i=pt_lengths
        figure();hold on;
        this_data = pt_total_fx(pt_total_fx.pretoneLength==i,:);
        for s=1:size(this_data.total_fx,2)
            plot(this_data.pretone_pos,this_data.total_fx(:,s));
        end
        xticks(this_data.pretone_pos);
        title(sprintf('length %d',i))
        xlabel('pretone pos');
        ylabel('log odds ratio');
        yline(0,'--');
	    legend(arrayfun(@num2str,SNR_leg,'UniformOutput',false));
    end
end
%}

end