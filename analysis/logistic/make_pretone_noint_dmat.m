function dat2Fit=make_pretone_noint_dmat(data,ptlen,omit_pt)
    dummy = ones(height(data), 1);
    
    ptmat = ragged_array(data.pretones,ptlen);
    if exist('omit_pt','var') && ~isempty(omit_pt)
        ptmat(:,omit_pt) = [];
    end
    ptmat(isnan(ptmat)) = 0; %zero out pretone positions w/ no pretones
    dat2Fit = [dummy data.zSNR ptmat data.choice01];

end