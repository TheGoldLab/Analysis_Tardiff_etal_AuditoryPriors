function dat2Fit=make_precue_noint_dmat(data)
    dummy = ones(height(data), 3); %repeat intercept for each condition
    dat2Fit = [dummy data.zSNR data.choice01]; %add choice data
    %zero out design mat columns that not part of a given condition
    dat2Fit(data.prior==-2, [2 3]) = 0;
    dat2Fit(data.prior== 0, [1 3]) = 0;
    dat2Fit(data.prior== 2, [1 2]) = 0;
end

