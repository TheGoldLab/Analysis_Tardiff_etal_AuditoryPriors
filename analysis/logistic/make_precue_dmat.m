function dat2Fit=make_precue_dmat(data)
    dummy = ones(height(data), 1);
    dat2Fit = repmat([dummy data.zSNR],1,3); %repeat design mat per condition
    dat2Fit = [dat2Fit data.choice01]; %add choice data
    %zero out design mat columns that' not part of a given condition
    dat2Fit(data.prior==-2, [3 4 5 6]) = 0;
    dat2Fit(data.prior== 0, [1 2 5 6]) = 0;
    dat2Fit(data.prior== 2, [1 2 3 4]) = 0;
end

