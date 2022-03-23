function dat2Fit=make_basic_dmat(data)
    dummy = ones(height(data), 1);
    dat2Fit = [dummy data.zSNR data.choice01]; %add choice data
end
