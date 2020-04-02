function freqout = dc_findfreq(freqin)

for freq = 1:length(freqin)
    switch freqin{freq}
        case 'delta'
            freqout(freq,:) = [1 4]; %LP 4/2/20 added 'delta' [1 4] 
        case 'theta'
            freqout(freq,:) = [4 8];
        case 'lowbeta'
            freqout(freq,:) = [12.5 21];
        case 'highbeta'
            freqout(freq,:) = [21 30];
        case 'beta'
            freqout(freq,:) = [12.5 30];
        case 'alpha'
            freqout(freq,:) =[8 14]; %LP 4/2/20 changed from [8 12.5] to [8 14] 
        case 'lowgamma'
            freqout(freq,:) =[30 55];
        case 'highgamma'
            freqout(freq,:) = [55 80];
        case 'all'
            freqout(freq,:) = [1 100];
    end
end

end