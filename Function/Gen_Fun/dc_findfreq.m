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
        case 'alpha2'
            freqout(freq,:) =[9 12];  %LP 4/27 added 'alpha2'
        case 'alpha3'
            freqout(freq,:) =[11 14];  %LP 4/28 added 'alpha3'
        case 'lowgamma'
            freqout(freq,:) =[30 55];
        case 'highgamma'
            freqout(freq,:) = [55 80];
        case 'gamma'
            freqout(freq,:) = [30 80];
        case 'all'
            freqout(freq,:) = [1 100];
        case 'beta2'
            freqout(freq,:) = [17 20]; % 4/27 LP added 'beta2'
        case 'lowf'
            freqout(freq,:) = [1 10]; %JC 4/29/20 added 'low'
    end
end

end