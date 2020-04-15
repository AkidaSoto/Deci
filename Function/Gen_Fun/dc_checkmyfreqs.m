function [out, vec] = dc_checkmyfreqs(freqoi,endnsample,fsample)

freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1; % is equivalent to: round(freqoi .* endtime) + 1;
freqboi   = unique(freqboi);

if length(freqboi) ~=  length(freqoi)
    out = 'not good';
else
    out = 'good';
end
    
end