function [out, vec] = dc_checkmyfreqs(freqs,length,samplingrate)

vec = unique(freqs/(samplingrate / length));

if length(vec) ~=  length(freqs)
    out = 'not good';
else
    out = 'good';
end
    
ende