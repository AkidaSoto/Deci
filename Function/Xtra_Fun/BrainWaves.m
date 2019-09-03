function BrainWaves(ax, Frequencies, Times, Amp, AmpNoise, PhaNoise)

if length(Frequencies) ~= length(Times) || length(Frequencies) ~= length(AmpNoise)
    error('incorrect number of input lengths')
end

Timelength = minmax(cell2mat(cellfun(@minmax,Times,'UniformOutput',false)));

Signal = nan(length(Frequencies),round(1000*diff(Timelength)));

Timing = Timelength(1)+.001:.001:Timelength(2);

for i = 1:length(Frequencies)
    
    CycleTime = 1/Frequencies(i);
    
    NumCycles = diff(Times{i})/CycleTime;
    
    
    
    Sig = sin(2*pi:[pi/[1000*diff(Times{i})]]*[[[NumCycles+1]*2]-2]:[NumCycles+1]*pi*2);
    
    Sig = Sig * Amp(i);
    
    if AmpNoise(i)
        Sig = arrayfun(@(c) c*rand,Sig);
        
    end
    
    Phase = Timing >= Times{i}(1) &  Timing <= Times{i}(2);
    
    if PhaNoise(i)
        Phase = circshift(Phase,randi([0 6]));
    end
    
    Signal(i,Phase) = Sig(1:end-1);

end

plot(ax,nansum(Signal,1))
end