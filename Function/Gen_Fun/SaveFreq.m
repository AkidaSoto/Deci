function SaveFreq(Fourier,Subject,Cond,Channels)

Fourier.freq = Fourier.oldfoi;

if ischar(Channels)
    Chan = Fourier.label;
elseif all(ismember(Channels,Fourier.label))
    Chan = Channels;
else
    error('Wrong Channel Selection in Analyis');
end

for i = 1:length(Chan)
    
    dcfg = [];
    dcfg.channel = Chan(i);
    freqplaceholder = ft_selectdata(dcfg,Fourier);
    
    mkdir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Subject filesep Cond]);
    mkdir([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep Subject filesep Cond]);
    
    label = freqplaceholder;
    label = rmfield(label,'fourierspctrm');
    label.label = Chan;
    label.dimord = 'chan_freq_time';
    
    freq = freqplaceholder;
    freq.dimord = 'chan_freq_time';
    freq.powspctrm      = permute(abs(mean(freq.fourierspctrm./abs(freq.fourierspctrm),1)),[2 3 4 1]);         % divide by amplitude
    freq  = rmfield(freq,'fourierspctrm');
    save([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep Subject filesep Cond filesep Chan{i}],'freq','label','-v7.3');
    
    freq = freqplaceholder;
    freq.powspctrm = permute(mean(abs(freq.fourierspctrm).^2 ,1),[2 3 4 1]);
    freq.dimord = 'chan_freq_time';
    freq  = rmfield(freq,'fourierspctrm');
    save([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Subject filesep Cond filesep Chan{i}],'freq','label','-v7.3');
    
end

end