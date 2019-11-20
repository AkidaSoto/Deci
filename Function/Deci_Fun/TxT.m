function TxT(Deci,info,freq,params)

k = 0;


if ismember(info.Channels{info.ChanNum},params.Channels)

mkdir([Deci.Folder.Analysis filesep 'Txt_avg_FreqBand_pow' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep info.Channels{info.ChanNum}]);

    
    for fh = 1:length(params.freq)
        switch params.freq{fh}
            case 'theta'
                HF{fh} = [4 8];
            case 'beta'
                HF{fh} = [12.5 30];
            case 'alpha'
                HF{fh} =[8 12.5];
            case 'lowgamma'
                HF{fh} =[30 60];
%             case 'highgamma'
%                 HF{fh} = [55 80];
        end
    end
    
    for f =  1:length(params.freq)
        cfg = [];
        cfg.frequency = [HF{f}];
        frequency = ft_selectdata(cfg,freq);
        trial_freq_time = squeeze(frequency.fourierspctrm);
        frequency.powspctrm = abs(trial_freq_time).^2 ;
        frequency.avg_FreqBand_pow = squeeze(mean(frequency.powspctrm,2));
        
        frequency = rmfield(frequency,'fourierspctrm');
        frequency = rmfield(frequency,'powspctrm');
        frequency = rmfield(frequency,'cumtapcnt');
        
        save([Deci.Folder.Analysis filesep 'Txt_avg_FreqBand_pow' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep info.Channels{info.ChanNum} filesep params.freq{f}],'frequency','-v7.3');
    end
end
end