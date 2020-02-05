function Deci_Corr(Deci,info,freq,params)
%% Load Data
foi = freq.freq;
label = freq.label;

cfg.resamplefs = params.DownSample;
evalc('freq = ft_resampledata(cfg,freq)');
freq.fourierspctrm = permute(cell2mat(permute(freq.trial,[3 1 2])),[3 4 1 2]);
freq.freq = foi;
freq.time = freq.time{1};
freq.label = label;
freq.dimord = 'rpt_chan_freq_time';


for brains = 1:length(params.Brain)
    
    switch params.Brain{brains}
        case 'Magnitude'
            brain = abs(freq.fourierspctrm);
        case 'Phase'
            brain = angle(freq.fourierspctrm);
    end
    
    for behaviors = 1:length(params.Behavior)
        
        behavior = load([Deci.Folder.Analysis filesep 'Extra' filesep params.Behavior{behaviors} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.CondTitle{info.Cond}],params.Behavior{behaviors});
        
        Type = fieldnames(behavior);
        behavior = behavior.(Type{1});
        
        parameter = zscore(behavior);
        
        
        time_window = params.window;
        toi = find(round(freq.time,4) >= round(params.toi(1),4) & round(freq.time,4) <= round(params.toi(2),4));
        
        for foi = 1:length(freq.freq)
            
            time_window_idx = round((1000/freq.freq(foi))*time_window(foi)/(1000*mean(diff(freq.time))));
            
            for ti = 1:length(toi)
                %compute phase snychronization
                b_time = squeeze(mean(brain(:,:,foi,toi(ti)-time_window_idx:toi(ti)+time_window_idx),4));
                
                switch params.Brain{brains}
                    case 'Magnitude'
                        
                        [r,p] = corrcoef(zscore(b_time),parameter);
                        R(1,foi,ti) = r(1,2);
                        P(1,foi,ti) = p(1,2);
                        
                    case 'Phase'
                        [R(1,foi,ti),P(1,foi,ti)] =  circ_corrcl(b_time, parameter);
                        
                end
            end
        end
        
        extracorr.label = freq.label;
        extracorr.freq = freq.freq;
        extracorr.time = freq.time;
        extracorr.dimord =  'chan_freq_time';
        
        extracorr.powspctrm = R;
        R = extracorr;

        extracorr.powspctrm = P;
        P = extracorr;
        
        mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Corr' filesep params.Brain{brains} '_' params.Behavior{behaviors}  filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}])
        save([Deci.Folder.Analysis filesep 'Extra' filesep 'Corr' filesep params.Brain{brains} '_' params.Behavior{behaviors}  filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep info.Channels{info.ChanNum}],'R','P');
        clear R P
        
    end
    
end
end
%%
