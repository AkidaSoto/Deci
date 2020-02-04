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
            brain = phase(freq.fourierspctrm);
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
            
            for ti = 1:length(freq.time)
                %compute phase snychronization
                b_time = mean(brain(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx),4);
                
                switch params.Brain{brains}
                    case 'Magnitude'
                        
                        [r,p] = corrcoef(zscore(b_time),parameter);
                        R(fois,tois) = r(1,2);
                        P(chns,fois,tois) = p(1,2);
                        
                    case 'Phase'
                        [R(chns,fois,tois),P(chns,fois,tois)] =  circ_corrcl(Subjects{subject_list,Conditions}.powspctrm(:,chns,fois,tois), parameter);
                        
                end
                
            end
        end
        
        extracorr.label = Subjects{subject_list,Conditions}.label;
        extracorr.freq = Subjects{subject_list,Conditions}.freq;
        extracorr.time = Subjects{subject_list,Conditions}.time;
        extracorr.dimord =  'chan_freq_time';
        
        extracorr.powspctrm = R;
        R = extracorr;
        
        RSub{subject_list,Conditions,Var}= R;
        
        extracorr.powspctrm = P;
        P = extracorr;
        
        mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Corr' filesep params.Freq{brains} '_' params.Behavior{behaviors}  filesep Deci.SubjectListinfo{subject_list}  filesep Deci.Analysis.LocksTitle{infoLock} filesep Deci.Analysis.CondTitle{Cond}])
        save([Deci.Folder.Analysis filesep 'Extra' filesep 'Corr' filesep params.Freq{brains} '_' params.Behavior{behaviors}  filesep Deci.SubjectListinfo{subject_list}  filesep Deci.Analysis.LocksTitle{infoLock} filesep Deci.Analysis.CondTitle{Cond} filesep info.Channels{info.ChanNum}],'R','P');
        
        
    end
    
end
end
%%
