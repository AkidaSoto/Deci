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
            
            if params.bsl.do
                Bsl = brain;
                toi = freq.time >= round(params.bsl.time(1),4) & freq.time <= round(params.bsl.time(2),4);
                
                    Bsl = nanmean(Bsl(:,:,:,toi),4);
                    Bsl = repmat(Bsl,[1 1 size(Bsl,4)]);
                
                switch params.bsl.type
                    case 'none'
                    case 'absolute'
                        brain =  brain - Bsl;
                    case 'relative'
                        brain=  brain ./ Bsl;
                    case 'relchange'
                        brain = ( brain - Bsl) ./ Bsl;
                    case 'db'
                        brain = 10*log10( brain ./ Bsl);
                end
                               
            end
            
        case 'Phase'
            brain = angle(freq.fourierspctrm);
    end
    
    

    
    
    
    for behaviors = 1:length(params.Behavior)
        
        behavior = load([Deci.Folder.Analysis filesep 'Extra' filesep params.Behavior{behaviors} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.CondTitle{info.Cond}]);
        
        Type = fieldnames(behavior);
        

        behavior = behavior.(Type{1});   
        parameter = zscore(behavior);
        
        
        time_window = params.window;
        toi = find(round(freq.time,4) >= round(params.toi(1),4) & round(freq.time,4) <= round(params.toi(2),4));
        
        for foi = 1:length(freq.freq)
            
            time_window_idx = round((1000/freq.freq(foi))*time_window(foi)/(1000*mean(diff(freq.time))));
            
            for ti = 1:length(toi)
                %compute phase snychronization
                
                if params.slidingwindow
                b_time = squeeze(mean(brain(:,:,foi,toi(ti)-time_window_idx:toi(ti)+time_window_idx),4));
                else
                b_time =  brain(:,:,foi,ti);
                end
                
                if params.zscorebrain
                   b_time = zscore(b_time); 
                end
                
                switch params.Brain{brains}
                    case 'Magnitude'
                        
                        params = Exist(params,'type',[]);
                        if strcmpi(params.type,'spearman')   
                        [R(1,foi,ti), P(1,foi,ti)] = corr(b_time,parameter,'Type','Spearman');
                        else
                        [r,p] = corrcoef(b_time,parameter);
                        R(1,foi,ti) = r(1,2);
                        P(1,foi,ti) = p(1,2);
                        end
                        
                        
                    case 'Phase'
                        [R(1,foi,ti),P(1,foi,ti)] =  circ_corrcl(b_time, parameter);
                        
                end
            end
        end
        
%         if any(any(arrayfun(@(c) iscomplex(c),R)))
%            k = 0; 
%         end
        
        extracorr.label = freq.label;
        extracorr.freq = freq.freq;
        extracorr.time = freq.time(toi);
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
