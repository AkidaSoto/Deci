function Analyzor(Deci)

if ~isfield(Deci.Analysis,'ArtifactReject')
    Deci.Analysis.ArtifactReject = 0;
    Warning('Parameter for ArtifactReject not found, presuming not wanted')
end

if ~isfield(Deci.Analysis,'Channels')
    Deci.Analysis.Channels = 'all';
    Warning('Parameter for Channels not found, presuming all')
end


if ~isfield(Deci.Analysis,'Laplace')
    Deci.Analysis.Laplace = 0;
    Warning('Parameter for Laplace not found, presuming not wanted')
end


if  isempty(Deci.Analysis.Freq) &&  ~Deci.Analysis.ERP
    error('No analysis step was called for.')
end


for subject_list = 1:length(Deci.SubjectList)
    
    data = load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list}]);
    label = data.label;
    data = data.data;
    
    
    if ~isempty(Deci.Analysis.Laplace)
        
        Eyecfg.channel = Deci.Analysis.Laplace.EyeChans;
        nonscalp_data = ft_preprocessing(Eyecfg,data);
        
        [elec.label, elec.elecpos] = elec_1020select(data.label);
        ecfg.elec = elec;
        scd_data = ft_scalpcurrentdensity(ecfg, data);
        
        data = ft_appenddata([],nonscalp_data,scd_data);
        clear scd_data ecfg nonscalp_data;
    end
    
    trialevents = unique(data.trialinfo,'stable');
    
    for j = 1:length(trialevents)
        
        
        if Deci.Analysis.ArtifactReject
            
            if exist([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} filesep num2str(j) '.mat']) == 2
                artifacts = [];
                load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} filesep num2str(j) '.mat'],'artifacts');
            else
                error(['artifacts not found for ' Deci.SubjectList{subject_list}]);
            end
        else
            artifacts = logical(ones([1 length(find(data.trialinfo==trialevents(j)))]))';
        end
        
        cfg = [];
        cfg.trials = find(data.trialinfo==trialevents(j));
        cfg.trials = cfg.trials(artifacts);
        
        if Deci.Analysis.ERP
            
            
            time =  ft_timelockanalysis(cfg, data);
            
            mkdir([Deci.Folder.Analysis filesep 'Volt_ERP' filesep Deci.SubjectList{subject_list}])
            label = rmfield(time,'avg');
            save([Deci.Folder.Analysis filesep 'Volt_ERP' filesep Deci.SubjectList{subject_list} filesep num2str(j)],'time','label');
            
        end
        
        if ~isempty(Deci.Analysis.Freq)
            
            if ~isfield(Deci.Analysis.Freq,'Toi')
                Deci.Analysis.Toi = [-inf inf];
                Warning('Parameter for Toi not found, presuming [-inf inf]')
            end
            
            fcfg = Deci.Analysis.Freq;
            fcfg.toi = Deci.Analysis.Toi(1):round(diff([data.time{1}(1) data.time{1}(2)]),5):Deci.Analysis.Toi(2);
            fcfg.output='fourier';
            fcfg.pad = 'maxperlen';
            fcfg.scc = 0;
            fcfg.keeptapers = 'no';
            fcfg.keeptrials = 'yes';
            fcfg.trials = cfg.trials;
            
            Fourier = rmfield(ft_freqanalysis(fcfg, data),'cfg');
            
            if ischar(Deci.Analysis.Channels)
                Chan = Fourier.label;
            elseif all(ismember(Deci.Analysis.Channels,Fourier.label))
                Chan = Deci.Analysis.Channels;
            else
                error('Wrong Channel Selection in Analyis');
            end
            
            for i = 1:length(Chan)
                
                dcfg = [];
                dcfg.channel = Chan(i);
                freqplaceholder = ft_selectdata(dcfg,Fourier);
                
                
                mkdir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list} filesep num2str(j)]);
                mkdir([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep Deci.SubjectList{subject_list} filesep num2str(j)]);
                
                label = freqplaceholder;
                label = rmfield(label,'fourierspctrm');
                label.label = Chan;
                label.dimord = 'chan_freq_time';
                
                freq = freqplaceholder;
                freq.dimord = 'chan_freq_time';
                freq.powspctrm      = permute(abs(mean(freq.fourierspctrm./abs(freq.fourierspctrm),1)),[2 3 4 1]);         % divide by amplitude
                freq  = rmfield(freq,'fourierspctrm');
                save([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep Deci.SubjectList{subject_list} filesep num2str(j) filesep Chan{i}],'freq','label','-v7.3');
                
                freq = freqplaceholder;
                freq.powspctrm = permute(mean(abs(freq.fourierspctrm).^2 ,1),[2 3 4 1]);
                freq.dimord = 'chan_freq_time';
                freq  = rmfield(freq,'fourierspctrm');
                save([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list} filesep num2str(j) filesep Chan{i}],'freq','label','-v7.3');
                
            end
            
        end
    end
    
    clear data
    
end

end


