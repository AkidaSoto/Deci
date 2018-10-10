function AllFreq = PCAnalyzor(Deci,subject_list)

AllFreq = [];

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

if ~isempty(Deci.Analysis.Freq.CFC.Within)
    
    if isfield(Deci.Analysis.Freq.CFC.Within,'freqhigh')
        if ischar(Deci.Analysis.Freq.CFC.Within.freqhigh)
            switch Deci.Analysis.Freq.CFC.Within.freqhigh
                case 'theta'
                    Deci.Analysis.Freq.CFC.Within.freqhigh = [4 8];
                case 'beta'
                    Deci.Analysis.Freq.CFC.Within.freqhigh = [12.5 30];
                case 'alpha'
                    Deci.Analysis.Freq.CFC.Within.freqhigh = [8 12.5];
                case 'gamma'
                     Deci.Analysis.Freq.CFC.Within.freqhigh = [30 50];
            end
        elseif isnumeric(Deci.Analysis.Freq.CFC.Within.freqhigh)
        else
            error(['cannot interrept freqhigh']);
        end
    else
           error(['cannot interrept freqhigh']);
    end
    
    if isfield(Deci.Analysis.Freq.CFC.Within,'freqlow')
        if ischar(Deci.Analysis.Freq.CFC.Within.freqlow)
            switch Deci.Analysis.Freq.CFC.Within.freqlow
                case 'theta'
                    Deci.Analysis.Freq.CFC.Within.freqlow = [4 8];
                case 'beta'
                    Deci.Analysis.Freq.CFC.Within.freqlow = [12.5 30];
                case 'alpha'
                    Deci.Analysis.Freq.CFC.Within.freqlow = [8 12.5];
                case 'gamma'
                    Deci.Analysis.Freq.CFC.Within.freqlow = [30 50];
            end
        elseif isnumeric(Deci.Analysis.Freq.CFC.Within.freqlow)
        else
           error(['cannot interrept freqlow']);
        end
    else
        error(['cannot interrept freqlow']);
    end
end

if ~isempty(Deci.Analysis.Freq.CFC.Between)
    
    if isfield(Deci.Analysis.Freq.CFC.Between,'freqhigh')
        if ischar(Deci.Analysis.Freq.CFC.Between.freqhigh)
            switch Deci.Analysis.Freq.CFC.Between.freqhigh
                case 'theta'
                    Deci.Analysis.Freq.CFC.Between.freqhigh = [4 8];
                case 'beta'
                    Deci.Analysis.Freq.CFC.Between.freqhigh = [12.5 30];
                case 'alpha'
                    Deci.Analysis.Freq.CFC.Between.freqhigh = [8 12.5];
                case 'gamma'
                     Deci.Analysis.Freq.CFC.Between.freqhigh = [30 50];
            end
        elseif isnumeric(Deci.Analysis.Freq.CFC.Between.freqhigh)
        else
            error(['cannot interrept freqhigh']);
        end
    else
           error(['cannot interrept freqhigh']);
    end
    
    if isfield(Deci.Analysis.Freq.CFC.Between,'freqlow')
        if ischar(Deci.Analysis.Freq.CFC.Between.freqlow)
            switch Deci.Analysis.Freq.CFC.Between.freqlow
                case 'theta'
                    Deci.Analysis.Freq.CFC.Between.freqlow = [4 8];
                case 'beta'
                    Deci.Analysis.Freq.CFC.Between.freqlow = [12.5 30];
                case 'alpha'
                    Deci.Analysis.Freq.CFC.Between.freqlow = [8 12.5];
                case 'gamma'
                    Deci.Analysis.Freq.CFC.Between.freqlow = [30 50];
            end
        elseif isnumeric(Deci.Analysis.Freq.CFC.Between.freqlow)
        else
           error(['cannot interrept freqlow']);
        end
    else
        error(['cannot interrept freqlow']);
    end
end

data = [];
load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list}],'data');

if Deci.Analysis.Laplace
    [elec.label, elec.elecpos] = elec_1020select(data.label);
    ecfg.elec = elec;
    data = ft_scalpcurrentdensity(ecfg, data);
end

trialevents = unique(data.trialinfo,'stable');

for Cond = 1:length(trialevents)

    cfg = [];
    cfg.trials = find(data.trialinfo==trialevents(Cond));
    
    redefine = 0;
    if exist([Deci.Folder.Version  filesep 'Redefine' filesep Deci.SubjectList{subject_list}  '.mat']) == 2
        redefine = 1;
        retrl = [];
        load([Deci.Folder.Version  filesep 'Redefine' filesep Deci.SubjectList{subject_list}  '.mat']);
    end
    
    
    if Deci.Analysis.ERP
        
        if redefine
            cfg.offset = retrl;
            cfg.shift =  Deci.Analysis.Redefine.ERPToi;
            datatime = ft_datashift(cfg,data);
        else
            datatime = data;
        end
        
        cfg.vartrllength = 2;
        time =  ft_timelockanalysis(cfg, datatime);
        clear datatime;
        
        mkdir([Deci.Folder.Analysis filesep 'Volt_ERP' filesep Deci.SubjectList{subject_list}])
        label = rmfield(time,'avg');
        save([Deci.Folder.Analysis filesep 'Volt_ERP' filesep Deci.SubjectList{subject_list} filesep num2str(Cond)],'time','label');
        
    end
    
    if ~isempty(Deci.Analysis.Freq)
        
        if ~isfield(Deci.Analysis.Freq,'Toi')
            Deci.Analysis.Toi = [-inf inf];
            Warning('Parameter for Toi not found, presuming [-inf inf]')
        end
        
        fcfg = Deci.Analysis.Freq;
        fcfg.toi = Deci.Analysis.Freq.Toi(1):round(diff([data.time{1}(1) data.time{1}(2)]),5):Deci.Analysis.Freq.Toi(2);
        fcfg.output='fourier';
        fcfg.pad = 'maxperlen';
        fcfg.scc = 0;
        fcfg.keeptapers = 'no';
        fcfg.keeptrials = 'yes';
        fcfg.trials = cfg.trials;
        
        Analysis = 1;
        
        if isfield(Deci.Analysis.Freq,'Redefine')
            
            retrl1 = retrl(find(data.trialinfo==trialevents(Cond)));
            
            fcfg.toi = [Deci.Analysis.Redefine.Bsl(1):diff([data.time{1}(1) data.time{1}(2)]):Deci.Analysis.Redefine.Bsl(2)];
            bsl = rmfield(ft_freqanalysis(fcfg, data),'cfg');
            bsl.powspctrm = permute(mean(mean(abs(bsl.fourierspctrm).^2 ,1),4),[2 3 1]);
            bsl.dimord = 'chan_freq_time';
            bsl = rmfield(bsl,'fourierspctrm');
            bsl.freq = bsl.oldfoi;
            
            Analysis = 0;
            if Deci.Analysis.Freq.Redefine  ~= 2
                begtim  = min(retrl1) + Deci.Analysis.Freq.Toi(1);
                endtim  = max(retrl1) + Deci.Analysis.Freq.Toi(2);
                fcfg.toi = [begtim:diff([data.time{1}(1) data.time{1}(2)]):endtim];
                
                Fourier = rmfield(ft_freqanalysis(fcfg, data),'cfg');
                
                shift_cfg.latency = Deci.Analysis.Freq.Toi(1):diff([data.time{1}(1) data.time{1}(2)]):Deci.Analysis.Freq.Toi(2);
                shift_cfg.offset = retrl1;
                shift_cfg.parameter = 'fourierspctrm';
                shift_cfg.keeptrials = 'no';
                
                Fourier = ft_freqshift(shift_cfg, Fourier);
                
                Analysis = 1;
            end
            
            mkdir([Deci.Folder.Version  filesep 'Redefine' filesep 'BSL' filesep Deci.SubjectList{subject_list}]);
            save([Deci.Folder.Version  filesep 'Redefine' filesep 'BSL' filesep Deci.SubjectList{subject_list} filesep num2str(Cond)],'bsl');
            clear bsl;
        else
            Fourier = rmfield(ft_freqanalysis(fcfg, data),'cfg');
        end
        
        
        if Analysis
            
            Fourier.freq = Fourier.oldfoi;
            
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
                
                mkdir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list} filesep num2str(Cond)]);
                mkdir([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep  Deci.SubjectList{subject_list} filesep num2str(Cond)]);
                
                label = freqplaceholder;
                label = rmfield(label,'fourierspctrm');
                label.label = Chan;
                label.dimord = 'chan_freq_time';
                
                freq = freqplaceholder;
                freq.dimord = 'chan_freq_time';
                freq.powspctrm      = permute(abs(mean(freq.fourierspctrm./abs(freq.fourierspctrm),1)),[2 3 4 1]);         % divide by amplitude
                freq  = rmfield(freq,'fourierspctrm');
                save([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep Deci.SubjectList{subject_list} filesep num2str(Cond) filesep Chan{i}],'freq','label','-v7.3');
                
                freq = freqplaceholder;
                freq.powspctrm = permute(mean(abs(freq.fourierspctrm).^2 ,1),[2 3 4 1]);
                freq.dimord = 'chan_freq_time';
                freq  = rmfield(freq,'fourierspctrm');
                save([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list} filesep num2str(Cond) filesep Chan{i}],'freq','label','-v7.3');
                
                
                
                
                 if ~isempty(Deci.Analysis.Freq.CFC.Within)
                     freq = freqplaceholder;
                     
                            
                           cscfg.latency =  Deci.Analysis.Freq.CFC.Within.latencylow;
                           cscfg.freqbin = Deci.Analysis.Freq.CFC.Within.freqlow;
                           cscfg.timebin = Deci.Analysis.Freq.CFC.Within.timebin;
                           freqlow = ft_binfreq(cscfg,freq);
                           
                           cscfg.latency =  Deci.Analysis.Freq.CFC.Within.latencyhigh;
                           cscfg.freqbin = Deci.Analysis.Freq.CFC.Within.freqhigh;
                           freqhigh = ft_binfreq(cscfg,freq);
                           
%                            cfccfg.channel = Deci.Analysis.Freq.CFC.Within.Channel;
%                            cfccfg.method = Deci.Analysis.Freq.CFC.Within.method;
%                            cfc =  ft_crossfrequencyanalysis(Deci.Analysis.Freq.CFC.Within,freqlow,freqhigh);


                            acfg.parameter = 'fourierspctrm';
                            acfg.appenddim = 'freq';
                            AllFreq{subject_list} = ft_appendfreq(acfg,freqlow,freqhigh);
                            
                 end
                 
                 if ~isempty(Deci.Analysis.Freq.CFC.Between)
                     
                     freq = freqplaceholder;

                     cscfg.latency =  Deci.Analysis.Freq.CFC.Between.latencylow;
                     cscfg.freqbin = Deci.Analysis.Freq.CFC.Between.Freqlow;
                     cscfg.timebin = Deci.Analysis.Freq.CFC.Between.timebin;
                     SubFreqLow{i} = ft_binfreq(cscfg,freq);
                     
                     cscfg.latency =  Deci.Analysis.Freq.CFC.Between.latencyhigh;
                     cscfg.freqbin = Deci.Analysis.Freq.CFC.Between.Freqhigh;
                     SubFreqHigh{i} = ft_binfreq(cscfg,freq);
                     
                 end
                 
            end
            
            if ~isempty(Deci.Analysis.Freq.CFC.Between)
                
                cacfg.parameter = 'fourierspctrm';
                cacfg.appenddim = 'chan';
                AllFreqLow{subject_list} = ft_appendfreq(cacfg,SubFreqLow{:});
                AllFreqHigh{subject_list} = ft_appendfreq(cacfg,SubFreqHigh{:});
                
                cfc=  ft_singlecfc(Deci.Analysis.Freq.CFC.Between,AllFreqLow{subject_list},AllFreqHigh{subject_list});
                
                mkdir([Deci.Folder.Analysis filesep 'Cross_' Deci.Analysis.Freq.CFC.Between.method filesep  Deci.SubjectList{subject_list}]);
                save([Deci.Folder.Analysis filesep 'Cross_' Deci.Analysis.Freq.CFC.Between.method filesep Deci.SubjectList{subject_list} filesep num2str(Cond)],'cfc','-v7.3');
                
            end
            
        end
        
    end
    
end


end



