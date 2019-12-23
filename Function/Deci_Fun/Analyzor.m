function Analyzor(Deci,subject_list)

data = [];

%% Load Data

% Load Definition Data if only using Extra.Once

if ~Deci.Analysis.Freq.do && ~Deci.Analysis.CFC.do && [Deci.Analysis.Extra.do && ~isempty(find(Deci.Analysis.Extra.Once))]
    noneeg_flag = 1;
    
else
    noneeg_flag = 0;
end

if noneeg_flag
    load([Deci.Folder.Definition filesep Deci.SubjectList{subject_list} '.mat'],'cfg');
    data =cfg;
    data.condinfo{1} = data.trl(:,end-length(Deci.DT.Locks)+1:end);
    data.condinfo{2} = data.event;
    data.condinfo{3} = data.trialnum;
    
    data.preart = data.condinfo;
    
else
    load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data');
end

condinfo = data.condinfo;
preart = data.preart;

%preart will be used to anaylzed data and then those trials will be
%selected from condinfo to only include not-artifact trials


%% Laplace Transformation
if Deci.Analysis.Laplace
    [elec.label, elec.elecpos] = CapTrakMake([Deci.Folder.Raw  filesep Deci.SubjectList{subject_list} '.bvct']);
    ecfg.elec = elec;
    data = ft_scalpcurrentdensity(ecfg, data);
end

%% HemifieldFlip

if Deci.Analysis.HemifieldFlip.do
    
    Hemifields = preart{2}(:,find(mean(ismember(preart{2},Deci.Analysis.HemifieldFlip.Markers),1)));
    
    FlipCfg.trials  = ismember(Hemifields,Deci.Analysis.HemifieldFlip.Markers(2));
    
    FlipData = ft_selectdata(FlipCfg,data);
    FlipData = hemifieldflip(FlipData);
    
    FlipCfg.trials = ~FlipCfg.trials;
    
    NotFlipData = ft_selectdata(FlipCfg,data);
    
    data = ft_appenddata([],FlipData,NotFlipData);
end

if ~strcmpi(Deci.Analysis.Channels,'all')
    cfg = [];
    cfg.channel = Deci.Analysis.Channels;
    data = ft_selectdata(cfg,data);
end
%% Downsample
if ~isempty(Deci.Analysis.DownSample)
    data = ft_resampledata(struct('resamplefs',Deci.Analysis.DownSample,'detrend','no'),data);
end


data.condinfo = condinfo;
data.preart = preart;

if Deci.Analysis.HemifieldFlip.do
    data.condinfo = cellfun(@(a,b) [a;b],FlipData.condinfo,NotFlipData.condinfo,'UniformOutput',false);
    data.preart = cellfun(@(a,b) [a;b],FlipData.preart,NotFlipData.preart,'UniformOutput',false);
end

%% Loop through Conditions
for Cond = 1:length(Deci.Analysis.Conditions)
    
    TimerCond = clock;
    
    %% do Extra.Once Analyses
    if Deci.Analysis.Extra.do
        
        info.subject_list = subject_list;
        info.Cond = Cond;
        
        if isfield(Deci.Analysis.Extra,'Once')
            
            for funs = find(Deci.Analysis.Extra.Once)
                
                if Deci.Analysis.Extra.list(funs)
                    feval(Deci.Analysis.Extra.Functions{funs},Deci,info,data,Deci.Analysis.Extra.Params{funs}{:});
                end
            end
        end
        
    end
    
    
    %% Find Relevant Trials from that Condition info
    maxt = length(find(cellfun(@(c) any(ismember(Deci.Analysis.Conditions{Cond},c)), Deci.DT.Markers)));
    info.alltrials = find(sum(ismember(preart{2},Deci.Analysis.Conditions{Cond}),2) == maxt);
    
    %ignore all locks that are missing
    minamountofnans = min(mean(isnan(preart{1}(info.alltrials,:)),2));
    info.allnonnans = mean(isnan(preart{1}(info.alltrials,:)),2) == minamountofnans;% & ~isnan(mean(condinfo{2},2));
    
    ccfg.trials =  info.alltrials(info.allnonnans);
    
    if Deci.Analysis.ApplyArtReject
        ccfg.trials = ccfg.trials(find(ismember(preart{3}(info.alltrials(info.allnonnans)),condinfo{3})));
    end
    
    dataplaceholder = ft_selectdata(ccfg,data);
    dataplaceholder.condinfo = dataplaceholder.preart; %preart no longer needed
    dataplaceholder = rmfield(dataplaceholder,'preart');
    
    %% Loop Through Locks
    
    if Deci.Analysis.Freq.do || Deci.Analysis.CFC.do || [Deci.Analysis.Extra.do && ~isempty(find(~Deci.Analysis.Extra.Once))] || Deci.Analysis.ERP.do
        
        for Lock = 1:length(Deci.Analysis.Locks)
            
            cfg.offset = preart{1}(ccfg.trials,Deci.Analysis.Locks(Lock));
            cfg.toilim = Deci.Analysis.Freq.Toilim;
            dat = ft_datashift2(cfg,dataplaceholder);
            
            %% Do ERP Analysis
            if Deci.Analysis.ERP.do
                mkdir([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list} filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond}]);
                ecfg.latency = Deci.Analysis.ERP.Toi;
                
                for chan = 1:length(dat.label)
                    ecfg.channel = dat.label(chan);
                    erp = ft_selectdata(ecfg,dat);
                    save([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list} filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond} filesep dat.label{chan}],'erp');
                end
                clear erp
            end
            
            
            %% Do Freq Analyses
            if Deci.Analysis.Freq.do || Deci.Analysis.CFC.do || [Deci.Analysis.Extra.do && ~isempty(find(~Deci.Analysis.Extra.Once))]
                
                skip = false;
                if Deci.Analysis.Freq.Skip_if_done && ~Deci.Analysis.Extra.do
                    if exist([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}  filesep filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond}])...
                            && exist([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep Deci.SubjectList{subject_list}  filesep filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond}])
                        
                        skip = true;
                    end
                end
                
                
                if ~skip
                    
                    if ~strcmp(Deci.Analysis.Freq.method,'hilbert')
                        fcfg = Deci.Analysis.Freq;
                        fcfg.output='fourier';
                        fcfg.pad = 'maxperlen';
                        fcfg.scc = 0;
                        fcfg.keeptapers = 'yes';
                        fcfg.keeptrials = 'yes';
                        fcfg.toi = Deci.Analysis.Freq.Toi(1):round(diff([data.time{1}(1) data.time{1}(2)]),5):Deci.Analysis.Freq.Toi(2);
                        
                        Fourier = rmfield(ft_freqanalysis(fcfg, dat),'cfg');
                        Fourier.condinfo = dat.condinfo;
                        trllength = size(Fourier.fourierspctrm,1);
                    else
                        
                        fcfg = Deci.Analysis.Freq;
                        nyquist = data.fsample/2;
                        
                        freqs = Deci.Analysis.Freq.foi;
                        
                        tempfreq = [];
                        
                        for foi = 1:length(freqs)
                            
                            hcfg = [];
                            hcfg.bpfilter2 = 'yes';  %Modified implementation to work with MikexCohen's formula
                            hcfg.bpfreq =[freqs(foi)-fcfg.width(foi) freqs(foi)+fcfg.width(foi)];
                            hcfg.bpfiltord = round(fcfg.order*(data.fsample/hcfg.bpfreq(1)));
                            hcfg.bpfilttype = 'firls';
                            hcfg.transition_width = fcfg.transition_width;
                            hcfg.hilbert = 'complex';
                            
                            evalc('hil = ft_preprocessing(hcfg,dat)');
                            
                            rcfg.latency = [Deci.Analysis.Freq.Toi];
                            Fo = ft_selectdata(rcfg,hil);
                            
                            tempfreq{foi}.fourierspctrm = permute(cell2mat(permute(Fo.trial,[3 1 2])),[3 1 4 2]);
                            tempfreq{foi}.label = Fo.label;
                            tempfreq{foi}.freq = freqs(foi);
                            tempfreq{foi}.trialinfo = Fo.trialinfo;
                            tempfreq{foi}.time = Fo.time{1}';
                            tempfreq{foi}.dimord = 'rpt_chan_freq_time';
                            
                        end
                        
                        acfg.parameter = 'fourierspctrm';
                        acfg.appenddim = 'freq';
                        
                        Fourier = rmfield(ft_appendfreq(acfg,tempfreq{:}),'cfg');
                        Fourier.dimord = 'rpt_chan_freq_time';
                        Fourier.condinfo = dat.condinfo;
                        trllength = size(Fourier.fourierspctrm,1);
                    end
                    
                    Chan = Fourier.label;
                    TimerChan = clock;
                    
                    %% Loop Through Channels
                    for i = 1:length(Chan)
                        dcfg = [];
                        dcfg.channel = Chan(i);
                        freq = ft_selectdata(dcfg,Fourier);
                        
                        warning('off', 'MATLAB:MKDIR:DirectoryExists');
                        mkdir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}  filesep filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond}]);
                        mkdir([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep  Deci.SubjectList{subject_list}  filesep filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond}]);
                        mkdir([Deci.Folder.Analysis filesep 'Freq_TotalPowerVar' filesep Deci.SubjectList{subject_list}  filesep filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond}]);
                        
                        freqplaceholder = freq;
                        
                        if Deci.Analysis.Freq.do
                            freq = freqplaceholder;
                            freq.dimord = 'chan_freq_time';
                            freq.powspctrm      = permute(abs(mean(freq.fourierspctrm./abs(freq.fourierspctrm),1)),[2 3 4 1]);         % divide by amplitude
                            freq  = rmfield(freq,'fourierspctrm');
                            freq.trllength = trllength;
                            save([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep Deci.SubjectList{subject_list}  filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond} filesep Chan{i}],'freq','-v7.3');
                            
                            
                            freq = freqplaceholder;
                            freq.powspctrm = permute(mean(abs(freq.fourierspctrm).^2 ,1),[2 3 4 1]);
                            freq.dimord = 'chan_freq_time';
                            freq  = rmfield(freq,'fourierspctrm');
                            freq.trllength = trllength;
                            save([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}  filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond} filesep Chan{i}],'freq','-v7.3');
                        end
                        
                        if isfield(Deci.Analysis.Extra,'do')
                            if Deci.Analysis.Extra.do
                                info.Channels = Chan;
                                info.ChanNum = i;
                                info.Lock = Lock;
                                if isfield(Deci.Analysis.Extra,'Once')
                                    for funs = find(~Deci.Analysis.Extra.Once)
                                        if Deci.Analysis.Extra.list(funs)
                                            feval(Deci.Analysis.Extra.Functions{funs},Deci,info,freqplaceholder,Deci.Analysis.Extra.Params{funs}{:});
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Do Connectivity Analysis
                
                if Deci.Analysis.Connectivity.do
                    info.Cond = Cond;
                    info.Lock = Lock;
                    info.subject_list = subject_list;
                    for funs = find(Deci.Analysis.Connectivity.list)
                        feval(Deci.Analysis.Connectivity.Functions{funs},Deci,info,dat,data,Deci.Analysis.Connectivity.Params{funs}{:});
                    end
                    
                end
                
                if ~skip
                    disp(['s:' num2str(subject_list) ' c:' num2str(Cond) ' Lock' num2str(Lock) ' time: ' num2str(etime(clock ,TimerChan))]);
                end
            end
        end
    end
    
    disp(['s:' num2str(subject_list) ' c:' num2str(Cond) ' Cond time: ' num2str(etime(clock ,TimerCond))]);
end
end


