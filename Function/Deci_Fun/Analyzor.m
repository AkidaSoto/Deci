function Analyzor(Deci,subject_list)
display('----------------------');
display(['Starting Analyzor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
display(' ')

%% Load Data
data = [];
info.subject_list = subject_list;

% Load Definition Data if only using Extra.Once
if Deci.Analysis.Behavioral
    load([Deci.Folder.Definition filesep Deci.SubjectList{subject_list} '.mat'],'cfg');
    data =cfg;
    data.condinfo{1} = data.trl(:,end-length(Deci.DT.Locks)+1:end);
    data.condinfo{2} = data.event;
    data.condinfo{3} = data.trialnum;
    data.preart = data.condinfo;
    display('Using Behavioral Data, Not EEG for current analysis')
else
    load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data');
end

if isfield(data,'condinfo')  %replacer starting 12/22, lets keep for ~4 months
    data.postart.locks = data.condinfo{1};
    data.postart.events = data.condinfo{2};
    data.postart.trlnum = data.condinfo{3};
    
    data.locks = data.preart{1};
    data.events = data.preart{2};
    data.trlnum = data.preart{3};
    
    data = rmfield(data,'condinfo');
    data = rmfield(data,'preart');
end

locks = data.locks;
events = data.events;
trlnum = data.trlnum;
postart.locks = data.postart.locks;
postart.events = data.postart.events;
postart.trlnum = data.postart.trlnum;

%% Laplace Transformation
if Deci.Analysis.Laplace
    [elec.label, elec.elecpos] = CapTrakMake([Deci.Folder.Raw  filesep Deci.SubjectList{subject_list} '.bvct']);
    ecfg.elec = elec;
    evalc('data = ft_scalpcurrentdensity(ecfg, data)');
    
    display('Laplace Transform Applied')
end

if ~strcmpi(Deci.Analysis.Channels,'all')
    cfg = [];
    cfg.channel = Deci.Analysis.Channels;
    evalc('data = ft_selectdata(cfg,data)');
    display('Channel Selection Applied')
end
%% Downsample
if ~isempty(Deci.Analysis.DownSample)
    rcfg = struct('resamplefs',Deci.Analysis.DownSample,'detrend','no');
    evalc('data = ft_resampledata(rcfg,data)');
    display('Downsampling Applied')
end

data.locks = locks;
data.events = events;
data.trlnum = trlnum;
data.postart.locks = postart.locks;
data.postart.events = postart.events;
data.postart.trlnum = postart.trlnum;
info.postart = data.postart;

%% HemifieldFlip
%check to see if postart is pull through selectdata
if Deci.Analysis.HemifieldFlip.do
    Hemifields = data.events(:,find(mean(ismember(data.events,Deci.Analysis.HemifieldFlip.Markers),1)));
    
    FlipCfg.trials  = ismember(Hemifields,Deci.Analysis.HemifieldFlip.Markers(2));
    evalc('FlipData = ft_selectdata(FlipCfg,data)');
    evalc('FlipData = hemifieldflip(FlipData)');
    
    nonFlipCfg.trials = ismember(Hemifields,Deci.Analysis.HemifieldFlip.Markers(1));
    evalc('NotFlipData = ft_selectdata(nonFlipCfg,data)');
    
    evalc('data = ft_appenddata([],FlipData,NotFlipData)');
    
    data.locks = locks([find(FlipCfg.trials);find(nonFlipCfg.trials)],:);
    data.events = events([find(FlipCfg.trials);find(nonFlipCfg.trials)],:);
    data.trlnum = trlnum([find(FlipCfg.trials);find(nonFlipCfg.trials)]);
    data.postart.locks = postart.locks([find(ismember(postart.trlnum,find(FlipCfg.trials)));find(ismember(postart.trlnum,find(nonFlipCfg.trials)))],:);
    data.postart.events = postart.events([find(ismember(postart.trlnum,find(FlipCfg.trials)));find(ismember(postart.trlnum,find(nonFlipCfg.trials)))],:);
    data.postart.trlnum = postart.trlnum([find(ismember(postart.trlnum,find(FlipCfg.trials)));find(ismember(postart.trlnum,find(nonFlipCfg.trials)))]);

    locks = data.locks;
    events = data.events;
    trlnum = data.trlnum;
    postart.locks = data.postart.locks;
    postart.events = data.postart.events;
    postart.trlnum = data.postart.trlnum;
    
    display('HemifieldFlip Applied')
end

if Deci.Analysis.Unique.do
    for funs = find(Deci.Analysis.Unique.list)
        if Deci.Analysis.Unique.list(funs)
            disp(['running Unique: ' Deci.Analysis.Unique.Functions{funs}])
            feval(Deci.Analysis.Unique.Functions{funs},Deci,info,data,Deci.Analysis.Unique.Params{funs}{:});
        end
    end
end

%% Loop through Conditions
for Cond = 1:length(Deci.Analysis.Conditions)
    info.Cond = Cond;
    display(' ')
    display(['---Starting Condition #' num2str(Cond) ': ' Deci.Analysis.CondTitle{Cond} '---'])
    display(' ')
    
    %% Find Relevant Trials from that Condition info
    
    if ~all(ismember(Deci.Analysis.Conditions{Cond},events))
        
        if ~strcmpi(Deci.DT.Type,'EEG_Polymerase')
            display(['Using unique Trial Def:' Deci.DT.Type])
        end
        error('1 or More Marker Codes not found in events. If using own DT, make sure Step 4 contains updated DT.Markers field')
    end
        
    maxt = length(find(cellfun(@(c) any(ismember(Deci.Analysis.Conditions{Cond},c)), Deci.DT.Markers)));
    info.alltrials = find(sum(ismember(events,Deci.Analysis.Conditions{Cond}),2) == maxt);
    %% ignore all locks with missing nans
    if Deci.Analysis.IgnoreNanLocks
        minamountofnans = min(mean(isnan(locks(info.alltrials,:)),2));
        info.nanlocks = mean(isnan(locks(info.alltrials,:)),2) ~= minamountofnans;
        
        if any(info.nanlocks)
            display(['ignoring ' num2str(length(find(info.nanlocks))) ' trials with missing locks'])
        end
    else
        info.nanlocks = logical(size(info.alltrials));
    end
    
        %% Reject Arts
    ccfg.trials =  info.alltrials(~info.nanlocks);
    if Deci.Analysis.ApplyArtReject
        ccfg.trials = ccfg.trials(ismember(trlnum(ccfg.trials),postart.trlnum));
        display('Applying Artifact Rejection')
        display(['Final trial count is ' num2str(length(ccfg.trials))])
    else
        display('Not Applying Artifact Rejection')
    end
    dataplaceholder = ft_selectdata(ccfg,data);
    
    if length(dataplaceholder.trial) < 40
       display('!!!Trial Count for this condition is < 40, which is abnormally low!!!') 
    end
    
    %% do Extra Analyses
    if Deci.Analysis.Extra.do
        for funs = find(Deci.Analysis.Extra.list)
            if Deci.Analysis.Extra.list(funs)
                disp(['running Extra: ' Deci.Analysis.Extra.Functions{funs}])
                feval(Deci.Analysis.Extra.Functions{funs},Deci,info,dataplaceholder,Deci.Analysis.Extra.Params{funs}{:});
            end
        end
    end
    %% Loop Through Locks
    
    if Deci.Analysis.Freq.do || Deci.Analysis.CFC.do || Deci.Analysis.ERP.do
        
        for Lock = 1:length(Deci.Analysis.Locks)
            display(' ')
            display(['---Starting Lock #' num2str(Lock) ': ' Deci.Analysis.LocksTitle{Lock} '---'])
            display(' ')
            info.Lock = Lock;
            
            cfg.offset = locks(ccfg.trials,Deci.Analysis.Locks(Lock));
            cfg.toilim = Deci.Analysis.Freq.Toilim;
            evalc('dat = ft_datashift2(cfg,dataplaceholder)');
            
            for lockstd = 1:size(dat.trialinfo,2)
              lockers(lockstd)  =  mean(dat.trialinfo(:,Lock) - dat.trialinfo(:,lockstd));
            end
            info.lockers = lockers;

            %% Do ERP Analysis
            if Deci.Analysis.ERP.do
                mkdir([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list} filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond}]);
                ecfg.latency = Deci.Analysis.ERP.Toi;
                
                for chan = 1:length(dat.label)
                    ecfg.channel = dat.label(chan);
                    erp = ft_selectdata(ecfg,dat);
                    erp.locker = lockers;
                    save([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list} filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond} filesep dat.label{chan}],'erp');
                end
                clear erp
            end
            
            %% Do Freq Analyses
                if ~strcmp(Deci.Analysis.Freq.method,'hilbert')
                    fcfg = Deci.Analysis.Freq;
                    fcfg.output='fourier';
                    fcfg.pad = 'maxperlen';
                    fcfg.scc = 0;
                    fcfg.keeptapers = 'yes';
                    fcfg.keeptrials = 'yes';
                    fcfg.toi = Deci.Analysis.Freq.Toi(1):round(diff([data.time{1}(1) data.time{1}(2)]),5):Deci.Analysis.Freq.Toi(2);
                    
                    Fourier = rmfield(ft_freqanalysis(fcfg, dat),'cfg');
                    trllength = size(Fourier.fourierspctrm,1);
                else
                    display('Applying Hilbert Transformation')
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
                    trllength = size(Fourier.fourierspctrm,1);
                end
                
                Chan = Fourier.label;
                
                %% Loop Through Channels
                for i = 1:length(Chan)
                    info.Channels = Chan;
                    info.ChanNum = i;
                    
                    dcfg = [];
                    dcfg.channel = Chan(i);
                    freq = ft_selectdata(dcfg,Fourier);
                    
                    warning('off', 'MATLAB:MKDIR:DirectoryExists');
                    mkdir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}  filesep filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond}]);
                    mkdir([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep  Deci.SubjectList{subject_list}  filesep filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond}]);

                    freqplaceholder = freq;
                    if Deci.Analysis.Freq.do
                        freq = freqplaceholder;
                        freq.dimord = 'chan_freq_time';
                        freq.powspctrm      = permute(abs(mean(freq.fourierspctrm./abs(freq.fourierspctrm),1)),[2 3 4 1]);         % divide by amplitude
                        freq  = rmfield(freq,'fourierspctrm');
                        freq.trllength = trllength;
                        freq.lockers = lockers;
                        
                        save([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep Deci.SubjectList{subject_list}  filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond} filesep Chan{i}],'freq','-v7.3');
                        
                        
                        freq = freqplaceholder;
                        freq.powspctrm = permute(mean(abs(freq.fourierspctrm).^2 ,1),[2 3 4 1]);
                        freq.dimord = 'chan_freq_time';
                        freq  = rmfield(freq,'fourierspctrm');
                        freq.trllength = trllength;
                        freq.lockers = lockers;
                        save([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}  filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond} filesep Chan{i}],'freq','-v7.3');
                    end
                    
                    if Deci.Analysis.Freq.Extra.do
                        for funs = find(Deci.Analysis.Freq.Extra.list)
                            if Deci.Analysis.Freq.Extra.list(funs)
                                feval(Deci.Analysis.Freq.Extra.Functions{funs},Deci,info,freqplaceholder,Deci.Analysis.Freq.Extra.Params{funs}{:});
                            end
                        end
                    end
                end
                
                % Do Connectivity Analysis
                if Deci.Analysis.Connectivity.do
                    for funs = find(Deci.Analysis.Connectivity.list)
                        feval(Deci.Analysis.Connectivity.Functions{funs},Deci,info,Fourier,Deci.Analysis.Connectivity.Params{funs}{:});
                    end
                end
        end
    end
end
end


