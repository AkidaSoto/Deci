
function dc_source(Deci,info,data,params)

%http://www.fieldtriptoolbox.org/workshop/oslo2019/forward_modeling/
%% Load Standards


%load('standard_sourcemodel3d8mm.mat','sourcemodel')

%Load Elec

elec = ft_read_sens('standard_1020.elc');

if exist([Deci.Folder.Raw  filesep Deci.SubjectList{info.subject_list} '.bvct']) == 2
    [elec.label, elec.elecpos] = CapTrakMake([Deci.Folder.Raw  filesep Deci.SubjectList{info.subject_list} '.bvct']);
end
eleccheck = find(ismember(elec.label,[data.label; {'Nasion'};{'LPA'};{'RPA'}]));

elec.elecpos = elec.elecpos(eleccheck,:);
elec.chanpos = elec.elecpos;
elec.chantype = elec.chantype(eleccheck);
elec.chanunit = elec.chanunit(eleccheck);
elec.label = elec.label(eleccheck);
elec = ft_convert_units(elec, 'cm');

%% Electro realignment

if exist('C:\Users\User\Documents\GitHub\Deci\Function\Xtra_Fun\Subject01_seg.mat') == 0
    mri = ft_read_mri('standard_seg.mat');
    mri = ft_convert_units(mri, 'cm');
    
    cfg = [];
    mrirs = ft_volumereslice(cfg,mri);
    
    cfg=[];
    cfg.output    = {'brain','skull','scalp'};
    mri_seg =ft_volumesegment(cfg,mrirs);
    
    save('C:\Users\User\Documents\GitHub\Deci\Function\Xtra_Fun\Subject01_seg.mat','mri_seg')
else
    load('C:\Users\User\Documents\GitHub\Deci\Function\Xtra_Fun\Subject01_seg')
end

if exist('C:\Users\User\Documents\GitHub\Deci\Function\Xtra_Fun\Subject01_bnd.mat') == 0
    cfg=[];
    cfg.tissues= {'scalp' 'skull' 'brain' };
    cfg.numvertices = [6000 4000 2000];
    bnd=ft_prepare_mesh(cfg,mri_seg);
    
    save('C:\Users\User\Documents\GitHub\Deci\Function\Xtra_Fun\Subject01_bnd.mat','bnd')
else
    load('C:\Users\User\Documents\GitHub\Deci\Function\Xtra_Fun\Subject01_bnd.mat')
    
end

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Source' filesep 'elec_realigned' ]);

if exist([Deci.Folder.Analysis filesep 'Extra' filesep 'Source' filesep 'elec_realigned' filesep Deci.SubjectList{info.subject_list} '.mat']) == 0
    
    cfg = [];
    cfg.method = 'fiducial';
    cfg.target = elec;
    cfg.target.pos(1,:) = elec.chanpos(ismember('Nasion',elec.label),:);     % location of the nose
    cfg.target.pos(2,:) = elec.chanpos(ismember('LPA',elec.label),:);     % location of the left ear
    cfg.target.pos(3,:) = elec.chanpos(ismember('RPA',elec.label),:);     % location of the right ear
    cfg.target.type =  'eeg1010';
    cfg.elec = elec;
    elec_realigned = ft_electroderealign(cfg);
    
    cfg = [];
    cfg.method = 'headshape';
    cfg.headshape = bnd(1); % scalp surface
    cfg.elec = elec;
    elec_realigned = ft_electroderealign(cfg);
    
    
    cfg = [];
    cfg.method = 'project';
    cfg.headshape = bnd(1); % scalp surface
    cfg.elec = elec_realigned;
    elec_realigned = ft_electroderealign(cfg);
    
    a = figure;
    hold on
    ft_plot_sens(elec_realigned, 'elecsize', 40);
    ft_plot_headshape(bnd, 'facealpha', 0.1);
    view(90, 0)
    title(Deci.SubjectList{info.subject_list})
    
    savefig(a, [Deci.Folder.Analysis filesep 'Extra' filesep 'Source' filesep 'elec_realigned' filesep Deci.SubjectList{info.subject_list}]);
    
    % cfg = [];
    % cfg.method = 'interactive';
    % cfg.headshape = bnd(1); % scalp surface
    % cfg.elec = elec_realigned;
    % ft_electroderealign(cfg);
    
    save([Deci.Folder.Analysis filesep 'Extra' filesep 'Source' filesep 'elec_realigned' filesep Deci.SubjectList{info.subject_list} '.mat'],'elec_realigned')
else
    load([Deci.Folder.Analysis filesep 'Extra' filesep 'Source' filesep 'elec_realigned' filesep Deci.SubjectList{info.subject_list} '.mat'])
end

%% headmodel
if exist('C:\Users\User\Documents\GitHub\Deci\Function\Xtra_Fun\Subject01_vol.mat') == 0
    
    cfg = [];
    cfg.method='openmeeg';
    vol = ft_prepare_headmodel(cfg,bnd);
    save('C:\Users\User\Documents\GitHub\Deci\Function\Xtra_Fun\Subject01_vol','vol', '-v7.3')
else
    load('C:\Users\User\Documents\GitHub\Deci\Function\Xtra_Fun\Subject01_vol.mat')
end
%% create grid source\

if exist('C:\Users\User\Documents\GitHub\Deci\Function\Xtra_Fun\Subject01_sm.mat') == 0
    % Load Grid
    grid = ft_read_headshape('cortex_5124.surf.gii');
    grid = ft_convert_units(grid, 'cm');
    
    cfg             = [];
    cfg.headmodel   = vol; % used to estimate extent of grid
    cfg.resolution  = .75; % a source per 0.01 m -> 1 cm
    cfg.elec = elec_realigned;
    cfg.inwardshift = 0.005; % moving sources 5 mm inwards from the skull, ...
    % since BEM models may be unstable her
    sourcemodel = ft_prepare_sourcemodel(cfg);
    
    % figure
    % ft_plot_mesh(sourcemodel, 'vertexsize', 20);
    % ft_plot_vol(vol, 'facealpha', 0.5)
    % view(90, 0)
    save('C:\Users\User\Documents\GitHub\Deci\Function\Xtra_Fun\Subject01_sm', 'sourcemodel')
else
    
    load('C:\Users\User\Documents\GitHub\Deci\Function\Xtra_Fun\Subject01_sm.mat')
end

%% lf
%
% figure
ft_plot_sens(elec_realigned)
ft_plot_vol(vol,'facealpha',.05)
 ft_plot_mesh(sourcemodel)

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Source' filesep 'lf' ]);
if exist([Deci.Folder.Analysis filesep 'Extra' filesep 'Source' filesep 'lf' filesep Deci.SubjectList{info.subject_list} '.mat']) == 0
    
    lcfg                 = [];
    lcfg.elec            = elec_realigned;
    lcfg.channel          =  data.label;
    lcfg.grid = sourcemodel;
    lcfg.headmodel    = vol;
    lcfg.senstype = 'EEG';
    lcfg.normalize = 'yes';
    %Fourier.label = cellfun(@lower,Fourier.label,'UniformOutput',false);
    
    lf = ft_prepare_leadfield(lcfg);
    lf.inside(cellfun(@rank,lf.leadfield) ~= 3) = 0;
    
    
    save([Deci.Folder.Analysis filesep 'Extra' filesep 'Source' filesep 'lf' filesep Deci.SubjectList{info.subject_list} ],'lf')
else
    load([Deci.Folder.Analysis filesep 'Extra' filesep 'Source' filesep 'lf' filesep Deci.SubjectList{info.subject_list} '.mat'])
end

for chan1 = 1:length(lf.label)
    for chan2 = 1:length(data.label)
        
        if strcmpi(lf.label(chan1),data.label(chan2))
            lf.label(chan1) = data.label(chan2);
            elec_realigned.label(chan1) = data.label(chan2);
        end
    end
end

%% Freq Analysis on All Conditions first
for Lock = 1:length(Deci.Analysis.Locks)
    
    %% ignore all locks with missing nans
    if Deci.Analysis.IgnoreNanLocks
        minamountofnans = min(mean(isnan(data.locks),2));
        info.nanlocks = mean(isnan(data.locks),2) ~= minamountofnans;
        
        if any(info.nanlocks)
            display(['ignoring ' num2str(length(find(info.nanlocks))) ' trials with missing locks'])
        end
    else
        info.nanlocks = logical(size(info.alltrials));
    end
    
    ccfg.trials =  data.trlnum(~info.nanlocks);
    if Deci.Analysis.ApplyArtReject
        ccfg.trials = data.trlnum(ismember(data.trlnum,data.postart.trlnum));
        display('Applying Artifact Rejection')
        display(['Final trial count is ' num2str(length(ccfg.trials))])
    else
        display('Not Applying Artifact Rejection')
    end
    
    
    dataplaceholder = ft_selectdata(ccfg,data);
    
    display(' ')
    display(['---Starting Lock #' num2str(Lock) ': ' Deci.Analysis.LocksTitle{Lock} '---'])
    display(' ')
    info.Lock = Lock;
    
    cfg.offset = data.locks(ccfg.trials,Deci.Analysis.Locks(Lock));
    
    cfg.toilim = Deci.Analysis.Toilim;
    evalc('dat = ft_datashift2(cfg,dataplaceholder)');
    
    for lockstd = 1:size(dat.trialinfo,2)
        lockers(lockstd)  =  mean(dat.trialinfo(:,Lock) - dat.trialinfo(:,lockstd));
    end
    info.lockers = lockers;
    
    %% Do ERP Analysis
    if Deci.Analysis.ERP.do
        mkdir([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list} filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond}]);
        ecfg.latency = Deci.Analysis.Toi;
        
        for chan = 1:length(dat.label)
            ecfg.channel = dat.label(chan);
            erp = ft_selectdata(ecfg,dat);
            
            if isfield(Deci.Analysis.ERP, 'filter')
                
            end
            
            evalc('erp = ft_timelockanalysis([],erp)');
            erp.lockers = lockers;
            erp.trllength = size(dat.trialinfo,1);
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
        fcfg.toi = Deci.Analysis.Toi(1):round(diff([data.time{1}(1) data.time{1}(2)]),5):Deci.Analysis.Toi(2);
        fcfg.gpu = Deci.GCom;
        fcfg.cpu = Deci.DCom;
        
        Fourier = dc_freqanalysis(fcfg, dat);
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
    
    for freq = 1:length(params.foi)
        
        LF = dc_findfreq(params.foi(freq));
        fcfg.avgoverfreq = 'no';
        fcfg.frequency = LF;
        SourceData = ft_selectdata(fcfg, Fourier);
        
        %% source
        
        % create cfg
        %cfg.lamda = .5;
        for type = 1:length(params.type)
            sacfg              = [];
            sacfg.method       = lower(params.type{type});
            sacfg.headmodel    = vol;
            sacfg.elec         = elec_realigned ;
            sacfg.channel = data.label;
            %sacfg.(lower(cfg.type)).lambda = cfg.lamda;
            sacfg.(lower(params.type{type})) = params.(lower(params.type{type}));
            
            sacfg.(lower(params.type{type})).keepfilter   = 'yes';
            sacfg.(lower(params.type{type})).fixedori     = 'yes';
            sacfg.grid         = lf;
            sacfg.keepmom       = 'yes';
            sacfg.(lower(params.type{type})).projectnoise = 'no';
            
            sourceAll       = ft_sourceanalysis_checkless(sacfg, SourceData);
            
            
            %% Find Relevant Trials from that Condition info
            for Cond = 1:length(Deci.Analysis.Conditions)
                
                maxt = length(find(cellfun(@(c) any(ismember(Deci.Analysis.Conditions{Cond},c)), Deci.DT.Markers)));
                info.alltrials = find(sum(ismember(data.events,Deci.Analysis.Conditions{Cond}),2) == maxt);
                
                %% ignore all locks with missing nans
                if Deci.Analysis.IgnoreNanLocks
                    minamountofnans = min(mean(isnan(data.locks(info.alltrials,:)),2));
                    info.nanlocks = mean(isnan(data.locks(info.alltrials,:)),2) ~= minamountofnans;
                    
                    if any(info.nanlocks)
                        display(['ignoring ' num2str(length(find(info.nanlocks))) ' trials with missing locks'])
                    end
                else
                    info.nanlocks = logical(size(info.alltrials));
                end
                
                %% Reject Arts
                ccfg.trials =  info.alltrials(~info.nanlocks);
                if Deci.Analysis.ApplyArtReject
                    ccfg.trials = ccfg.trials(ismember(data.trlnum(ccfg.trials),data.postart.trlnum));
                    display('Applying Artifact Rejection')
                    display(['Final trial count is ' num2str(length(ccfg.trials))])
                else
                    display('Not Applying Artifact Rejection')
                end
                
                ccfg.trials = find(ismember(dat.trlnum,ccfg.trials));
                cond = ft_selectdata(ccfg,SourceData);
                
                cacfg = sacfg;
                cacfg.sourcemodel.filter = sourceAll.avg.filter;
                source = ft_sourceanalysis_checkless(cacfg, cond);
                
                mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Source' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{Cond}]);
                save([Deci.Folder.Analysis filesep 'Extra' filesep 'Source' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{Cond} filesep params.type{type} '_' params.foi{freq}],'source','-v7.3');
                
                
            end
        end
    end
    
end
end