

function PreProcessor(Deci,subject_list)

disp('----------------------');
disp(['Starting PreProcessor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);

%Load and Set-up
cfg = [];
load([Deci.Folder.Definition filesep Deci.SubjectList{subject_list}],'cfg');
TrlDefs = cfg;

feedback = Deci.PP.feedback;

%% File check

if ~strcmpi(TrlDefs.datafile,[Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.eeg']) || ~strcmpi(TrlDefs.headerfile,[Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.vhdr']) || ~strcmpi(TrlDefs.dataset,[Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.eeg'])
TrlDefs.datafile = [Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.eeg'];
TrlDefs.headerfile = [Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.vhdr'];
TrlDefs.dataset = [Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.eeg'];
end

%% Detrending and Filtering the FullData
fullcfg = rmfield(TrlDefs,'trl');
evalc('data_eeg = ft_preprocessing(fullcfg)');

if ~isempty(Deci.PP.filter)
    filter_cfg = Deci.PP.filter;
    filter_cfg.feedback = feedback;
    evalc('data_eeg = ft_preprocessing(filter_cfg,data_eeg)');
    disp('Full Data Filtering Applied');
end

% Once Filtered, Apply trial definition to full data.
Trialcfg.trl = TrlDefs.trl;
evalc('data_eeg = ft_redefinetrial(Trialcfg,data_eeg)');

locks = data_eeg.trialinfo;
events = TrlDefs.event;
trlnum = TrlDefs.trialnum;

if Deci.PP.demean 
   dcfg.demean = 'yes';
   dcfg.baselinewindow = 'all';
   evalc('data_eeg = ft_preprocessing(dcfg,data_eeg)');
end



%%
if ~isempty(Deci.PP.ScalingFactor)
    disp('Data Scaled');
    data_eeg.trial = cellfun(@(c) c*Deci.PP.ScalingFactor,data_eeg.trial,'un',0);
end

if ~isempty(Deci.PP.Imp)
    Imp = strsplit(Deci.PP.Imp,':');
    
    if ~ismember(Imp{2},data_eeg.label)
        error('invalid Implicit channels for reference')
    end
    Imp_cfg.reref = 'yes';
    Imp_cfg.channel  = 'all';
    Imp_cfg.implicitref = Imp{1};
    Imp_cfg.refchannel = Imp;
    Imp_cfg.feedback = feedback;
    evalc('data_eeg = ft_preprocessing(Imp_cfg,data_eeg)');
    disp('---Implicit Rereference applied---');
end

% if ~isempty(Deci.PP.Ocu)
%     Ocu = cellfun(@(c) strsplit(c,':'), strsplit(Deci.PP.Ocu,','),'UniformOutput',false);
%     allOcu = [Ocu{:}];
%
%     cfg = [];
%
%     if ~all(ismember(allOcu,data_eeg.label))
%         error('invalid ocular channels for reference')
%     end
%
%     for i = 1:length(Ocu)
%         cfg.channel = Ocu{i};
%         cfg.refchannel = Ocu{i}(1);
%         cfg.feedback = feedback;
%         evalc('data_eog(i) = ft_preprocessing(cfg,data_eeg)');
%         Hcfg.channel = Ocu{i}(2);
%         cfg.feedback = feedback;
%         evalc('data_eog(i)   = ft_preprocessing(Hcfg, data_eog(i))'); % nothing will be done, only the selection of the interesting channel
%     end
%
%     cfg.channel = [{'all'} arrayfun(@(c) strjoin(['-' c],''),allOcu,'un',0)] ;
%     cfg.feedback = feedback;
%     evalc('data_noeog = ft_selectdata(cfg,data_eeg)');
%
%     arraydata = arrayfun(@(c) {c},[data_noeog, data_eog]);
%     clear data_noeog data_eog
%     cfg.feedback = feedback;
%     evalc('data_eeg = ft_appenddata([],arraydata{:})');
%     clear arraydata
%     disp('Ocular Rereference');
% end

if Deci.PP.CleanLabels
    if ~isempty( data_eeg.label(strcmp(data_eeg.label,'OL')))
        data_eeg.label(strcmp(data_eeg.label,'OL')) = {'PO7'};
    end
    if  ~isempty( data_eeg.label(strcmp(data_eeg.label,'OR')))
        data_eeg.label(strcmp(data_eeg.label,'OR')) = {'PO8'};
    end
    if ~isempty( data_eeg.label(strcmp(data_eeg.label,'T3')))
        data_eeg.label(strcmp(data_eeg.label,'T3')) = {'T7'};
    end
    if  ~isempty( data_eeg.label(strcmp(data_eeg.label,'T4')))
        data_eeg.label(strcmp(data_eeg.label,'T4')) = {'T8'};
    end
    if ~isempty( data_eeg.label(strcmp(data_eeg.label,'T5')))
        data_eeg.label(strcmp(data_eeg.label,'T5')) = {'P7'};
    end
    if  ~isempty( data_eeg.label(strcmp(data_eeg.label,'T6')))
        data_eeg.label(strcmp(data_eeg.label,'T6')) = {'P8'};
    end
    if  ~isempty( data_eeg.label(strcmp(data_eeg.label,'VEM')))
        data_eeg.label(strcmp(data_eeg.label,'VEM')) = {'BVEOG'};
    end
    if  ~isempty( data_eeg.label(strcmp(data_eeg.label,'HEM')))
        data_eeg.label(strcmp(data_eeg.label,'HEM')) = {'RHEOG'};
    end
    if  ~isempty( data_eeg.label(strcmp(data_eeg.label,'LM')))
        data_eeg.label(strcmp(data_eeg.label,'LM')) = {'TP9'};
    end
    if  ~isempty( data_eeg.label(strcmp(data_eeg.label,'RM')))
        data_eeg.label(strcmp(data_eeg.label,'RM')) = {'TP10'};
    end
    
    if  ~isempty( data_eeg.label(strcmp(data_eeg.label,'32')))
        rm32.channel = data_eeg.label(~strcmp(data_eeg.label,'32'));
        data_eeg = ft_selectdata(rm32,data_eeg);
    end
end

if  ~isempty(data_eeg.label(strcmp(data_eeg.label,'StimTrak')))
    rm32.channel = data_eeg.label(~strcmp(data_eeg.label,'StimTrak'));
    data_eeg = ft_selectdata(rm32,data_eeg);
    
    disp('StimTrak Removed')
end

% if ~isempty(Deci.PP.DownSample)
%     data_eeg = ft_resampledata(struct('resamplefs',Deci.PP.DownSample,'detrend','no'),data_eeg);
% end

if ~isempty(Deci.PP.More)
    More_cfg = Deci.PP.More;
    More_cfg.feedback = feedback;
    evalc('data_eeg = ft_preprocessing(More_cfg,data_eeg)');
    disp('Additional Preprocessing');
end

data = data_eeg;
data.locks = locks;
data.events = events;
data.trlnum = trlnum;

if ~isempty(find(cellfun(@(c) any(any(isnan(c))), data.trial) == 1)) && Deci.PP.RejectNans
    nantrials = find(cellfun(@(c) any(any(isnan(c))), data.trial) == 1);
    Reject_cfg.trials = logical(ones([size(data.trial)]));
    Reject_cfg.trials(nantrials) = false;
    
    data = ft_selectdata(Reject_cfg,data);
    
    locks = data.locks;
    events = data.events;
    trlnum = data.trlnum;
    
    warning(['Found trial(s) containing nan in rawdata for ' Deci.SubjectList{subject_list} '. Revise Data and then use .RejectNans']);
elseif ~isempty(find(cellfun(@(c) any(any(isnan(c))), data.trial) == 1)) && ~Deci.PP.RejectNans
    error(['Found trial(s) containing nan in rawdata for ' Deci.SubjectList{subject_list} '. Revise Data and then use .RejectNans']);
end

 %% Manual Trial Rejection
 
if Deci.PP.Manual_Trial_Rejection
    cfg =[];
    cfg.method = 'trial';
    cfg.alim = 100;
    tcfg.toilim = [abs(nanmax(locks,[],2)/1000)+Deci.Art.crittoilim(1) abs(nanmin(locks,[],2)/1000)+Deci.Art.crittoilim(2)]; 
    evalc('data_rej = ft_rejectvisual(cfg,ft_redefinetrial(tcfg,data))');
    
    postart.locks = locks(logical(ismember(data.trlnum,data_rej.saminfo)),:);
    postart.events = events(logical(ismember(postart.trlnum,data_rej.saminfo)),:);
    postart.trlnum = trlnum(logical(ismember(postart.trlnum,data_rej.saminfo)),:);
    
    display(' ')
    display('---Manual Trial Rejection Applied---')
    display(['Rejected ' num2str(length(find(~logical(data_rej.saminfo)))) ' trials'])
    display(['Remaining ' num2str(length(postart.trlnum)) ' trials'])
    display(['Remaining ' num2str([length(postart.trlnum)/length(trlnum)]*100) '% trials'])
    display(' ')
    pause(.05);
else
    postart.locks = locks;
    postart.events = events;
    postart.trlnum = trlnum;
end

%% ICA

display(' ')

if Deci.ICA.do

disp(['---Starting ICA---']);
display(' ')
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = Deci.ICA.bpfreq;
evalc('data_bp = ft_preprocessing(cfg,data)');

cfg = [];
cfg.method  = 'runica';
cfg.numcomponent= Deci.ICA.CompNum;
cfg.feedback = feedback;
cfg.demean     = 'no';
data_musc = ft_componentanalysis(cfg, data_bp);

cfg           = [];
cfg.numcomponent= Deci.ICA.CompNum;
cfg.unmixing  =data_musc.unmixing;
cfg.topolabel = data_musc.topolabel;

end

data.locks = locks;
data.events = events;
data.trlnum = trlnum;
data.postart = postart;

data = rmfield(data,'cfg');

mkdir([Deci.Folder.Preproc])
save([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list}],'data','cfg','-v7.3')

disp(['----------------------']);
end