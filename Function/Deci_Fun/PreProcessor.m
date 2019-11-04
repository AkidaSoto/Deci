function PreProcessor(Deci,subject_list)

disp('----------------------');
disp(['Starting PreProcessor for ' Deci.SubjectList{subject_list}]);
tic;
%Load and Set-up
cfg = [];
load([Deci.Folder.Definition filesep Deci.SubjectList{subject_list}],'cfg');
TrlDefs = cfg;



feedback = Deci.PP.feedback;
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

condinfo = {data_eeg.trialinfo TrlDefs.event TrlDefs.trialnum};
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
    disp('Implicit Rereference');
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
data.condinfo = condinfo;

if ~isempty(find(cellfun(@(c) any(any(isnan(c))), data.trial) == 1)) && Deci.PP.RejectNans
    nantrials = find(cellfun(@(c) any(any(isnan(c))), data.trial) == 1);
    Reject_cfg.trials = logical(ones([size(data.trial)]));
    Reject_cfg.trials(nantrials) = false;
    
    data = ft_selectdata(Reject_cfg,data);
    condinfo = data.condinfo;
    warning(['Found trial(s) containing nan in rawdata for ' Deci.SubjectList{subject_list} '. Revise Data and then use .RejectNans']);
elseif ~isempty(find(cellfun(@(c) any(any(isnan(c))), data.trial) == 1)) && ~Deci.PP.RejectNans
    error(['Found trial(s) containing nan in rawdata for ' Deci.SubjectList{subject_list} '. Revise Data and then use .RejectNans']);
end

preart   = condinfo;

disp(['Finished PreProcessor at ' num2str(toc)]);
disp('----------------------');
%% ICA
disp(['Starting ICA at ' num2str(toc)]);

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = Deci.ICA.bpfreq;
data_bp = ft_preprocessing(cfg,data);

cfg = [];
cfg.method  = 'runica';
cfg.numcomponent= 20;
cfg.feedback = feedback;
cfg.demean     = 'no';
data_musc = ft_componentanalysis(cfg, data_bp);

cfg           = [];
cfg.numcomponent= 20;
cfg.unmixing  =data_musc.unmixing;
cfg.topolabel = data_musc.topolabel;
cfg.feedback = feedback;
cfg.demean     = 'no';
data_all     = rmfield(ft_componentanalysis(cfg, data),'cfg');

figure;
cfg.component = [1:20];
cfg.viewmode = 'component';

clear cfg.method
cfg.channel = 'all';

cfg.component = [];


comps = [data_musc.trial{:}];
eyes = [data.trial{:}];
eyechan = eyes(ismember(data.label,Deci.ICA.eog),:);

for eye = 1:size(eyechan,1)
    for comp = 1:size(comps,1)
        [compcorr, p] = corrcoef(eyechan(eye,:),comps(comp,:));
        corr(eye,comp,1) = compcorr(1,2);
        corr(eye,comp,2) = p(1,2);
    end
    
    component{eye} = find(abs(corr(eye,:,1)) >= Deci.ICA.cutoff);
end



if ~Deci.ICA.Automatic
    
    cfg.component = [1:20];
    cfg.viewmode = 'component';
    cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
    
    cfg.channelcolormap = zeros(2,3);
    cfg.colorgroups = ones(20,1)+1;
    
    cfg.channelcolormap(1,1) = 1;
    cfg.colorgroups(unique([component{:}]),1) = 1;
    
    
    cfg.channel = 'all';
    
    fakeUI = figure;
    select_labels(fakeUI,[],sort(data_musc.label));
    fakeUI.Visible =  'off';
    ft_databrowser(cfg,data_all);
    suptitle(Deci.SubjectList{subject_list});
    waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
    
    if isempty(fakeUI.UserData)
        cfg.component = [];
    else
        cfg.component = find(ismember(data_musc.label,fakeUI.UserData));
    end
    close(fakeUI)
    corr = [];
    
%     if ~isempty(artf.artfctdef.visual.artifact)
%     
%        data_art = ft_rejectartifact(artf,data_all)
%         
%     end
else
    cfg.component = unique([component{:}]);
end

cfg.demean = 'yes';
data = ft_rejectcomponent(cfg, data_all);
data.unmixing = cfg.unmixing;
data.topolabel = cfg.topolabel;

data.condinfo = condinfo;
data.preart = preart;

data = rmfield(data,'cfg');

mkdir([Deci.Folder.Preproc])
save([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list}],'data','-v7.3')
data = rmfield(data,'trial');
save([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '_info'],'data','corr','-v7.3');

end