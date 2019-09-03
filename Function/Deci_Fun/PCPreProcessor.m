function PCPreProcessor(Deci,subject_list)

disp('----------------------');
disp(['Starting PreProcessor for ' Deci.SubjectList{subject_list}]);
tic;

feedback  = 'none';

cfg = load([Deci.Folder.Definition filesep Deci.SubjectList{subject_list}]);
cfg = cfg.cfg;

cfg.datafile = filesyntax(cfg.datafile);
cfg.headerfile = filesyntax(cfg.headerfile);
cfg.dataset = filesyntax(cfg.dataset);

cfg.feedback = feedback;
evalc('data_eeg = ft_preprocessing(cfg)');

condinfo = {data_eeg.trialinfo cfg.event cfg.trialnum};

% for test = 1:length(data_eeg.trial)
%     data_eeg.trial{test}(10,:) = 0;
%     data_eeg.trial{test}(end,:) = 0;
% end

if ~isempty(Deci.PP.ScalingFactor)
    disp('Data Scaled');
    data_eeg.trial = cellfun(@(c) c*Deci.PP.ScalingFactor,data_eeg.trial,'un',0);
end


if ~isempty(Deci.PP.Imp)
    Imp = strsplit(Deci.PP.Imp,':');
    
    if ~ismember(Imp{2},data_eeg.label)
        error('invalid Implicit channels for reference')
    end
    
    cfg.reref = 'yes';
    cfg.channel  = 'all';
    cfg.implicitref = Imp{1};
    cfg.refchannel = Imp;
    cfg.feedback = feedback;
    evalc('data_eeg = ft_preprocessing(cfg,data_eeg)');
    disp('Implicit Rereference');
end

if ~isempty(Deci.PP.Ocu)
    Ocu = cellfun(@(c) strsplit(c,':'), strsplit(Deci.PP.Ocu,','),'UniformOutput',false);
    allOcu = [Ocu{:}];
    
    cfg = [];
    
    if ~all(ismember(allOcu,data_eeg.label))
        error('invalid ocular channels for reference')
    end
    
    for i = 1:length(Ocu)
        cfg.channel = Ocu{i};
        cfg.refchannel = Ocu{i}(1);
        cfg.feedback = feedback;
        evalc('data_eog(i) = ft_preprocessing(cfg,data_eeg)');
        Hcfg.channel = Ocu{i}(2);
        cfg.feedback = feedback;
        evalc('data_eog(i)   = ft_preprocessing(Hcfg, data_eog(i))'); % nothing will be done, only the selection of the interesting channel
    end
    
    cfg.channel = [{'all'} arrayfun(@(c) strjoin(['-' c],''),allOcu,'un',0)] ;
    cfg.feedback = feedback;
    evalc('data_noeog = ft_selectdata(cfg,data_eeg)');
    
    arraydata = arrayfun(@(c) {c},[data_noeog, data_eog]);
    clear data_noeog data_eog
    cfg.feedback = feedback;
    evalc('data_eeg = ft_appenddata([],arraydata{:})');
    clear arraydata
    disp('Ocular Rereference');
end

if ~isempty(Deci.PP.HBP)
    cfg =[];
    cfg.hpfreq = Deci.PP.HBP;
    cfg.hpfilter      = 'yes';
    cfg.feedback = feedback;
    evalc('data_eeg = ft_preprocessing(cfg,data_eeg)');
    disp('Highbandpass Filter');
end
% 
% if ~isempty(Deci.PP.Demean)
%     cfg = [];
%     cfg.demean = 'yes';
%     cfg.baselinewindow = Deci.PP.Demean;
%     cfg.feedback = feedback;
%     evalc('data_eeg = ft_preprocessing(cfg,data_eeg)');
%     disp('Baseline Correction');
% end


if ~isempty(Deci.PP.Repair)
    vcfg.viewmode = 'vertical';
    
    repaircheck = 1;
    while repaircheck
        
        if strcmpi(Deci.PP.Repair.Type,'Manual')
            fakeUI = axes;
            select_labels(fakeUI,[],sort(data_eeg.label));
            fakeUI.Parent.Visible =  'off';
            ft_databrowser(vcfg,data_eeg);
            suptitle(Deci.SubjectList{subject_list});
            waitfor(gcf,'BeingDeleted','on');
            
            if isempty(fakeUI.UserData)
                repaircheck = 0;
                disp('Satisfied with all channels, continuing...')
            else
                
                scfg.badchannel = fakeUI.UserData;
                neighbours       = load('easycap_neighbours','neighbours');
                scfg.neighbours = neighbours.neighbours;
                scfg.method = 'spline';
                scfg.elec = ft_read_sens('standard_1020.elc');
                
                eyecfg.channel = data_eeg.label(~ismember(data_eeg.label,scfg.elec.label));
                eyedata = ft_selectdata(eyecfg,data_eeg);
                data_eeg = ft_channelrepair(scfg,data_eeg);
                
                data_eeg = ft_appenddata([],eyedata,data_eeg);
                disp('Fixed Channels, Check new data if statisfied')
            end
            
            close(fakeUI.Parent)
        else
            if ~isempty(Deci.PP.Repair.Auto{subject_list})
                
                scfg.badchannel = Deci.PP.Repair.Auto{subject_list};
                neighbours       = load('easycap_neighbours','neighbours');
                scfg.neighbours = neighbours.neighbours;
                scfg.method = 'spline';
                scfg.elec = ft_read_sens('standard_1020.elc');
                
                eyecfg.channel = data_eeg.label(~ismember(data_eeg.label,scfg.elec.label));
                eyedata = ft_selectdata(eyecfg,data_eeg);
                data_eeg = ft_channelrepair(scfg,data_eeg);
                
                data_eeg = ft_appenddata([],eyedata,data_eeg);
                
            end
            clear eyedata
            repaircheck = 0;
        end
        
    end
end


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
    
    %
    %          if ~isempty( layout.label(strcmp(layout.label,'PO7')))
    %             layout.label(strcmp(layout.label,'PO7')) = {'OL'};
    %         end
    %         if  ~isempty( layout.label(strcmp(layout.label,'PO8')))
    %             layout.label(strcmp(layout.label,'PO8')) = {'OR'};
    %         end
    %
    %         if ~isempty( layout.label(strcmp(layout.label,'T7')))
    %             layout.label(strcmp(layout.label,'T7')) = {'T3'};
    %         end
    %         if  ~isempty( layout.label(strcmp(layout.label,'T8')))
    %             layout.label(strcmp(layout.label,'T8')) = {'T4'};
    %         end
    %
    %         if ~isempty( layout.label(strcmp(layout.label,'P7')))
    %             layout.label(strcmp(layout.label,'P7')) = {'T5'};
    %         end
    %         if  ~isempty( layout.label(strcmp(layout.label,'P8')))
    %             layout.label(strcmp(layout.label,'P8')) = {'T6'};
    %         end
    % if  ~isempty( layout.label(strcmp(layout.label,'BVEOG')))
    %     layout.label(strcmp(layout.label,'BVEOG')) = {'VEM'};
    % end
    % if  ~isempty( layout.label(strcmp(layout.label,'RHEOG')))
    %     layout.label(strcmp(layout.label,'RHEOG')) = {'HEM'};
    % end
    %         if  ~isempty( layout.label(strcmp(layout.label,'TP9')))
    %             layout.label(strcmp(layout.label,'TP9')) = {'LM'};
    %         end
    %         if  ~isempty( layout.label(strcmp(layout.label,'TP10')))
    %             layout.label(strcmp(layout.label,'TP10')) = {'RM'};
    %         end
    %
    
end


if isfield(Deci.PP,'HemiFlip')
    if ~isempty(Deci.PP.HemiFlip)
        eachcond = unique(Deci.PP.Conditions);
        for k = 1:length(eachcond)
            
            cond = (ismember(Deci.PP.Conditions,eachcond(k)));
            condvalue = find(Deci.PP.HemiFlip & cond);
            
            if ~isempty(condvalue)
                
                hemi.trials = any(data_eeg.trialinfo == [eachcond(k)+find(ismember(find(cond),condvalue))*.01],2);
                
                nonhemi.trials = all(data_eeg.trialinfo ~= [eachcond(k)+find(ismember(find(cond),condvalue))*.01],2);
                
                hemidata = ft_selectdata(hemi,data_eeg);
                hemidata = hemifieldflip(hemidata);
                
                nonhemidata = ft_selectdata(nonhemi,data_eeg);
                
                data_eeg = ft_appenddata([],hemidata,nonhemidata);
            end
            
        end
    end
end

% 
% [data_eeg.trialinfo,i] = sort(floor(data_eeg.trialinfo(:,1)));
% data_eeg.sampleinfo = data_eeg.sampleinfo(i,:);
% data_eeg.trial = data_eeg.trial(:,i);
% data_eeg.time = data_eeg.time(:,i);


if ~isempty(Deci.PP.DownSample)
    data_eeg = ft_resampledata(struct('resamplefs',Deci.PP.DownSample,'detrend','no'),data_eeg);
end

if ~isempty(Deci.PP.More)
    cfg = Deci.PP.More;
    cfg.feedback = feedback;
    evalc('data_eeg = ft_preprocessing(cfg,data_eeg)');
    disp('Additional Preprocessing');
end

data = data_eeg;
data.condinfo = condinfo;


disp('Saving Preprocessing Data');
mkdir([Deci.Folder.Preproc])

% for test = 1:length(data.trial)
%     
%     data.trial{test}(find(ismember(data.label,'FCz')),:) = 0;
% end

data = rmfield(data,'cfg');

save([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list}],'data','-v7.3');
%save([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list}],'data');


disp(['Finished PreProcessor at ' num2str(toc)]);
disp('----------------------');

end