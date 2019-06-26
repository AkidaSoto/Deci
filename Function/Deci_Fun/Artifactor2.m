function Artifactor2(Deci)

%% Trial Rejection
for subject_list = 1:length(Deci.SubjectList)
    if ~isempty(Deci.Art.TR.Eye) || ~isempty(Deci.Art.TR.Muscle)
        
        data = [];
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
        
           
        cfg = [];
        cfg.demean = 'yes';
        cfg.baselinewindow = Deci.PP.Demean;
        data_eeg = ft_preprocessing(cfg,data_eeg);
        
        redefine = 0;
        if exist([Deci.Folder.Version  filesep 'Redefine' filesep Deci.SubjectList{subject_list} '.mat']) == 2
            redefine = 1;
            
            retrl = [];
            load([Deci.Folder.Version  filesep 'Redefine' filesep Deci.SubjectList{subject_list} '.mat']);
            cfg = [];
            cfg.offset = retrl;
            cfg.shift = Deci.Art.TR.rToi(1);
            
            EEG_redefine = ft_datashift(cfg,data);
            %             testcfg.latency = Deci.Art.TR.rToi;
            %             EEG_redefine = ft_selectdata(testcfg,EEG_redefine);
        end
        
        
        artif = [];
        
        cfg = [];
        cfg.preproc.bpfilter    = 'yes';
        cfg.preproc.bpfreq      = [110 140];
        cfg.preproc.bpfiltord   =  8;
        cfg.preproc.bpfilttype  = 'but';
        cfg.preproc.rectify     = 'yes';
        cfg.preproc.boxcar      = 0.2;
        cfg.latency = Deci.Art.TR.Toi;
        cfg.method = 'summary';
        
        
        [mus_data] = ft_rejectvisual(cfg, data);
        
        if ~all(mus_data.saminfo)
            artif = [artif;mus_data.saminfo];
        end
        
        if redefine
            cfg.latency = Deci.Art.TR.rToi;
            [rmus_data] = ft_rejectvisual(cfg, data);
            
            if ~all(rmus_data.saminfo)
                artif = [artif;rmus_data.saminfo];
            end
        end
        
        artif = all(artif,1);
        
        if Deci.Art.TR.View
        cfg = [];
        cfg.viewmode = 'vertical';
        cfg.trl = ~artif;
        cfg = ft_databrowser(cfg,data);
        uiwait(cfg.h);
        end
        
        if ~isempty(artif)
            conds =  unique(floor(cfg.trl(:,4)),'stable');
            
            mkdir([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}]);
            artif = all(artif,1);
            
            for con = 1:length(conds)
                artifacts = artif(floor(data.trialinfo) ==conds(con));
                arties(subject_list,con,:) =[length(find(artifacts)) length(artifacts)] ;
                save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} filesep num2str(conds(con))],'artifacts');
            end
        end
        
    end
end


for subject_list = 1:length(Deci.SubjectList)
    
    
    
    %% ICA
    
    if Deci.Art.ICA
        
        data = [];
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
        if exist([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} filesep num2str(conds(con)) '.mat'])
            artifacts = [];
            load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} filesep num2str(conds(con))],'artifacts');
        else
            artifacts = ones([size(data.trial));
        end
        
        cfg = [];
        cfg.channel =  'all';
        
        cfg.method  = 'runica';
        cfg.numcomponent = 20;
        cfg.trials = artifacts;
        datacomp = ft_componentanalysis(cfg, data);
        
        figure;
        cfg.component = [1:20];
        cfg.viewmode = 'component';
        
        cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
        
        cfg.comment   = 'no';
        ft_topoplotIC(cfg, datacomp);
        
        clear cfg.method
        
        cfg.channel = 'all';
        
        fakeUI = figure;
        select_labels(fakeUI,[],sort(datacomp.label));
        fakeUI.Visible =  'off';
        ft_databrowser(cfg,datacomp);
        suptitle(Deci.SubjectList{subject_list});
        waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
        
        if isempty(fakeUI.UserData)
            cfg.component = [];
        else
            cfg.component = find(ismember(fakeUI.UserData,datacomp.label));
            if any(cfg.component <= 0 & cfg.component >= 20)
                error('component numbers outside of limits of 1-20');
            end
        end
        close(fakeUI)
        
        data = ft_rejectcomponent(cfg, datacomp);
        clear datacomp
        
        mkdir([Deci.Folder.Preproc])
        save([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list}],'data','-v7.3')
    end
    
end


for subject_list = 1:length(Deci.SubjectList)
    
    if ~isempty(Deci.Art.Manual)
        
        data = [];
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
        artif = [];
        
        testcfg.latency = Deci.Art.TR.Toi;
        EEG_data = ft_selectdata(testcfg,data);
        
        redefine = 0;
        if exist([Deci.Folder.Version  filesep 'Redefine' filesep Deci.SubjectList{subject_list} '.mat']) == 2
            redefine = 1;
            
            retrl = [];
            load([Deci.Folder.Version  filesep 'Redefine' filesep Deci.SubjectList{subject_list} '.mat']);
            cfg.offset = retrl;
            cfg.shift = Deci.Art.Manual.rToi(1);
            
            EEG_redefine = ft_datashift(cfg,data);
            testcfg.latency = Deci.Art.Manual.rToi;
            EEG_redefine = ft_selectdata(testcfg,EEG_redefine);
        end
        
        if isdir([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}])
            
            ArtFolders = CleanDir([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}]);
            artifacts = [];
            for Arts = 1:length(ArtFolders)
                load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} filesep ArtFolders{Arts}],'artifacts');
                
                artif = [artif artifacts];
            end
            
        end
        
        if redefine
            for trl = 1:length(data.trial)
                data.trial{trl} = cat(2,EEG_data.trial{trl},EEG_redefine.trial{trl});
                data.time{trl} = cat(2,EEG_data.time{trl},EEG_redefine.time{trl});
            end
        else
            
            data = EEG_data;
        end
        
        cfg          = [];
        cfg.method   = 'summary';
        cfg.trials = find(artif);
        dummy        = ft_rejectvisual(cfg,data);
        
        logic = ismember(data.sampleinfo(:,1),dummy.sampleinfo(:,1));
        
        conds =  unique( data.trialinfo,'stable');
        
        for con = 1:length(conds)
            
            artifacts = zeros(size(data.trialinfo(data.trialinfo == conds(con))))';
            artifacts(logic(data.trialinfo == conds(con))) = true;
            save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} filesep num2str(conds(con))],'artifacts');
        end
        
    end
    
end


end
