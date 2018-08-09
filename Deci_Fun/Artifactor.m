function Artifactor(Deci)

for subject_list = 1:length(Deci.SubjectList)
    

    
    %% ICA
    
    if Deci.Art.ICA
        
        data = [];
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);

        
        cfg = [];
        cfg.channel =  'all';
        
        cfg.method  = 'runica';
        cfg.runica.pca = 20;
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
%% Trial Rejection
for subject_list = 1:length(Deci.SubjectList)
    if ~isempty(Deci.Art.TR.Eye) || ~isempty(Deci.Art.TR.Muscle)
        
        data = [];
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);

        
        cfg = [];
        load([Deci.Folder.Definition filesep Deci.SubjectList{subject_list}],'cfg');
        
%         testcfg.latency = Deci.Art.TR.Toi;
%         EEG_data = ft_selectdata(testcfg,data);
        
        redefine = 0;
        if exist([Deci.Folder.Version  filesep 'Redefine' filesep Deci.SubjectList{subject_list} '.mat']) == 2
            redefine = 1;
            
            retrl = [];
            load([Deci.Folder.Version  filesep 'Redefine' filesep Deci.SubjectList{subject_list} '.mat']);
            cfg.offset = retrl;
            cfg.shift = Deci.Art.TR.rToi(1);
            
            EEG_redefine = ft_datashift(cfg,data);
%             testcfg.latency = Deci.Art.TR.rToi;
%             EEG_redefine = ft_selectdata(testcfg,EEG_redefine);
        end
        
        cfg.channel = 'all';
        cfg.datafile = filesyntax(cfg.datafile);
        cfg.headerfile = filesyntax(cfg.headerfile);
        cfg.dataset = filesyntax(cfg.dataset);
        
        artif = [];
        
        if ~isempty(Deci.Art.TR.Eye)
            
            EOGpadval = [0 .1 0];
            EOGcutoffval = 4;
            EOGbpval = [1 15];
            
            % EOG artifact detection
            acfg = cfg;
            acfg.artfctdef.zvalue.channel     =  Deci.Art.TR.Eye.Chans;
            acfg.artfctdef.zvalue.cutoff      =  EOGcutoffval;
            acfg.artfctdef.zvalue.trlpadding  =  EOGpadval(1);
            acfg.artfctdef.zvalue.artpadding  =  EOGpadval(2);
            acfg.artfctdef.zvalue.fltpadding  =  EOGpadval(3);
            
            % algorithmic parameters
            acfg.artfctdef.zvalue.bpfilter   = 'yes';
            acfg.artfctdef.zvalue.bpfilttype = 'but';
            acfg.artfctdef.zvalue.bpfreq     = [EOGbpval(1) EOGbpval(2)];
            acfg.artfctdef.zvalue.bpfiltord  = 4;
            acfg.artfctdef.zvalue.hilbert    = 'yes';
            
            % feedback
            if Deci.Art.TR.Eye.Interactive
                acfg.artfctdef.zvalue.interactive = 'yes';
            else
                acfg.artfctdef.zvalue.interactive = 'no';
            end
            
            
            [~, artifact_EOG_cond] = ft_artifact_zvalue_checkless(acfg,data);
            
            if ~isempty(artifact_EOG_cond)
                ecfg = cfg;
                ecfg.artfctdef.reject = 'complete';
                ecfg.artfctdef.eog.artifact = artifact_EOG_cond;
                ecfg.artfctdef.crittoilim = Deci.Art.TR.Toi;
                
                artifact = ft_rejectartifact(ecfg,data);
                artif = [artif;artifact.saminfo];
            end
            
            if redefine
                [~, artifact_EOGr_cond] = ft_artifact_zvalue_checkless(acfg,EEG_redefine);
                
                if ~isempty(artifact_EOGr_cond)
                    ercfg = cfg;
                    ercfg.artfctdef.reject = 'complete';
                    ercfg.artfctdef.eog.artifact = artifact_EOGr_cond;
                    ercfg.artfctdef.crittoilim = Deci.Art.TR.rToi;
                    
                    artifact = ft_rejectartifact(ercfg,EEG_redefine);
                    artif = [artif;artifact.saminfo];
                end
            end
            
        end
        
        if ~isempty(Deci.Art.TR.Muscle)
            
            hdr = ft_read_header(cfg.headerfile);
            Nyq = [hdr.Fs/2] -1;
            
            MUSpadval = [0 .1 0];
            MUScutoffval = 25;
            MUSbpval = [110 Nyq];
            
            bcfg = cfg;
            bcfg.artfctdef.zvalue.channel     = 'all';
            bcfg.artfctdef.zvalue.cutoff      = MUScutoffval;
            bcfg.artfctdef.zvalue.trlpadding   =MUSpadval(1);
            bcfg.artfctdef.zvalue.artpadding  = MUSpadval(2);
            bcfg.artfctdef.zvalue.fltpadding  = MUSpadval(3);
            
            % algorithmic parameters
            bcfg.artfctdef.zvalue.bpfilter              = 'yes';
            bcfg.artfctdef.zvalue.bpfreq                = MUSbpval;
            bcfg.artfctdef.zvalue.bpfiltord             = 1;
            bcfg.artfctdef.zvalue.bpfilttype            = 'but';
            bcfg.artfctdef.zvalue.hilbert               = 'yes';
            bcfg.artfctdef.zvalue.boxcar                = 0.2;
            bcfg.artfctdef.zvalue.bpinstabilityfix      = 'reduce';
            
            % make the process interactive
            if Deci.Art.TR.Muscle.Interactive
                bcfg.artfctdef.zvalue.interactive = 'yes';
            else
                bcfg.artfctdef.zvalue.interactive = 'no';
            end
            
            
            
            [~, artifact_muscle_cond] = ft_artifact_zvalue_checkless(bcfg,data);
            
            if ~isempty(artifact_muscle_cond)
                mcfg =cfg;
                mcfg.artfctdef.reject = 'complete';
                mcfg.artfctdef.muscle.artifact = artifact_muscle_cond;
                mcfg.artfctdef.crittoilim = Deci.Art.TR.Toi;
                
                artifact = ft_rejectartifact(mcfg,data);
                artif = [artif;artifact.saminfo];
            end
            
            if redefine
                [~, artifact_muscler_cond] = ft_artifact_zvalue_checkless(bcfg,EEG_redefine);
                
                if ~isempty(artifact_muscler_cond)
                    emcfg = cfg;
                    emcfg.artfctdef.reject = 'complete';
                    emcfg.artfctdef.eog.artifact = artifact_muscler_cond;
                    emcfg.artfctdef.crittoilim = Deci.Art.TR.rToi;
                    
                    artifact = ft_rejectartifact(emcfg,EEG_redefine);
                    artif = [artif;artifact.saminfo];
                end
            end
            
            
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
        
        logic = find(ismember(data.sampleinfo(:,1),dummy.sampleinfo(:,1)));
        condition = data.trialinfo(ismember(data.sampleinfo(:,1),dummy.sampleinfo(:,1)));
        
        conds =  unique( data.trialinfo,'stable');
        
        for con = 1:length(conds)
            
            artifacts = zeros(size(data.trialinfo(data.trialinfo == conds(con))))';
            artifacts(logic(condition == conds(con))) = true;
            save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} filesep num2str(conds(con))],'artifacts');
        end
        
    end
    
end


end
