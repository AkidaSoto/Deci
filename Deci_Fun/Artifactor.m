function Artifactor(Deci)

for subject_list = 1:length(Deci.SubjectList)
    
    if Deci.Art.ICA  || ~isempty(Deci.Art.TR.Eye) || ~isempty(Deci.Art.TR.Muscle)
        data = load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list}]);
        label = data.label;
        data = data.data;
        
    end
    
    %% ICA
    
    if Deci.Art.ICA
        
        cfg = [];
        cfg.channel = label;
        
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
        ft_databrowser(cfg, datacomp)
        
        
        cfg.component = input('Components to reject? (Delimit by space)','s');
        
        if isempty(cfg.component)
            cfg.component = [];
        else
            cfg.component = str2num(cfg.component);
            if any(cfg.component <= 0 & cfg.component >= 20)
                error('component numbers outside of limits of 1-20');
            end
        end
        
        data = ft_rejectcomponent(cfg, datacomp);
        clear datacomp
        
        mkdir([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list}])
        save([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list}],'data','label','-v7.3')
    end
    
    %% Trial Rejection
    
    if ~isempty(Deci.Art.TR.Eye) || ~isempty(Deci.Art.TR.Muscle)
        cfg = [];
        load([Deci.Folder.Definition filesep Deci.SubjectList{subject_list}],'cfg');
        
        
        cfg.datafile = filesyntax(cfg.datafile);
        cfg.headerfile = filesyntax(cfg.headerfile);
        cfg.dataset = filesyntax(cfg.dataset);
        
        artif = [];
        
        if ~isempty(Deci.Art.TR.Eye)
            
            EOGpadval = [0 .1 0];
            EOGcutoffval = 6;
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
            
            testcfg.latency = Deci.Art.TR.Eye.Toi;
            EOG_data = ft_selectdata(testcfg,data);
            
            [~, artifact_EOG_cond] = ft_artifact_zvalue_checkless(acfg,EOG_data);
            
            if ~isempty(artifact_EOG_cond)
                ecfg = cfg;
                ecfg.artfctdef.reject = 'complete';
                ecfg.artfctdef.eog.artifact = artifact_EOG_cond;
                
                artifact = ft_rejectartifact(ecfg,EOG_data);
                artif = [artif;artifact.saminfo];
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
            
            testcfg.latency = Deci.Art.TR.Muscle.Toi;
            MUS_data = ft_selectdata(testcfg,data);
            
            [~, artifact_muscle_cond] = ft_artifact_zvalue_checkless(bcfg,MUS_data);
            
            if ~isempty(artifact_muscle_cond)
                mcfg =cfg;
                mcfg.artfctdef.reject = 'complete';
                mcfg.artfctdef.muscle.artifact = artifact_muscle_cond;
                
                artifact = ft_rejectartifact(mcfg,MUS_data);
                artif = [artif;artifact.saminfo];
            end
            
        end
        
       
        
        if ~isempty(artif)
            conds =  unique(floor(cfg.trl(:,4)),'stable');
            
            mkdir([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}]);
            artif = all(artif,1);
            
            for con = 1:length(conds)
                artifacts = artif(floor(data.trialinfo) ==conds(con));
                save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} filesep num2str(conds(con))],'artifacts');
            end
        end
        
    end
    
end


end
