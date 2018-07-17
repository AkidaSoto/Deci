function PreProcessor(Deci)
    
for subject_list = 1:length(Deci.SubjectList)
     
    cfg = load([Deci.Folder.Definition filesep Deci.SubjectList{subject_list}]);
    cfg = cfg.cfg;
    
    cfg.datafile = filesyntax(cfg.datafile);
    cfg.headerfile = filesyntax(cfg.headerfile);
    cfg.dataset = filesyntax(cfg.dataset);
    
    
    
     data_eeg = ft_preprocessing(cfg);
    
    if ~isempty(Deci.PP.ScalingFactor)
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
        data_eeg = ft_preprocessing(cfg);
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
            data_eog(i) = ft_preprocessing(cfg,data_eeg);
            Hcfg.channel = Ocu{i}(2);
            data_eog(i)   = ft_preprocessing(Hcfg, data_eog(i)); % nothing will be done, only the selection of the interesting channel
        end
        
        cfg.channel = [{'all'} arrayfun(@(c) strjoin(['-' c],''),allOcu,'un',0)] ;
        data_noeog = ft_selectdata(cfg,data_eeg);

        arraydata = arrayfun(@(c) {c},[data_noeog, data_eog]);
        data = ft_appenddata([],arraydata{:});
    else
        
        data = data_eeg;
    end
   
    if ~isempty(Deci.PP.hbp)
        cfg =[];
        cfg.hpfreq = Deci.PP.hbp;
        cfg.hpfilter      = 'yes';
        data = ft_preprocessing(cfg,data);
    end
    
     if ~isempty(Deci.PP.Demean)
    cfg = [];
    cfg.demean = 'yes';
    cfg.baselinewindow = Deci.PP.Demean;
    data = ft_preprocessing(cfg,data);
     end
     
    [data.trialinfo,i] = sort(floor(data.trialinfo));
    data.sampleinfo = data.sampleinfo(i,:);
    data.trial = data.trial(i);
    
    mkdir([Deci.Folder.Preproc])
    label = data.label;
    save([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list}],'data','label','-v7.3');
    
end


end