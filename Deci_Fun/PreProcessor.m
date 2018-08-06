function PreProcessor(Deci)

savelog = parfeval('',0);

for subject_list = 1:length(Deci.SubjectList)
    
    finishedjobs = find(ismember({savelog.State},'finished'));
    delete(savelog(finishedjobs));
    savelog = savelog(~finishedjobs);
    
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
        data_eeg = ft_preprocessing(cfg,data_eeg);
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
        clear data_noeog data_eog
        data_eeg = ft_appenddata([],arraydata{:});
        clear arraydata
    end
    
    if ~isempty(Deci.PP.HBP)
        cfg =[];
        cfg.hpfreq = Deci.PP.HBP;
        cfg.hpfilter      = 'yes';
        data_eeg = ft_preprocessing(cfg,data_eeg);
    end
    
    if ~isempty(Deci.PP.Demean)
        cfg = [];
        cfg.demean = 'yes';
        cfg.baselinewindow = Deci.PP.Demean;
        data_eeg = ft_preprocessing(cfg,data_eeg);
    end
    
    
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
    
    
    [data_eeg.trialinfo,i] = sort(floor(data_eeg.trialinfo(:,1)));
    data_eeg.sampleinfo = data_eeg.sampleinfo(i,:);
    data_eeg.trial = data_eeg.trial(:,i);
    data_eeg.time = data_eeg.time(:,i);
     
    mkdir([Deci.Folder.Preproc])
    savelog(subject_list) = parfeval(@parsave,0,[Deci.Folder.Preproc filesep Deci.SubjectList{subject_list}],data_eeg,'data','-v7.3');
    
end


end