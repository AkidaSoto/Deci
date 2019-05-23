function Artifactor(Deci)


disp('----------------------');
disp('Starting Artifactor');
tic;

for subject_list = 1:length(Deci.SubjectList)
    
    
    feedback =  'none';
    %% ICA
    
    if Deci.Art.do
        
        data = [];
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
        
        
        if ~isempty(find(cellfun(@(c) any(any(isnan(c))), data.trial) == 1)) && Deci.Art.RejectNans
            nantrials = find(cellfun(@(c) any(any(isnan(c))), data.trial) == 1);
            
            cfg = [];
            
            cfg.trials = logical(ones([size(data.trial)]));
            cfg.trials(nantrials) = false;
            
            data = ft_selectdata(cfg,data);
            warning(['Found trial(s) containing nan in rawdata for ' Deci.SubjectList{subject_list} '. Revise Data and then use .RejectNans']);
        elseif ~isempty(find(cellfun(@(c) any(any(isnan(c))), data.trial) == 1)) && ~Deci.Art.RejectNans
            
            error(['Found trial(s) containing nan in rawdata for ' Deci.SubjectList{subject_list} '. Revise Data and then use .RejectNans']);
            
        end
        condinfo = data.condinfo;
        preart   = condinfo;
        cfg = [];
        
        cfg.artfctdef.muscle = Deci.Art.muscle;
        
        [cfg, cfg.artfctdef.muscle.artifact] = ft_artifact_muscle(cfg, data);
        
        cfg.artfctdef.eog = Deci.Art.eog;
        
        [cfg, cfg.artfctdef.eog.artifact] = ft_artifact_eog(cfg, data);
        
        cfg.artfctdef.crittoilim = [[abs(condinfo{1}(:,1))/1000]+Deci.Art.crittoilim(1) [abs(condinfo{1}(:,end))/1000]+Deci.Art.crittoilim(2)];
        data_musc = ft_rejectartifact(cfg, data);
        
        tcfg.toilim = cfg.artfctdef.crittoilim;
        cfg =[];
        cfg.method = 'summary';
        cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
        
        
        data_musc = ft_rejectvisual(cfg,ft_redefinetrial(tcfg,data_musc));
        
        cfg.method  = 'runica';
        cfg.numcomponent= 20;
        cfg.feedback = feedback;
        data_musc = ft_componentanalysis(cfg, data_musc);
        
        cfg           = [];
        cfg.numcomponent= 20;
        cfg.unmixing  =data_musc.unmixing;
        cfg.topolabel = data_musc.topolabel;
        cfg.feedback = feedback;
        data_musc     = ft_componentanalysis(cfg, data);
        
        figure;
        cfg.component = [1:20];
        cfg.viewmode = 'component';
        
        cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
        
        cfg.comment   = 'no';
        ft_topoplotIC(cfg, data_musc);
        
        clear cfg.method
        
        cfg.channel = 'all';
        
        fakeUI = figure;
        select_labels(fakeUI,[],sort(data_musc.label));
        fakeUI.Visible =  'off';
        ft_databrowser(cfg,data_musc);
        suptitle(Deci.SubjectList{subject_list});
        waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
        
        if isempty(fakeUI.UserData)
            cfg.component = [];
        else
            cfg.component = find(ismember(data_musc.label,fakeUI.UserData));
        end
        close(fakeUI)
        
        %         cfg = [];
        %         cfg.method = 'summary';
        %         cfg.layout    = Deci.Layout.eye;
        %         data_musc = ft_rejectvisual(cfg,data_musc);
        %
        data_musc = ft_rejectcomponent(cfg, data_musc);
        
        if ~isempty(Deci.Art.PostICAbpf)
            bpcfg.bpfreq = Deci.Art.PostICAbpf;
            bpcfg.bpfilter = 'yes';
            data_musc = ft_preprocessing(bpcfg,data_musc);
        end
        
        cfg.artfctdef.muscle = Deci.Art.muscle;
        [cfg, cfg.artfctdef.muscle.artifact] = ft_artifact_muscle(cfg, data_musc);
        
        cfg.artfctdef.eog = Deci.Art.eog;
        [cfg, cfg.artfctdef.eog.artifact] = ft_artifact_eog(cfg, data_musc);
        
        cfg.artfctdef.crittoilim = [[abs(condinfo{1}(:,1))/1000]+Deci.Art.crittoilim(1) [abs(condinfo{1}(:,end))/1000]+Deci.Art.crittoilim(2)];
        data = ft_rejectartifact(cfg, data_musc);
        
        condinfo{1} = condinfo{1}(logical(data.saminfo),:);
        condinfo{2} = condinfo{2}(logical(data.saminfo),:);
        
        if length(condinfo) > 2
            condinfo{3} = condinfo{3}(logical(data.saminfo));
        end
        
        cfg =[];
        cfg.method = 'summary';
        cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
        cfg.channel = Deci.Art.eog.channel;
        cfg.keepchannel = 'yes';
        tcfg.toilim = [[abs(condinfo{1}(:,1))/1000]+Deci.Art.crittoilim(1) [abs(condinfo{1}(:,end))/1000]+Deci.Art.crittoilim(2)];
        data_rej = ft_rejectvisual(cfg,ft_redefinetrial(tcfg,data));
        
        cfg =[];
        cfg.trials = data_rej.saminfo;
        data = ft_selectdata(cfg,data);
        data.saminfo = data_rej.saminfo;
        condinfo{1} = condinfo{1}(logical(data_rej.saminfo),:);
        condinfo{2} = condinfo{2}(logical(data_rej.saminfo),:);
        
        if length(condinfo) > 2
            condinfo{3} = condinfo{3}(logical(data_rej.saminfo));
        end
        
        cfg = [];
        cfg.channel = 'all';
        data_rej = ft_rejectvisual(cfg,ft_redefinetrial(tcfg,data));
        
        cfg = [];
        cfg.trials = data_rej.saminfo;
        data = ft_selectdata(cfg,data);
        
        condinfo{1} = condinfo{1}(logical(data_rej.saminfo),:);
        condinfo{2} = condinfo{2}(logical(data_rej.saminfo),:);
        if length(condinfo) > 2
            condinfo{3} = condinfo{3}(logical(data_rej.saminfo));
        end
        data.condinfo = condinfo;
        data.preart = preart;
        mkdir([Deci.Folder.Artifact])
        save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data','-v7.3')
    else
        mkdir([Deci.Folder.Artifact])
        copyfile([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat'], [Deci.Folder.Artifact]);
    end
    
end
end

