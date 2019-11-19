function Artifactor(Deci)

disp('----------------------');
disp('Starting Artifactor');
tic;

for subject_list = 1:length(Deci.SubjectList)
    
    if Deci.Art.do
        
        
        %% load Data
        data = [];
        cfg = [];
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
        condinfo = data.condinfo;
        preart   = data.preart;
        
        %% Do ICA
        
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
        
        
        %% Artifact Reject
        cfg =[];
        cfg.method = 'summary';
        cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
        cfg.eog = Deci.Art.eog;
        cfg.keepchannel = 'yes';
        tcfg.toilim = [abs(nanmax(condinfo{1},[],2)/1000)+Deci.Art.crittoilim(1) abs(nanmin(condinfo{1},[],2)/1000)+Deci.Art.crittoilim(2)];
        cfg.channel = 'all';
        
        
        data_rej = ft_rejectvisual(cfg,ft_redefinetrial(tcfg,data));
        
        cfg = [];
        cfg.trials = data_rej.saminfo;
        
        condinfo{1} = condinfo{1}(logical(data_rej.saminfo),:);
        condinfo{2} = condinfo{2}(logical(data_rej.saminfo),:);
        if length(condinfo) > 2
            condinfo{3} = condinfo{3}(logical(data_rej.saminfo));
        end
        
        if Deci.Art.AddComponents
            cfg =[];
            cfg.viewmode = 'vertical';
            
            scfg.trials = condinfo{3};
            data_comp = ft_selectdata(scfg,data_all);
            
            tcfg.toilim = [abs(nanmax(condinfo{1},[],2)/1000)+Deci.Art.crittoilim(1) abs(nanmin(condinfo{1},[],2)/1000)+Deci.Art.crittoilim(2)];
            
            artf = ft_databrowser(cfg,ft_redefinetrial(tcfg,data_comp));
            
            datacomp_rej = ft_rejectartifact(artf,ft_redefinetrial(tcfg,data_comp));
            
            condinfo{1} = condinfo{1}(logical(datacomp_rej.saminfo),:);
            condinfo{2} = condinfo{2}(logical(datacomp_rej.saminfo),:);
            if length(condinfo) > 2
                condinfo{3} = condinfo{3}(logical(datacomp_rej.saminfo));
            end
        end
        
        data.condinfo = condinfo;
        data.preart = preart;
        
        mkdir([Deci.Folder.Artifact])
        save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data','-v7.3')
        data = rmfield(data,'trial');
        save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} '_info'],'data','-v7.3')
        
    else
        mkdir([Deci.Folder.Artifact])
        
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
        
        data = rmfield(rmfield(data,'unmixing'),'topolabel');
        save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data','-v7.3')
        
        
        data = [];
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
        data = rmfield(data,'trial');
        save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} '_info'],'data','-v7.3')
    end
    
    
    
end