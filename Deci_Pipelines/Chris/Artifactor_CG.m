function Artifactor_CG(Deci)


for subject_list = 1:length(Deci.SubjectList)
    
    if Deci.Art.do
        

%% load Data
        data = [];
        cfg = [];
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
        
        if isfield(data,'condinfo')  
            data.postart.locks = data.condinfo{1};
            data.postart.events = data.condinfo{2};
            data.postart.trlnum = data.condinfo{3};
            
            data.locks = data.preart{1};
            data.events = data.preart{2};
            data.trlnum = data.preart{3};
            
            data = rmfield(data,'condinfo');
            data = rmfield(data,'preart');
        end
        
        locks = data.locks;
        events = data.events;
        trlnum = data.trlnum;
        postart.locks = data.postart.locks;
        postart.events = data.postart.events;
        postart.trlnum = data.postart.trlnum;
        
        
        
%         condinfo = data.condinfo;
%         preart   = data.preart;
%         ica = cfg;
        
%% Do ICA
        feedback = 'no';
        cfg.feedback = feedback;
        cfg.demean = 'no';
        data_ica = rmfield(ft_componentanalysis(cfg, data),'cfg');
        
        figure;
        cfg.component = [1:20];
        cfg.viewmode = 'component';
        
        clear cfg.method
        cfg.channel = 'all';
        
        cfg.component = [];
        
        comps = [data_ica.trial{:}];
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
            select_labels(fakeUI,[],sort(data_ica.label));
            fakeUI.Visible =  'off';
            ft_databrowser(cfg,data_ica);
            suptitle(Deci.SubjectList{subject_list});
            waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
            
            if isempty(fakeUI.UserData)
                cfg.component = [];
            else
                cfg.component = find(ismember(data_ica.label,fakeUI.UserData));
            end
            close(fakeUI)
            corr = [];

        else
            cfg.component = unique([component{:}]);
        end
        
        cfg.demean = 'yes';
        data = ft_rejectcomponent(cfg, data_ica);
        
%% Manual Trial Rejection
        
        if Deci.Art.Manual_Trial_Rejection
            cfg =[];
            cfg.method = 'trial';
            tcfg.toilim = [abs(nanmax(locks,[],2)/1000)+Deci.Art.crittoilim(1) abs(nanmin(locks,[],2)/1000)+Deci.Art.crittoilim(2)];

            cfg.viewmode = 'vertical';
            artf = ft_databrowser(cfg,ft_redefinetrial(tcfg,data));
            
            datacomp_rej = ft_rejectartifact(artf,ft_redefinetrial(tcfg,data));
            
            postart.locks = postart.locks(ismember(postart.trlnum,find(datacomp_rej.saminfo)),:);
            postart.events = postart.events(ismember(postart.trlnum,find(datacomp_rej.saminfo)),:);
            postart.trlnum = postart.trlnum(ismember(postart.trlnum,find(datacomp_rej.saminfo)));
                      
            cfg = [];
            cfg.trials = postart.trlnum;
            data_copy = ft_selectdata(cfg,data);
            
            display(' ')
            display('---Manual Trial Rejection Applied---')
            display(['Rejected ' num2str(length(find(~logical(datacomp_rej.saminfo)))) ' trials'])
            display(' ')
            pause(.05);
            
        else
            data_copy = data;
        end

%% Artifact Reject
        cfg =[];
        cfg.method = 'summary';
        cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
        cfg.eog = Deci.Art.eog;
        cfg.keepchannel = 'yes';
        tcfg.toilim = [abs(nanmax(locks,[],2)/1000)+Deci.Art.crittoilim(1) abs(nanmin(locks,[],2)/1000)+Deci.Art.crittoilim(2)];
        cfg.channel = 'all';
        
        
        data_rej = ft_rejectvisual(cfg,ft_redefinetrial(tcfg,data_copy));
        
        postart.locks = postart.locks(ismember(postart.trlnum,find(data_rej.saminfo)),:);
        postart.events = postart.events(ismember(postart.trlnum,find(data_rej.saminfo)),:);
        postart.trlnum = postart.trlnum(ismember(postart.trlnum,find(data_rej.saminfo)));
        
        display(' ')
        disp('---Trial Summary Rejection Applied---')
        disp(['Rejected ' num2str(length(find(~logical(data_rej.saminfo)))) ' trials'])
        display(' ')
        pause(.05);
        
        
        if Deci.Art.AddComponents
            cfg =[];
            cfg.viewmode = 'vertical';
            
            scfg.trials = condinfo{3};
            
            data_ica     = rmfield(ft_componentanalysis(ica, data),'cfg');
            data_comp = ft_selectdata(scfg,data_ica);
            
            tcfg.toilim = [abs(nanmax(condinfo{1},[],2)/1000)+Deci.Art.crittoilim(1) abs(nanmin(condinfo{1},[],2)/1000)+Deci.Art.crittoilim(2)];
            
            artf = ft_databrowser(cfg,ft_redefinetrial(tcfg,data_comp));
            
            datacomp_rej = ft_rejectartifact(artf,ft_redefinetrial(tcfg,data_comp));
            
            condinfo{1} = condinfo{1}(logical(datacomp_rej.saminfo),:);
            condinfo{2} = condinfo{2}(logical(datacomp_rej.saminfo),:);
            if length(condinfo) > 2
                condinfo{3} = condinfo{3}(logical(datacomp_rej.saminfo));
            end
        end
        
        data.locks = locks;
        data.events = events;
        data.trlnum = trlnum;
        data.postart = postart;
        
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

