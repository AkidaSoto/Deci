function Artifactor(Deci)

disp('----------------------');
disp('Starting Artifactor');
tic;

for subject_list = 1:length(Deci.SubjectList)
    
    if Deci.Art.do
        
        data = [];
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
        
        
        condinfo = data.condinfo;
        preart   = data.preart;
        
        cfg =[];
        cfg.method = 'summary';
        cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
        cfg.eog = Deci.Art.eog;
        cfg.keepchannel = 'yes';
        tcfg.toilim = [abs(nanmax(condinfo{1},[],2)/1000)+Deci.Art.crittoilim(1) abs(nanmin(condinfo{1},[],2)/1000)+Deci.Art.crittoilim(2)];
        cfg.channel = 'all';
        
        if Deci.Art.AddComponents
            %non functional for now
            ccfg           = [];
            ccfg.numcomponent= 20;
            ccfg.unmixing  =data.unmixing;
            ccfg.topolabel = data.topolabel;
            ccfg.feedback = 'no';
            ccfg.demean     = 'no';
            data = rmfield(rmfield(data,'unmixing'),'topolabel');
            data_comp     = rmfield(ft_componentanalysis(ccfg, data),'cfg');
            cfg.components = data_comp.label;
        else
            data = rmfield(rmfield(data,'unmixing'),'topolabel');
        end
        
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
            data_comp = ft_selectdata(scfg,data_comp);
            
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
        copyfile([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat'],[Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} '.mat'])
        
        
        data = [];
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
        data = rmfield(data,'trial');
        save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} '_info'],'data','-v7.3')
    end
    
    
    
end