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
        data_rej = ft_rejectvisual(cfg,ft_redefinetrial(tcfg,data));
        
        cfg = [];
        cfg.trials = data_rej.saminfo;
        %data = ft_selectdata(cfg,data);
        
        condinfo{1} = condinfo{1}(logical(data_rej.saminfo),:);
        condinfo{2} = condinfo{2}(logical(data_rej.saminfo),:);
        if length(condinfo) > 2
            condinfo{3} = condinfo{3}(logical(data_rej.saminfo));
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