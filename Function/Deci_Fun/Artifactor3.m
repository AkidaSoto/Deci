function Artifactor2(Deci)


disp('----------------------');
disp('Starting Artifactor');
tic;

for subject_list = 1:length(Deci.SubjectList)
    
    data = [];
    load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
    
    condinfo = data.condinfo;
    preart   = condinfo;
    
    cfg =[];
    cfg.method = 'summary';
    cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
    cfg.eog = Deci.Art.eog.channel;
    cfg.keepchannel = 'yes';
    tcfg.toilim = [[abs(condinfo{1}(:,1))/1000]+Deci.Art.crittoilim(1) [abs(condinfo{1}(:,end))/1000]+Deci.Art.crittoilim(2)];
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
    
    data = rmfield(data,'cfg');
    
    mkdir([Deci.Folder.Artifact])
    save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data','-v7.3')
    data = rmfield(data,'trial');
    save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} '_info'],'data','-v7.3')
    
end