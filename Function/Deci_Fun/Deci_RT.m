function out =  Deci_RT(Deci,info,dat,params)

RT_Trials = dat.locks;

RT = -[RT_Trials(:,params.Lock(1)) - RT_Trials(:,params.Lock(2))];

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'RT' filesep Deci.SubjectList{info.subject_list} ])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'RT' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.CondTitle{info.Cond}],'RT');
end