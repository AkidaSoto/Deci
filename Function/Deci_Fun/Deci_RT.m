function out =  Deci_RT(Deci,info,freq,params)

All_Trials = info.alltrials & info.allnonnans;

RT_Trials = freq.condinfo{1}(All_Trials,:);

RT = -[RT_Trials(:,params.Lock(2)) - RT_Trials(:,params.Lock(1))];

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'RT' filesep Deci.SubjectList{info.subject_list}  filesep filesep Deci.Analysis.LocksTitle{info.Lock}])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'RT' filesep Deci.SubjectList{info.subject_list}  filesep filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}],'RT');


end