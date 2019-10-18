function out =  Deci_RT(Deci,info,dat,params)

maxt = max(sum(ismember(dat.preart{2},Deci.Analysis.Conditions{info.Cond}),2));
info.alltrials = find(sum(ismember(dat.preart{2},Deci.Analysis.Conditions{info.Cond}),2) == maxt);

%ignore all locks that are missing
minamountofnans = min(mean(isnan(dat.preart{1}(info.alltrials,:)),2));
info.allnonnans = mean(isnan(dat.preart{1}(info.alltrials,:)),2) == minamountofnans;% & ~isnan(mean(condinfo{2},2));

RT_Trials = dat.preart{1}(info.alltrials(info.allnonnans),:);

RT = -[RT_Trials(:,params.Lock(2)) - RT_Trials(:,params.Lock(1))];

if Deci.Analysis.ApplyArtReject
    RT = RT(find(ismember(dat.preart{3}(info.alltrials(info.allnonnans)),dat.condinfo{3})));
end

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'RT' filesep Deci.SubjectList{info.subject_list} ])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'RT' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.CondTitle{info.Cond}],'RT');


end