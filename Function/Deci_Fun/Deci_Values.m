function out =  Deci_Values(Deci,info,dat,params)

Values = repmat(0,[1 length(dat.trial)]);

for Val = 1:length(params.Values)
maxt = length(find(cellfun(@(c) any(ismember(params.Markers{Val},c)), Deci.DT.Markers)));
Values(find(sum(ismember(dat.events,params.Markers{Val}),2) == maxt)) = params.Values(Val);
end

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Values' filesep Deci.SubjectList{info.subject_list} ])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'Values' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.CondTitle{info.Cond}],'Values');
end