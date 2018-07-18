
function Deci = Checkor(Deci)

Deci.Folder.Definition   = [Deci.Folder.Version 'Definition'];
Deci.Folder.Preproc      = [Deci.Folder.Version 'PreProc'];
Deci.Folder.Artifact     = [Deci.Folder.Version 'Artifact'];
Deci.Folder.Analysis     = [Deci.Folder.Version 'Analysis'];
Deci.Folder.Plot         = [Deci.Folder.Version 'Plot'];

if strcmp(Deci.SubjectList,'all')
    Deci.SubjectList = cellfun(@(c) strsplit(c,'.'),CleanDir(Deci.Folder.Raw),'un',0);
    Deci.SubjectList = unique(cellfun(@(c) c{1},Deci.SubjectList,'un',0));
end

steps2 = {'Definition','PreProc','PreProc'};
steps1 = [2 3 4];

if Deci.Step == 1
   mkdir(Deci.Folder.Version) 
end

if any(ismember(steps1,[Deci.Step]))
    
    for subject_list = 1:length(Deci.SubjectList)

        if exist([Deci.Folder.Version steps2{ismember(steps1,[Deci.Step])} filesep Deci.SubjectList{subject_list} '.mat']) ~= 2
            error([steps2{ismember(steps1,[Deci.Step])} ' step not found for ' Deci.SubjectList{subject_list}]);
        end
    end
    
    
end

copyfile(which('MainMenu'),[Deci.Folder.Version 'Parameter.text'])

end