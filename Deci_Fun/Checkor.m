
function Deci = Checkor(Deci)

Deci.Folder.Definition   = [Deci.Folder.Version filesep 'Definition'];
Deci.Folder.Preproc      = [Deci.Folder.Version filesep 'PreProc'];
Deci.Folder.Artifact     = [Deci.Folder.Version filesep 'Artifact'];
Deci.Folder.Analysis     = [Deci.Folder.Version filesep 'Analysis'];
Deci.Folder.Plot         = [Deci.Folder.Version filesep 'Plot'];

if strcmp(Deci.SubjectList,'all')
    Deci.SubjectList = cellfun(@(c) strsplit(c,'.'),CleanDir(Deci.Folder.Raw),'un',0);
    Deci.SubjectList = unique(cellfun(@(c) c{1},Deci.SubjectList,'un',0));
elseif strcmp(Deci.SubjectList,'gui')
    Deci.SubjectList = cellfun(@(c) strsplit(c,'.'),CleanDir(Deci.Folder.Raw),'un',0);
    Deci.SubjectList = unique(cellfun(@(c) c{1},Deci.SubjectList,'un',0));
    
    fakeUI = figure;
    fakeUI.UserData = Deci.SubjectList;
    fakeUI.Visible =  'off';
    select_labels(fakeUI,[],Deci.SubjectList);
    waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
    Deci.SubjectList = fakeUI.UserData;
    close all
end

steps2 = {'Definition','PreProc','PreProc'};
steps1 = [2 3 4];

if Deci.Step == 1
   mkdir(Deci.Folder.Version) 
end

if any(ismember(steps1,[Deci.Step]))
    
    for subject_list = 1:length(Deci.SubjectList)

        if exist([Deci.Folder.Version filesep steps2{ismember(steps1,[Deci.Step])} filesep Deci.SubjectList{subject_list} '.mat']) ~= 2
            error([steps2{ismember(steps1,[Deci.Step])} ' step not found for ' Deci.SubjectList{subject_list}]);
        end
    end
    
    
end

copyfile(which('MainMenu'),[Deci.Folder.Version 'Parameter.text'])

end