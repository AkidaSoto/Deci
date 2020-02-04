
function Deci = Checkor(Deci)


warning('off', 'MATLAB:MKDIR:DirectoryExists');
Deci.Folder.Definition   = [Deci.Folder.Version filesep 'Definition'];
Deci.Folder.Preproc      = [Deci.Folder.Version filesep 'PreProc'];
Deci.Folder.Artifact     = [Deci.Folder.Version filesep 'Artifact'];

if strcmp(Deci.SubjectList,'all')
    Deci.SubjectList = cellfun(@(c) strsplit(c,'.'),CleanDir(Deci.Folder.Raw),'un',0);
    Deci.SubjectList = unique(cellfun(@(c) c{1},Deci.SubjectList,'un',0));
    
    
elseif strcmp(Deci.SubjectList,'gui')
    Deci.SubjectList = cellfun(@(c) strsplit(c,'.'),CleanDir(Deci.Folder.Raw),'un',0);
    Deci.SubjectList = unique(cellfun(@(c) c{1},Deci.SubjectList,'un',0));
    
    fakeUI = figure;
    fakeUI.UserData = Deci.SubjectList;
    fakeUI.Visible =  'off';
    dc_select_labels(fakeUI,[],Deci.SubjectList);
    waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
    Deci.SubjectList = fakeUI.UserData;
    close(fakeUI);
end
disp(['Running Deci for ' num2str(length(Deci.SubjectList)) ' subjects'])
display(' ');

Pronouns = {'He','She','They'};

Deci = Exist(Deci,'Person',[]);

Deci.Person = Exist(Deci.Person,'Type','Concise'); % Concise
Deci.Person.Pro = [Pronouns(Randi(3))];
end