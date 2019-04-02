function DefineTrialor(Deci)

disp('----------------------');
disp('Starting DefineTrialor');
tic;
for subject_list = 1:length(Deci.SubjectList)
     
    % Step 0: file check 
    if ~isempty(dir([Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.mat']))
        files_ending = {'.mat'};
    elseif   ~isempty(dir([Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.dat']))
        files_ending = {'.dat','.vmrk','.vhdr'};
    elseif   ~isempty(dir([Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.bdf']))
        files_ending = {'.bdf'};
    elseif ~isempty(dir([Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.eeg']))
        files_ending = {'.eeg','.vmrk','.vhdr'};
    elseif  ~isempty(dir([Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.cnt']))
        files_ending = {'.cnt'};
    elseif  ~isempty(dir([Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.rdf']))
        files_ending = {'.rdf','.vmrk','.vhdr'};
    elseif ~isempty(dir([Deci.Folder.Raw filesep Deci.SubjectList{subject_list} '.vhdr']))
        files_ending = {'.vhdr','.vmrk'};
    end
     
    for file_ending = 1:length(files_ending)
        if ~exist([Deci.Folder.Raw filesep Deci.SubjectList{subject_list} files_ending{file_ending}],'file')
            error([Deci.SubjectList{subject_list} ' does not have files for ' files_ending{file_ending} ' in raw data folder']);
        end
    end
    
    cfg = [];
    cfg.dataset = [Deci.Folder.Raw filesep Deci.SubjectList{subject_list} files_ending{1}];
    cfg.DT = Deci.DT;
    cfg.file_ending = files_ending{1};
    
    if strcmpi(Deci.DT.Type,'Manual') 
    cfg.trialfun = 'expfunor2';
    else
    cfg.trialfun = Deci.DT.Type;
    end
    
    cfg.Raw = Deci.Folder.Version;
    cfg.Subject = Deci.SubjectList{subject_list};
    evalc('cfg = ft_definetrial(cfg)'); % will return cfg.trl, the segmented data
    
    trllength = size(cfg.trl,1);
    disp(['Found ' num2str(trllength) ' trials for ' Deci.SubjectList{subject_list}]);
    
%     [~,i] = sort(cfg.trl(:,4));
%     cfg.trl = cfg.trl(i,:);
%     
%     if ~all(ismember([1:length(unique(cfg.DT.Locks))],unique(floor(cfg.trl(:,4)))))
%         error(['No trials were defined for condition(s) ' num2str(find(~ismember([1:length(unique(cfg.DT.Locks))],unique(floor(cfg.trl(:,4))))))]);
%     end
    
    mkdir([Deci.Folder.Definition]);
    disp('Saving Preprocessing Data');
    save([Deci.Folder.Definition filesep Deci.SubjectList{subject_list}],'cfg')

    
end

disp(['Finished DefineTrial at ' num2str(toc)]);
disp('----------------------');
end