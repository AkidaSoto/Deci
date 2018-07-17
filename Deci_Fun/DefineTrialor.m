function DefineTrialor(Deci)

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
    end
     
    for file_ending = 1:length(files_ending)
        if ~exist([Deci.Folder.Raw filesep Deci.SubjectList{subject_list} files_ending{file_ending}],'file')
            error([Deci.SubjectList{subject_list} ' does not have files for ' files_ending{file_ending} ' in raw data folder']);
        end
    end
    
    cfg.dataset = [Deci.Folder.Raw filesep Deci.SubjectList{subject_list} files_ending{1}];
    cfg.DT = Deci.DT;
    cfg.trialfun = 'expfunor';
    cfg = ft_definetrial(cfg); % will return cfg.trl, the segmented data
    
    [~,i] = sort(cfg.trl(:,4));
    cfg.trl = cfg.trl(i,:);
    
    if ~all(ismember([1:length(cfg.DT.Markers)],unique(floor(cfg.trl(:,4)))))
        error(['No trials were defined for condition(s) ' num2str(find(~ismember([1:length(cfg.DT.Markers)],unique(floor(cfg.trl(:,4))))))]);
    end
    
    mkdir([Deci.Folder.Definition]);
    save([Deci.Folder.Definition filesep Deci.SubjectList{subject_list}],'cfg')

    
end


end