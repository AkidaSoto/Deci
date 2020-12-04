
function SplitEEGs(Dir,datatype,code,time)

% Dir = Directory of files to split, New folder called [Dir '_new'] will be
% created

% datatype, 'eeg'

% Code to cut data by. preferably block start marker code
% buffer window to move away from beginning of the trial.

AllFiles = CleanDir(Dir);

if exist(strcat(Dir,'_new'),'dir')
 AllFiles = AllFiles(~ismember(AllFiles,CleanDir(strcat(Dir,'_new'))));
end

AllFiles = unique(cellfun(@(d) d(1),cellfun(@(c) strsplit(c,'.'),AllFiles,'un',0)));

mkdir(strcat(Dir,'_new'));
for Each = 1:length(AllFiles)
    
    disp('new subj starting')
    hdr   = ft_read_header([Dir filesep AllFiles{Each} '.' datatype]);
    event = ft_read_event([Dir filesep AllFiles{Each} '.' datatype]); 
    dat   = ft_read_data([Dir filesep AllFiles{Each} '.' datatype]);
   
    event = StandardizeEventMarkers(event);
    
    codetimes = [event(ismember({event.value},code)).sample];
    codeindex = find(ismember({event.value},code));
    cuttime = codetimes([end/2]+1*(time < 0))+time;
    cutindex = codeindex([end/2]+1*(time < 0));
   
    dat1 = dat(:,1:cuttime);
    event1 = event(1:cutindex);
    
    dat2 = dat(:,cuttime+1:end);
    event2 = event(cutindex+1:end);
    
    ft_write_data([Dir '_new' filesep  AllFiles{Each} '_1'],dat1,'header',hdr,'dataformat','brainvision_eeg','event',event1);
    ft_write_data([Dir '_new' filesep  AllFiles{Each} '_2'],dat2,'header',hdr,'dataformat','brainvision_eeg','event',event2);
    disp('new subj done')
    disp('----------------------------')
end

end