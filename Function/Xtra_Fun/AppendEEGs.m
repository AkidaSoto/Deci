function AppendEEGs(Dir,datatype)

AllFiles = CleanDir(Dir);

if exist(strcat(Dir,'_new'),'dir')
 AllFiles = AllFiles(~ismember(AllFiles,CleanDir(strcat(Dir,'_new'))));
end

AllFiles = unique(cellfun(@(d) d(1),cellfun(@(c) strsplit(c,'.'),AllFiles,'un',0)));

IsCopy = cellfun(@(c) isstrprop(c{end},'digit'),cellfun(@(c) strsplit(c,'_'),AllFiles,'UniformOutput',false));

BaseFiles = AllFiles(~IsCopy);
CopyFiles = AllFiles(IsCopy);

CopysBaseName = cellfun(@(c) strjoin(c(1:end-1),'_'),cellfun(@(c) strsplit(c,'_'),AllFiles,'UniformOutput',false));


mkdir(strcat(Dir,'_new'));
for Each = 1:length(BaseFiles)
    
    Copies = CopyFiles(ismember(CopysBaseName,BaseFiles{Each}));
    
    hdr   = ft_read_header([Dir filesep BaseFiles{Each} '.' datatype]);
    event = ft_read_event([Dir filesep BaseFiles{Each} '.' datatype]); 
    dat   = ft_read_data([Dir filesep BaseFiles{Each} '.' datatype]);
    
    for cop = 1:length(Copies)
        
        event2 = ft_read_event([Dir filesep Copies{cop} '.' datatype]);
        dat2   = ft_read_data([Dir filesep Copies{cop} '.' datatype]);
        
        event = [event arrayfun(@(c) setfield(c,'sample',c.sample + event(end).sample),event2)];
        dat = cat(2,dat,dat2);
    end
   
    ft_write_data([Dir '_new' filesep BaseFiles{Each}],dat,'header',hdr,'dataformat','brainvision_eeg','event',event)
    
end

end