function AppendEEGs(Dir,datatype,varargin)

if length(varargin) < 1
    base = 3;
else
    base = varargin{1};
end

AllFiles = CleanDir(Dir);

if exist(strcat(Dir,'_new'),'dir')
 AllFiles = AllFiles(~ismember(AllFiles,CleanDir(strcat(Dir,'_new'))));
end

AllFiles = unique(cellfun(@(d) d(1),cellfun(@(c) strsplit(c,'.'),AllFiles,'un',0)));

frags = cellfun(@(c) strsplit(c,'_'),AllFiles,'UniformOutput',false);

IsCopy = cellfun(@(c) length(c),frags);
IsCopy = IsCopy == max(IsCopy);

BaseFiles = AllFiles(~IsCopy);
CopyFiles = AllFiles(IsCopy);

BaseName  = cellfun(@(a) a(base),cellfun(@(c) strsplit(c,'_'),BaseFiles,'UniformOutput',false));
CopysBaseName = cellfun(@(a) a(base),cellfun(@(c) strsplit(c,'_'),CopyFiles,'UniformOutput',false));


mkdir(strcat(Dir,'_new'));
for Each = 1:length(BaseFiles)
    
    Copies = CopyFiles(ismember(CopysBaseName,BaseName{Each}));
    
    if ~isempty(Copies)
    hdr   = ft_read_header([Dir filesep BaseFiles{Each} '.' datatype]);
    event = ft_read_event([Dir filesep BaseFiles{Each} '.' datatype]); 
    dat   = ft_read_data([Dir filesep BaseFiles{Each} '.' datatype]);
    
    event = StandardizeEventMarkers(event);
    for cop = 1:length(Copies)
        
        event2 = ft_read_event([Dir filesep Copies{cop} '.' datatype]);
        dat2   = ft_read_data([Dir filesep Copies{cop} '.' datatype]);
        event2 = StandardizeEventMarkers(event2);

      event = [event arrayfun(@(c) setfield(c,'sample',c.sample + event(end).sample+1),event2)];
   %   event = [event arrayfun(@(c) setfield(c,'sample',c.sample + size(dat,2)+1open ),event2)];
      dat = cat(2,dat,dat2);
    end
   
    ft_write_data([Dir '_new' filesep BaseFiles{Each}],dat,'header',hdr,'dataformat','brainvision_eeg','event',event)
    else
       copyfile( [Dir filesep BaseFiles{Each} '.' datatype], [Dir '_new' filesep BaseFiles{Each} '.' datatype]);
       copyfile( [Dir filesep BaseFiles{Each} '.vhdr' ], [Dir '_new' filesep BaseFiles{Each} '.vhdr']);
       copyfile( [Dir filesep BaseFiles{Each} '.vmrk'], [Dir '_new' filesep BaseFiles{Each} '.vmrk']);
    end
end

end