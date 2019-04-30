function AppendEEGs(Dir,datatype)

b = CleanDir(Dir);
b = unique(cellfun(@(d) d(1),cellfun(@(c) strsplit(c,'.'),b,'un',0)));

a = cellfun(@(c) isequal('2',c{end}(1)),cellfun(@(c) strsplit(c,'_'),b,'UniformOutput',false));

c = b(~a);
d = b(a);

d2 = cellfun(@(e) e([1:end-2]),d,'UniformOutput',false);

[commons,notcommons] = intersect(d2,c);

mkdir([Dir  '_new']);


for j = 1:length(commons)
    
    
    hdr   = ft_read_header([Dir filesep commons{j} '.' datatype]);
    event1 = ft_read_event([Dir filesep commons{j} '.' datatype]); 
    dat1   = ft_read_data([Dir filesep commons{j} '.' datatype]);
    
    event2 = ft_read_event([Dir filesep commons{j} '.' datatype]); 
    dat2   = ft_read_data([Dir filesep commons{j} '.' datatype]);
    
    event = [event1 arrayfun(@(c) setfield(c,'sample',c.sample + event1(end).sample),event2)];
    dat = cat(2,dat1,dat2);
    
    ft_write_data([Dir '_new' filesep commons{j}],dat,'header',hdr,'dataformat','brainvision_eeg','event',event)
    
end

movefile(notcommons,[Dir  '_new'])

end