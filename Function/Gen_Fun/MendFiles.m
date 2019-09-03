Files = CleanDir(['*.vhdr']);
mended = [];
for f = 1:length(Files)
    
    mendme = find(cellfun(@(c) isequal(Files{f}(1:end-5),c(1:end-7)),Files));
    
    mainhdr = ft_read_header(Files{f});
    maindat = ft_read_data(Files{f});
    mainevent = ft_read_event(Files{f});
    
    if ~isempty(mendme) 
        
        for i = 1:length(mendme)
            xdat{i}  = ft_read_data(Files{mendme(i)});
            xevent{i} = ft_read_event(Files{mendme(i)});
            
            newsamps = -xevent{i}(2).sample + mainevent(end).sample + 1;
            
            for j = 2:length(xevent{i})
                xevent{i}(j).sample = xevent{i}(j).sample + newsamps;
            end
            
            mainevent = [mainevent xevent{i}];
            maindat   = [maindat xdat{i}];
            mainhdr.nSamples = size(maindat,2);
            mainhdr.orig.nSamples = size(maindat,2);
        end
        
        
        mended = [mended mendme];
    end
    
    if ~ismember(f,mended)
    mkdir('Mended');
    ft_write_data(['Mended' filesep Files{f}],maindat,'dataformat','brainvision_eeg','event',mainevent,'header',mainhdr);
    end
    
end