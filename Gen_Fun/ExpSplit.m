function ExpSplit(begcode,endcode,file)

    hdr   = ft_read_header(file);
    event = ft_read_event(file); %Mk<Marker number>=<Type>,<Description>,<Position in data points>,
    dat   = ft_read_data(file);
    
    event = StandardizeEventMarkers(event);
    
    
    
    
    ft_write_data([folder '_new' filesep  string{1}],dat,'header',hdr,'dataformat','brainvision_eeg','event',event2)



end