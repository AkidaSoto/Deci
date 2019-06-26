function ExpSplit(begcode,endcode,file)

hdr   = ft_read_header(file);
event = ft_read_event(file); %Mk<Marker number>=<Type>,<Description>,<Position in data points>,
dat   = ft_read_data(file);

event = StandardizeEventMarkers(event);

if isnumeric(begcode)
    begcode = num2str(begcode);
end

if isnumeric(endcode)
    endcode = num2str(endcode);
end


begc = find(ismember({event.value},begcode));
endc = find(ismember({event.value},endcode));

if length(begc) ~= length(endc)
    error('unequal amounts of end and start codes');
end

data = [];
for i = 1:length(begc)
    
    sample = event(endc(i)).sample + hdr.Fs*10;
    
    if sample > size(dat,2)
       sample = size(dat,2); 
    end
    
    data = dat(:,1:sample);
    even = event(1:endc(i));
    
    ft_write_data([file(1:end-5) '_' num2str(i) '.vhdr'],data,'header',hdr,'dataformat','brainvision_eeg','event',even)
    
   
    
    if ~isempty(event)
        
        dat = dat(:,event(endc(i)).sample:end);
        sam =  event(endc(i)).sample;
        for j = 1:length(event)
        event(j).sample = [event(j).sample] - sam;
        end
        event = event(endc(i)+1:end);
           
        endc = endc - endc(i);
    end
end





end