function trl = expfunor(cfg)

startstop ={cfg.DT.Starts{:} cfg.DT.Ends{:}};
startstop = cellfun(@num2str,startstop,'un',0);
sstime = cfg.DT.Toi;

event = ft_read_event(cfg.dataset);
hdr   = ft_read_header(cfg.dataset);

trl = [];

event =  StandardizeEventMarkers(event);

if rem(length(find(ismember({event.value},startstop)')),2) ~= 0
    error('invalid number of start and end pairings, check data!')
end

startstopseg = reshape(find(ismember({event.value},startstop)'),[2 length(find(ismember({event.value},startstop)'))/2]);
for j = 1:length(startstopseg)
    
    value = {event(startstopseg(1,j):startstopseg(2,j)).value};
    value = cellfun(@str2num,value);
    sample = {event(startstopseg(1,j):startstopseg(2,j)).sample};
    
    eachcond = unique(cfg.DT.Conditions);
    for i = 1:length(eachcond)
        
        markers = cfg.DT.Markers(ismember(cfg.DT.Conditions,eachcond(i)));
        
        posvalue = cell2mat(cellfun(@(c) all(ismember(c(c >= 0),value)),markers,'un',0));
        negvalue = cell2mat(cellfun(@(c) all(~ismember(abs(c(c <= 0)),value)),markers,'un',0));
        
        condvalue = find(posvalue & negvalue);
        
        if length(condvalue) == 1
          
            condlock = [];
            condlock = sample{ismember(value,cfg.DT.Locks(i))};
            
            begsample = condlock + sstime(1)*hdr.Fs;
            endsample = condlock + sstime(2)*hdr.Fs;
            offset        = sstime(1)*hdr.Fs;
            
            trl(end + 1,:) = [begsample endsample offset  i+condvalue*.01];
            
            
        elseif length(condvalue) > 1
            error(['trial ' num2str(j) ' satisfies 2 or more trial definitions']);
        end
    end
end
