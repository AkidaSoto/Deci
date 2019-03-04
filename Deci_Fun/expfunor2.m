function [trl,trialinfo] = expfunor(cfg)

startstop ={cfg.DT.Starts{:} cfg.DT.Ends{:}};
startstop = cellfun(@num2str,startstop,'un',0);
sstime = cfg.DT.Toi;

cfg.DT = Exist(cfg.DT,'Beha',[]);
cfg.DT.Beha = Exist(cfg.DT.Beha,'Acc',[]);
cfg.DT.Beha = Exist(cfg.DT.Beha,'RT',[]);

cfg.DT = Exist(cfg.DT,'Conditions',1:length(cfg.DT.Markers));

if ~strcmp(cfg.file_ending,'.mat')
    event = ft_read_event(cfg.dataset);
    hdr   = ft_read_header(cfg.dataset);
else
    events = [];
    load(cfg.dataset,'events');
    event = events;
    hdr.Fs = 1;
end

trl = [];

event =  StandardizeEventMarkers(event);

if rem(length(find(ismember({event.value},startstop)')),2) ~= 0
    error('invalid number of start and end pairings, check data!')
end

startstopseg = reshape(find(ismember({event.value},startstop)'),[2 length(find(ismember({event.value},startstop)'))/2]);

if ~isempty(cfg.DT.Block)
    bstartstop = {cfg.DT.Block.Start{:} cfg.DT.Block.End{:}};
    bstartstop = cellfun(@num2str,bstartstop,'un',0);
    
    bstartstopseg = reshape(find(ismember({event.value},bstartstop)'),[2 length(find(ismember({event.value},bstartstop)'))/2]);

end

for j = 1:length(startstopseg)
    
    value = {event(startstopseg(1,j):startstopseg(2,j)).value};
    value = cellfun(@str2num,value);
    sample = {event(startstopseg(1,j):startstopseg(2,j)).sample};
    
    
    if all(ismember(cfg.DT.Locks,value))
        
        begsample = sample{1} + sstime(1)*hdr.Fs;
        endsample = sample{end} + sstime(2)*hdr.Fs;
        
        offsets = begsample - [sample{ismember(value,cfg.DT.Locks)}];
        
        trl(end + 1,:) = [begsample endsample 0 offsets];
        
        
        for i = 1:length(cfg.DT.Markers)
            
            if any(ismember([cfg.DT.Markers{i}],value))
                
                trialinfo(size(trl,1),i) = cfg.DT.Markers{i}(ismember([cfg.DT.Markers{i}],value));
                
            end
        end
        
    end
    
end

end


