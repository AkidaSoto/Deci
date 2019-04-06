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
    
    if rem(length(find(ismember({event.value},bstartstop)')),2)
        blockvec = find(ismember({event.value},bstartstop)');
        blockvec = blockvec(1:end-1);
        
        bstartstopseg = reshape(blockvec,[2 length(blockvec)/2]);
        
    else
        bstartstopseg = reshape(find(ismember({event.value},bstartstop)'),[2 length(find(ismember({event.value},bstartstop)'))/2]);
    end
    
    if ~isempty(cfg.DT.Block.Markers)
        bmarkers = cellfun(@num2str,cfg.DT.Block.Markers,'un',0);
        
        bmarkers = find(ismember({event.value},bmarkers)');
        bmarkers = cellfun(@str2num,{event(bmarkers).value});
    else
        
    end
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
        
        if any(ismember([cfg.DT.Markers{3}],value)) && any(ismember([cfg.DT.Markers{4}],value))
            
            new = abs(diff([cfg.DT.Markers{3}(ismember([cfg.DT.Markers{3}],value)) cfg.DT.Markers{4}(ismember([cfg.DT.Markers{4}],value))])-10);
            
            if new > 4
                new = [8 - new];
            end
            
            trialinfo(size(trl,1),length(cfg.DT.Markers)+1) =new;
        end
        
        if ~isempty(cfg.DT.Block)
            
            if find([event(startstopseg(1,j)).sample] >  [event(bstartstopseg(1,:)).sample],1,'last')
                k = 0;
            end
            
            if ~isempty(cfg.DT.Block.Markers)
                trialinfo(size(trl,1),size(trialinfo,2)+1) = bmarkers(find([event(startstopseg(1,j)).sample] >  [event(bstartstopseg(1,:)).sample],1,'last'));
            else
                trialinfo(size(trl,1),size(trialinfo,2)+1) = find([event(startstopseg(1,j)).sample] >  [event(bstartstopseg(1,:)).sample],1,'last');
            end
            
        end
        
    end
    
end

end


