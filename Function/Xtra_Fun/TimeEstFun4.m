function [trl,trialinfo] = EEG_Polymerase(cfg)

%EEG_Polymerase will create
%trl: Relevant times, [TrialStart-Toi(1) TrialEnd+Toi(2) 0 trialnumber Lock1 Lock2 ....]
%trialinfo: Relevant Markers
%trialnumber will be transfered to it's own field during DefineTrialor

% Each Row represents 1 trial

startstop ={cfg.DT.Starts{:} cfg.DT.Ends{:}};
startstop = cellfun(@num2str,startstop,'un',0);
sstime = cfg.DT.Toi;

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
    
        bstartstop = [cfg.DT.Block.Markers{:}];
        bstartstop = arrayfun(@num2str,bstartstop,'un',0);
        bstartstop = [find(ismember({event.value},bstartstop))];
end

for j = 1:length(startstopseg)
    
    value = {event(startstopseg(1,j):startstopseg(2,j)).value};
    value = cellfun(@str2num,value);
    sample = {event(startstopseg(1,j):startstopseg(2,j)).sample};
    
    if ~all(ismember(cfg.DT.Locks,value))
        continue;
    else
        begsample = sample{ismember(value,cfg.DT.Locks(1))} + sstime(1)*hdr.Fs;
    end
    
    if isempty(find(ismember(value,cfg.DT.Locks(end))))
        continue;
    else
        endsample = sample{ismember(value,cfg.DT.Locks(end))} + sstime(2)*hdr.Fs;
    end
    
    offsets = begsample - [sample{ismember(value,cfg.DT.Locks)}];
    
    trl(end + 1,1:4+length(cfg.DT.Locks)) = nan;
    trl(end ,1:4) = [begsample endsample 0 size(trl,1)];
    trl(end,logical([0 0 0 0 ismember(cfg.DT.Locks,value)])) = offsets;
    
    
    for i = 1:length(cfg.DT.Markers)
        
        if any(ismember([cfg.DT.Markers{i}],value))
            
            
            trialinfo(size(trl,1),i) = cfg.DT.Markers{i}(ismember([cfg.DT.Markers{i}],value));
        else
            trialinfo(size(trl,1),i) = nan;
        end
    end
    
    
    trialinfo(size(trl,1),length(cfg.DT.Markers)+1) = nan;
    if ~isempty(cfg.DT.Block)   
        
        
            trialinfo(size(trl,1),length(cfg.DT.Markers)+1) = [-1*find([event(startstopseg(1,j)).sample] >  [event(bstartstop).sample],1,'last')];
    else
    trialinfo(size(trl,1),length(cfg.DT.Markers)+1) = -1; 
    end
    
    
end

   
end


