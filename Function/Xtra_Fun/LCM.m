function [trl,trialinfo] = expfunor(cfg)

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
%%

for j = 1:length(startstopseg)
    
    value = {event(startstopseg(1,j):startstopseg(2,j)).value};
    value = cellfun(@str2num,value);
    sample = {event(startstopseg(1,j):startstopseg(2,j)).sample};
    
    % Add in lock for First Stim onset
    
    stim1 = find(ismember(value,[10 11]));
    value = [value(1:stim1) 9 value(stim1+1:end)];
    sample = [sample(1:stim1) sample(stim1) sample(stim1+1:end)];

    stim2 = find(ismember(value,[100 101 102]));
    value = [value(1:stim2) 99 value(stim2+1:end)];
    sample = [sample(1:stim2) sample(stim2) sample(stim2+1:end)];
    
    begsample = sample{ismember(value,cfg.DT.Locks(1))} + sstime(1)*hdr.Fs;
    
    if isempty(find(ismember(value,cfg.DT.Locks(end))))
        continue;
    else
        endsample = sample{ismember(value,cfg.DT.Locks(end))} + sstime(2)*hdr.Fs;
    end
    
    offsets = begsample - [sample{ismember(value,cfg.DT.Locks)}];
    
    trl(end + 1,1:4+length(cfg.DT.Locks)) = nan;
    trl(end ,1:4) = [begsample endsample 0 j];
    
    trl(end,logical([0 0 0 0 ismember(cfg.DT.Locks,value)])) = offsets;
    
    
    for i = 1:length(cfg.DT.Markers)
        if any(ismember([cfg.DT.Markers{i}],value))
            trialinfo(size(trl,1),i) = cfg.DT.Markers{i}(ismember([cfg.DT.Markers{i}],value));
        else
            trialinfo(size(trl,1),i) = nan;
        end
    end
end

k = 0;



