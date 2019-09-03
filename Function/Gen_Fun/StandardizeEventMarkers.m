function events = StandardizeEventMarkers(events)




cell_value = {events.value};

if any(cellfun(@isempty,cell_value))
    events = events(~cellfun(@isempty,cell_value));
end

cell_value = {events.value};
if isnumeric([cell_value{:}])
    events = arrayfun(@(c) setfield(c,'value', num2str(c.value)), events,'UniformOutput',0);
    events = [events{:}];
end

cell_value = {events.value};
if cell_value{1}(1) == 'S'
    events = arrayfun(@(c) setfield(c,'value',strtrim(c.value(2:end))),  events,'UniformOutput',0);
    events = [events{:}];
end

if isfield(events,'type')
events = events(~ismember({events.type},'Comment'));
end

[uniquesamples uniqueindex] = unique([[events.sample]' cellfun(@(c) str2num(c),{events.value})'],'rows');
events = events(uniqueindex);


end