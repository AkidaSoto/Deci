function events = StandardizeEventMarkers(events)




cell_value = {events.value};

if any(cellfun(@isempty,cell_value))
    if  any(cellfun(@(c) isequal(c(1),'S'),cell_value(find(cellfun(@(c) ~isempty(c),cell_value)))))
        events(cellfun(@isempty,cell_value)) = arrayfun(@(c) setfield(c,'value','S  0'), events(cellfun(@isempty,cell_value)));
        cell_value(cellfun(@isempty,cell_value)) = arrayfun(@(c) 'S  0', events(cellfun(@isempty,cell_value)),'un',0);
        
        
        if ~isempty(events(~cellfun(@(c) isequal(c(1),'S'),cell_value)))
         
            if ismember('Buffer Overflow',cell_value)
                error('Buffer Overflow found, please reject subject')
            end
            
            
        events(~cellfun(@(c) isequal(c(1),'S'),cell_value)) = arrayfun(@(c) setfield(c,'value','S  0'), events(~cellfun(@(c) isequal(c(1),'S'),cell_value)))
        end
        
    else
        events(cellfun(@isempty,cell_value)) = arrayfun(@(c) setfield(c,'value','0'), events(cellfun(@isempty,cell_value)));
    end
end

cell_value = {events.value};
if isnumeric([cell_value{:}])
    events = arrayfun(@(c) setfield(c,'value', num2str(c.value)), events,'UniformOutput',0);
    events = [events{:}];
end

cell_value = {events.value};
if any(cellfun(@(c) isequal(c(1),'S'),cell_value))
    events = events(cellfun(@(c) isequal(c(1),'S'),cell_value));
    events = arrayfun(@(c) setfield(c,'value',strtrim(c.value(2:end))),  events,'UniformOutput',0);
    events = [events{:}];
end

if isfield(events,'type')
    events = events(~ismember({events.type},'Comment'));
end

[uniquesamples uniqueindex] = unique([[events.sample]' cellfun(@(c) str2num(c),{events.value})'],'rows');
events = events(uniqueindex);


end