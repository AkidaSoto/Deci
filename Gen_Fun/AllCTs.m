function Markers = AllCTs(Syntax,Variations)

if length(find(Syntax == 'X')) ~= length(Variations)
    error('unequal amount of X and Variations')
end

Marker = {};

for S = 1:length(find(Syntax == 'X'))
    
    tempS = Syntax;
    
    for V = 1:length(Variations)
        
        for inV = 1:length(Variations{V})
           placement = find(tempS == 'X',1,'first');
           tempS = [tempS(1:placement-1)  Variations{V}{inV} tempS(placement+1:end)];
        end
    end
    
    tempS = str2num(strsplit(tempS));
    Marker{end+1} = tempS;
end


end