function Markers = AllCTs(Markers,Variations)

X = find(Markers{1} == 'X');
A = find(Markers{1} == 'A');
B = find(Markers{1} == 'B');



if length([X A B]) ~= length(Variations)
    error('unequal amount of X and Variations')
end





for S = 1:length([X A B])
    
   
    
    mTemp = [];
    
    for V = 1:length(Variations{S})
        
        
        for M = 1:length(Markers)
            
           placement = regexp(Markers{M},'[A-Z]');
            
           mTemp{end+1} = [Markers{M}(1:placement-1)  Variations{S}{V} Markers{M}(placement+1:end)];
        end
        
    end
    
     Markers = mTemp;
    
end
Markers = cellfun(@(c) str2num(strjoin(strsplit(c,'_'),' ')),Markers,'un',0);
    
end