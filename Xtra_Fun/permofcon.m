function [mat]= permofcon(l)
%enter a cell with 4 matrices each filled with the number of markers per
%condition. Will modify at some point to take a variable number of
%arguments, but not today.

i = 0; mat={};  l1 = [l{1}]; l2=[l{2}]; l3=[l{3}]; l4=[l{4}]; 

for A = 1:length(l{1})
    for S = 1:length(l{2})
        for D = 1:length(l{3})
            for F = 1:length(l{4})
                i = i+1;
                mat(i) = {[l1(A) l2(S) l3(D) l4(F)]};%%Add brackets and it works but isnt a valid output for deci
            end
        end
    end
end
     
end