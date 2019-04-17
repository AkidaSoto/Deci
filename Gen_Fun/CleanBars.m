function CleanBars(data)

 A =gca;
if length(find(size(data) ~= 1)) == 1
   
    
    bar(A,[zeros([1 length(data)]);data]);
    arrayfun(@(c) set(c,'YData',c.YData(2)),A.Children)

else
    bar(A,data);
end

end