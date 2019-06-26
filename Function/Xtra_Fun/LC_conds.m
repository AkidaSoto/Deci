function [cond,math] = LC_conds(Conditions)

new = [];
for i = 1:length(Conditions)
    new(i) = abs(diff([Conditions{i}(2) Conditions{i}(3)])-10);
    
    if new(i) > 4
        new(i) = 8 - new(i);
    end
    
    if i ~= 1
        if new(i-1) == new(i)
            new(i) = new(i) + 5;
        end
    end
end

j = unique(new);

cond = [];
math = [];
for k = 1:length(j)
   
   cond{k} =  find(new == j(k));
   math{k} = ['[' strjoin(arrayfun(@(c) ['x' num2str(c)],find(new == j(k)),'UniformOutput',false),'+') ']/' num2str(length(find(new == j(k))))];
    
end

end