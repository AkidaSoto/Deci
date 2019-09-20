function happyfuntimes(dirr)

cdr = CleanDir(dirr);
for i = 1:length(cdr)
    files = CleanDir(strcat(dirr+cdr{i}));
    if exist(strcat(dirr+cdr{i}),'dir')
        files = files(~ismember(files,CleanDir(strcat(dirr+cdr{i}))));
    end
    
end

end