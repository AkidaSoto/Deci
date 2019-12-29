function [  ] = checksemicolons( functionnames )
%CHECKSEMICOLONS takes a cell argument with the filename(s) of an .m file to
% be checked, and will display to the command window a list of the lines
% that are missing semicolons
%
% It will not check subfunctions. I'm not that smart.
% E.g. > checksemicolons({'PDM9', 'PAAM9'}) 
%
% Harry Smith, University of Glasgow 2013
for f = 1:length(functionnames)
    functionname = functionnames{f};
    
    I = mlint(functionname, '-struct');
    j = 0;
    lineout = [];
    for i = 1:length(I)
        isErrormessage(i) = strcmp(I(i).message, 'Terminate statement with semicolon to suppress output (within a script).');
        
        if isErrormessage(i)
            j = j + 1;
            lineout(j) = I(i).line;
            
        end
        
    end
    
    
    if length(lineout) > 0
        disp(['Semicolons missing from line(s) ' num2str(lineout) ' in ' functionname])
    else
        disp(['No Semicolons missing from ' functionname])
    end
end
end