function [new_file] = filesyntax(filename)
macsyntax = '/';
pcsyntax = '\';

if ismac && ismember(pcsyntax,filename)
    new_file = strrep(filename,pcsyntax,macsyntax);
elseif ispc && ismember(macsyntax,filename)
    new_file = strrep(filename,macsyntax,pcsyntax);
elseif isunix && ismember(pcsyntax,filename)
     new_file = strrep(filename,pcsyntax,macsyntax);
else
    new_file = filename;
end
end