function [elec,pos] = CapTrakMake(file)

fid = fopen(file, 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
 
C=C{1,1};
 
%Find the line containing 'Impedance'

D = strfind(C, '</CapTrakElectrodeList>');
rows = find(~cellfun('isempty', D));

D = C(1:rows);

rows = strfind(D, '<Name>');
rows = D(~cellfun('isempty', rows));
 
%Read channelnames and impedance

elec = {};
for i = 1:length(rows)
     pat = '<[^>]*>';
    elec{i,1} = regexprep(rows{i}, pat, '');
end


rows = strfind(D, '<X>');
rows = D(~cellfun('isempty', rows));

X = {};
for i = 1:length(rows)
     pat = '<[^>]*>';
    X{i} = str2num(regexprep(rows{i}, pat, ''));
end


rows = strfind(D, '<Y>');
rows = D(~cellfun('isempty', rows));
Y = {};
for i = 1:length(rows)
     pat = '<[^>]*>';
    Y{i} =  str2num(regexprep(rows{i}, pat, ''));
end

rows = strfind(D, '<Z>');
rows = D(~cellfun('isempty', rows));
Z = {};
for i = 1:length(rows)
     pat = '<[^>]*>';
    Z{i} =  str2num(regexprep(rows{i}, pat, ''));
end

pos = cell2num(cat(2,X',Y',Z'));

end


