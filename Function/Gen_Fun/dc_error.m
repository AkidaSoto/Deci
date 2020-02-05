function dc_error(Deci, str)
fid = fopen('adverb.txt', 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');
Starter = ['Deci ' C{1}{ismember({'Concise','Precise'},Deci.Person.Type)} ' says: '];

error([Starter str]);
end