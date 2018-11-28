%Read header file
fid = fopen('Prob_EEGTest_006PS.vhdr', 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
 
C=C{1,1};
 
%Find the line containing 'Impedance'
D = strfind(C, 'Impedance [kOhm]');
rows = find(~cellfun('isempty', D));
 
%Read channelnames and impedance
temp = {};

timer = textscan(char(C(rows)),'%s%s%s%s%s');
timer = timer{4};

for i = rows+1:length(C)
    
    temp = [temp; textscan(char(C(i)),'%s%d')];
end

%Remove colon from channel name
chantemp = temp(:,2);
for j = 1:length(chantemp)
    chan = char(chantemp{j});
    chantemp{j} = chan(1:end-1);
end
 
imp             = [];
imp.label       = chantemp;
imp.imp         = cell2mat(temp(:,2));
imp.time        = 1;
imp.dimord      = 'chan_time';
imp.timer = timer;

