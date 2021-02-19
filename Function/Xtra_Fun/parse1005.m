
layout = 'C:\Users\User\Documents\GitHub\OurFieldTrip\Toolboxes\fieldtrip\template\layout\easycapM15.mat';
load(layout);

layout = lay;

fid = fopen('EEG1005.lay', 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');

pos = [];
for line = 1:numel(C{1})
    
    D = C{1}{line};
    D = textscan(D, '%s %s %s %s %s %s');
    
    pos.label{line} = D{6};
    pos.pos(line,:) = cellfun(@str2num,[D{2:3}]);
    pos.width(1) = lay.width(1);
    pos.height(1) = lay.height(1);
end

pos.outline = lay.outline;
pos.mask = lay.mask;

figure;
ft_plot_lay(pos);

% outline 4*4, x + .1, y -.1

lay = pos

