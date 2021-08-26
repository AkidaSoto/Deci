clearvars;close all; clc;
macfg = [];
cfg.dataset     = '/Users/ReinhartLab/Desktop/Debra/t9.vhdr';
data_eeg        = ft_preprocessing(cfg);

hdr             = ft_read_header([cfg.dataset]);
events          = ft_read_event([cfg.dataset]);

values          = {events.value};
samples         = {events.sample};
t_starts        = find(strcmp(values,'S100'));
t_ends          = find(strcmp(values,'S101'));
for i = 1:length(t_starts)
    trials{i}.values    = values(t_starts(i):t_ends(i));
    trials{i}.samples   = samples(t_starts(i):t_ends(i));
end