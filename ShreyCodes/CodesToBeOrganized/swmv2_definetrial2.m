function [trl, trialinfo] = swmv2_definetrial(cfg)

hdr = ft_read_header([cfg.dataset]);
events = ft_read_event([cfg.dataset]);

values = {events.value};
samples = {events.sample};
t_starts = find(strcmp(values,'S100'));
t_ends = find(strcmp(values,'S101'));
for i = 1:length(t_starts)
    trials{i}.values = values(t_starts(i):t_ends(i));
    trials{i}.samples = samples(t_starts(i):t_ends(i));
end

for i = 1:length(trials)
    stimOnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 11'}))); %1
    imp1OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 13'}))); %2
    imp2OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 18'}))); %4
    del1OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 12'}))); %3
    del2OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 14'}))); %3
    del3OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 17'}))); %3
    del4OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 19'}))); %3
    resp1Sample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 71'})));  %5
    resp2Sample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 81'})));  %6
    fd1OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 91'})));  %7
    fd2OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 95'})));  %8
    resp1Type(i) = find(ismember({'S 74','S 75'},trials{i}.values));                  %9
    resp2Type(i) = find(ismember({'S 84','S 85'},trials{i}.values));                  %10
    leftOrRight(i) = find(ismember({'S 51','S 52'},trials{i}.values));                %11
    Pb1OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 15'})));  %12
    Pb2OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 20'})));  %13
end

prestim    = 2;
poststim   = 3;

pretrig = -round(prestim * hdr.Fs);
posttrig = round(poststim * hdr.Fs);

for i = 1:length(trials)
    begsample(i) = stimOnSample(i) + pretrig;
    endsample(i) = fd2OnSample(i) + posttrig;
    offset(i) = pretrig;
    zerostri(i) = 0;
end

trl = [begsample' endsample' zerostri' [1:length(trials)]' [begsample - stimOnSample]' [begsample - imp1OnSample]' [begsample - imp2OnSample]' [begsample - del1OnSample]' [begsample - del2OnSample]' [begsample - del3OnSample]' [begsample - del4OnSample]' [begsample - fd1OnSample]' [begsample - fd2OnSample]' [begsample-resp1Sample]' [begsample-resp2Sample]' [begsample - Pb1OnSample]' [begsample - Pb2OnSample]'];

trialinfo = [resp1Type'+73 resp2Type'+84 leftOrRight'+50];

%stimOnSample' imp1OnSample' del3OnSample' imp2OnSample' resp1Sample' resp2Sample' fd1OnSample' fd2OnSample' resp1Type' resp2Type' leftOrRight' Pb1OnSample' Pb2OnSample'];

trl = trl(11:end,:);
trialinfo = trialinfo(11:end,:);

if endsample(end) > samples{end}
   fprintf(['Removing last trial from ' cfg.dataset ' because epoch length greater than the number of samples.']); 
   trl(end,:) = [];
   trialinfo(end,:) = [];
end

end
    
