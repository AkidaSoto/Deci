function trl = subroutine_definetrial(cfg)

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
    del3OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 17'}))); %3
    imp2OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 18'}))); %4
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
end
    
trl = [begsample' endsample' offset' stimOnSample' imp1OnSample' del3OnSample' imp2OnSample' resp1Sample' resp2Sample' fd1OnSample' fd2OnSample' resp1Type' resp2Type' leftOrRight' Pb1OnSample' Pb2OnSample'];

if endsample(end) > samples{end}
   fprintf(['Removing last trial from ' cfg.dataset ' because epoch length greater than the number of samples.']); 
   trl(end,:) = [];
end

end
    
