%% Event Descriptors

% Trn_beg     = 71;         % beginning of training
% Trn_end     = 72;         % end of training
% 
% Exp_beg     = 1;          % Beginning of experiment
% Exp_end     = 2;          % End of experiment
% 
% Block_beg   = 99;         % Beginning of a new block
% Block_end   = 102;        % End of a block
% 
% Trial_beg   = 100;        % Beginning of a trial
% Trial_end   = 101;        % End of a trial
% 
% stimPresent = 11;         % two gratings presented in the trial
% delay1      = 12;         % gratings disappear, first delay period begins
% cueLeft     = 13;         % first delay period ends, cue is presented, if cue is on left this is the marker code...
% cueRight    = 14;         % ...or if cue is on right this is the marker code
% delay2      = 15;         % cue is removed, second delay begins
% impPresent  = 16;         % second delay ends, impulse stimulus is presented
% delay3      = 17;         % impulse is removed, third delay begins
% probePresent = 18;        % third delay ends, probe is presented
% respAntiCC  = 19;         % subjects says anticlockwise
% respCC      = 20;         % subject says clockwise
% noresp      = 21;         % subject doesn't give a response in the allotted time
% fdbckCorr   = 22;         % subject is presented 'Correct' feedback
% fdbckIncorr = 23;         % subject is presented 'Incorrect' feedback
% fdbckNoResp = 24;         % subject is presented feedback that 'No Response' was given
% respAcquired = 25;        % additional marker confirming response was acquired, to distinguish from no response trials
% 
% changeDir   = [61 62];    % whether actual change was counterclockwise or clockwise respectively
% changeMag   = 50+(1:length(des.change_mag)); % what was the change magnitude out of [3 7 12 18 25 33 42] respectively 


function [trl, trialinfo] = exp1_definetrial(cfg)

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
    del1OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 12'}))); %2
    cueOnSample(i)  = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 13'}) | strcmp(trials{i}.values,{'S 14'}))); %3
    cueType(i)      = find(ismember({'S 13','S 14'},trials{i}.values));               %4
    del2OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 15'}))); %5
    impOnSample(i)  = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 16'}))); %6
    del3OnSample(i) = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 17'}))); %7
    PbOnSample(i)   = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 18'}))); %8
%     respSample(i)   = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 25'}))); %9
    fdOnSample(i)   = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 22'}) | strcmp(trials{i}.values,{'S 23'}) | strcmp(trials{i}.values,{'S 24'}))); %3
    fdType(i)       = find(ismember({'S 22','S 23', 'S 24'},trials{i}.values));       %10
end

prestim    = 2;
poststim   = 3;

pretrig = -round(prestim * hdr.Fs);
posttrig = round(poststim * hdr.Fs);

for i = 1:length(trials)
    begsample(i) = stimOnSample(i) + pretrig;
    endsample(i) = fdOnSample(i) + posttrig;
    offset(i) = pretrig;
    zerostri(i) = 0;
end

trl = [begsample' endsample' zerostri' [1:length(trials)]' [begsample - stimOnSample]'  [begsample - cueOnSample]' [begsample - impOnSample]'];

trialinfo = [cueType'+73 fdType'+83];

if endsample(end) > samples{end}
   fprintf(['Removing last trial from ' cfg.dataset ' because epoch length greater than the number of samples.']); 
   trl(end,:) = [];
   trialinfo(end,:) = [];
end

end
    
