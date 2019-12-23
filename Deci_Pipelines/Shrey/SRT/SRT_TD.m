function [trl, trialinfo] = SRT_Trial_Definition(cfg)

%% Read EEG Data

hdr     = ft_read_header(cfg.dataset);
event   = ft_read_event(cfg.dataset);

%% Divide data into sequence and random blocks

block_seq = 11;
block_rnd = 12;
block_end = 13;
seq_beg = sprintf('S %d',block_seq);
rnd_beg = sprintf('S %d',block_rnd);
block_end = sprintf('S %d',block_end);

values = {event.value};
idx_beg = []; idx_end = []; idx_seq_beg=[];
for i = 1:length(values)
    if any(strcmp(values{i}, {seq_beg rnd_beg}))
        idx_beg = [idx_beg i];
    end
end

for i = 1:length(values)
    if strcmp(values{i}, block_end)
        idx_end = [idx_end i];
    end
end

for i = 1:length(values)
    if (strcmp(values{i}, seq_beg))
        idx_seq_beg = [idx_seq_beg i];
    end
end

if length(idx_beg)~=length(idx_end)
    idx_beg = idx_beg(1:96:length(idx_beg));
    idx_seq_beg = idx_seq_beg(1:96:length(idx_seq_beg));
end

seq_block_indices = find(ismember(idx_beg,idx_seq_beg));

for i = 1:length(idx_beg)
    blocks{i} = event(idx_beg(i):idx_end(i));
    if any(ismember(seq_block_indices,i))
        blocktype(i) = 1;
    else
        blocktype(i) = 2;
    end
end

%% Divide blocks into trials

% Declarations
trial_start = 20; trial_end = 30; 
Trial_end   = sprintf('S %d',trial_end);
Trial_beg   = sprintf('S %d',trial_start);

seq_count = 0; rnd_count = 0;

unit_rep_basis_start = 1:8:96;
unit_rep_basis_end = 8:8:96;
for i = 1:length(unit_rep_basis_start)
    unit_rep_basis(i,:) = unit_rep_basis_start(i):unit_rep_basis_end(i);
end
unit_basis = repmat(1:8, 12,1);
rep_basis = repmat(linspace(1,12,12)',1,8);

fdbckCorr   = sprintf('S %d',42);
fdbckIncorr = sprintf('S %d',43);
fdbckNoResp = sprintf('S %d',41);
fdbck       = {fdbckCorr fdbckIncorr fdbckNoResp};

% Division into trials
for p = 1:length(blocks)
    trl_idx_beg = []; trl_idx_end = [];
    values = {blocks{p}.value};
    for i = 1:length(values)
        if (strcmp(values{i}, Trial_beg))
            trl_idx_beg = [trl_idx_beg i];
        end
        if (strcmp(values{i}, Trial_end))
            trl_idx_end = [trl_idx_end i];
        end
    end
    for i = 1:length(trl_idx_beg)
        trials{p,i} = blocks{p}(trl_idx_beg(i):trl_idx_end(i));             % Trial defined
        if blocktype(p) == 1
            if i == 1 
                seq_count = seq_count+1;
            end
            seqorrnd(p,i) = 1;                                              % Whether it is a sequence block trial or a random block trial is defined
            blocknum(p,i) = seq_count;                                      % Of all the sequence or of all the random blocks, which position was this block at - e.g., trial belongs to the third sequence block tested
        else
            if i==1
                rnd_count = rnd_count+1;
            end
            seqorrnd(p,i) = 2;
            blocknum(p,i) = rnd_count;
        end
        repnum(p,i) = rep_basis(i == unit_rep_basis);                       % Within a block, which repition of the sequence is it (of course, valid only in sequence blocks but also applicable in random blocks for comparison)
        unitnum(p,i) = unit_basis(i == unit_rep_basis);                     % Within a sequence, what is the unit number of that trial, e.g., the second unit of a sequence
        
        subvalues = {trials{p,i}.value};
        if any(ismember(fdbck,subvalues))
            responseType(p,i) = find(ismember(fdbck,subvalues));
        else
            responseType(p,i) = NaN;
        end
         clear subvalues
    end   
    clear values
end

% Concatenate
m = 1;
for w = 1:size(trials,1)
    for v = 1:size(trials,2)
        alltrials{m} = trials{w,v};
        allseqorrnd(m) = seqorrnd(w,v);
        allblocknum(m) = blocknum(w,v);
        allrepnum(m) = repnum(w,v);
        allunitnum(m) = unitnum(w,v);
        allresptype(m) = responseType(w,v);
        m = m+1;
    end
end
        

%% Make the trl matrix consisting of sample starts and ends based on prestim, poststim, and type of trigs

% trig        = cfg.trig;
 stimonset   = {'S 21' 'S 22' 'S 23' 'S 24'};
 response    = {'S 25' 'S 26' 'S 27' 'S 28' 'S 29'}; 
% 
% if  strcmp(trig,'StimOnset')
%     trig_event = stimonset;
%     cfg.trialdef.prestim    = 1.5;
%     cfg.trialdef.poststim   = 3.5;
% elseif strcmp(trig,'Response')
%     trig_event = response;
%     cfg.trialdef.prestim    = 3.5;
%     cfg.trialdef.poststim   = 2.5;
% else
%     error('Please check the time t=0 string entered.');
% end

% Define trial boundaries with respect to the trigger t=0
% pretrig = -round(cfg.trialdef.prestim * hdr.Fs);
% posttrig = round(cfg.trialdef.poststim * hdr.Fs);

rjct = [];
trig_values = [];
trig_samples = [];
for i = 1:length(alltrials)
    values = {alltrials{i}.value};
    samples = {alltrials{i}.sample};
    
    disp(i)
    trig_values(i,:) = [find(ismember(stimonset,values))+20 find(ismember(response,values))+24];
    trig_samples(i,:) = [samples{ismember(values,stimonset)} samples{ismember(values,response)}];

    if length(trig_samples(i,:)) <= 1
        trig_samples(i,2) = NaN;
        trig_values(i,2) = NaN;
        rjct = [rjct; i];
    end
    
    begsample(i) = trig_samples(i,1) -round(2 * hdr.Fs);
    endsample(i) = trig_samples(i,2) + round(3 * hdr.Fs);
    offset(i) = 0;
    trlnum(i) = i;
end

trl = [begsample' endsample' offset' trlnum' begsample'-trig_samples];
    
%     4th column: Sequence block (1) or a random block (2)
%     5th column: Which sequence or random block?
%     6th column: Which repetition of the sequence? (1:12)
%     7th column: Which unit of the sequence? (1:8)
%     8th column: Accurate/inaccurate/no response? (1/2/3)
trialinfo = [trig_values allseqorrnd'+10 allblocknum'+100 allrepnum'+200 allunitnum'+300 allresptype'+50];

trl(rjct,:)         = [];
trialinfo(rjct,:)   = [];

end



