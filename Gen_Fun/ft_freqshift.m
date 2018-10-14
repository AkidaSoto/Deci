function freq = ft_freqshift(cfg,freq)

%cfg.offset Nxtrl size for redefinition
%cfg.parameter for which parameter to shift (usually fourierspctrm)
%cfg.bsl [beg end] for saving precue baseline
%cfg.keeptrials for if you want to keep trial info for bsl


if ~isfield(cfg,'offset')
    error('needs offset')
end

if ~isfield(cfg,'parameter')
    error('needs parameter')
end


if size(cfg.offset,1) ~= size(freq.(cfg.parameter),1)
    error('offset size needs to match trl size')
end

center = find(freq.time == min(abs(freq.time)));
% 
% bsl = freq.time >= cfg.bsl(1) & freq.time <= cfg.bsl(2);
% 
% bsltime = mean(freq.(cfg.parameter)(:,:,:,bsl),4);
% 
% if strcmp(cfg.keeptrials,'yes')
% bsltime = squeeze(mean(bsltime,1));
% end

datasize = size(freq.(cfg.parameter));

ft_progress('init', 'text',     'Please wait...');
for trl = 1:size(freq.(cfg.parameter),1)
    ft_progress(trl/size(freq.(cfg.parameter),1), ['freqshift %d, of %d'], trl,size(freq.(cfg.parameter),1));
    
    offset = round(freq.time,4) - round(cfg.offset(trl),4);
    offset = round(offset,4) >= min(round(cfg.latency,4)) & round(offset,4) <= max(round(cfg.latency,4));
    
    freq.shift(trl,:,:,:) = freq.(cfg.parameter)(trl,:,:,offset);
end
ft_progress('close');
freq.(cfg.parameter) = freq.shift;
freq.time = cfg.latency;
freq = rmfield(freq,'shift');
end