function tdata = ft_datashift(cfg,data)


if ~isfield(cfg,'offset')
    error('needs offset')
end


if ~isfield(cfg,'shift')
    cfg.shift = 0;
end


if size(cfg.offset,1) ~= size(data.trialinfo,1)
    error('offset size needs to match trl size')
end

cfg.offset = round(cfg.offset,4) - rem(round(cfg.offset,4),round(diff([data.time{1}(1) data.time{1}(2)]),4));

for k = 1:size(cfg.offset,1)
    tcfg.begsample(k,1) = find(data.time{k} <= cfg.offset(k) + cfg.shift,1,'last');
    tcfg.endsample(k,1) = length(data.time{k});
    time{k} = cfg.shift:round(diff([data.time{1}(1) data.time{1}(2)]),4):[tcfg.endsample(k,1) - tcfg.begsample(k,1)]*round(diff([data.time{1}(1) data.time{1}(2)]),4) + cfg.shift;
end

    tdata = ft_redefinetrial(tcfg,data);
    tdata.time = time;
end