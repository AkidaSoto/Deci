
function tdata = ft_datashift2(cfg,data)


if ~isfield(cfg,'offset')
    error('needs offset')
end


if ~isfield(cfg,'shift')
    cfg.shift = 0;
end


if size(cfg.offset,1) ~= size(data.trialinfo,1)
    error('offset size needs to match trl size')
end

cfg.offset = [round(cfg.offset,4)/1000] - rem(round(cfg.offset,4)/1000,round(diff([data.time{1}(1) data.time{1}(2)]),4));

for k = 1:size(cfg.offset,1)
    time{k} = round(data.time{k},4) + round(cfg.offset(k),4);
    %data.trial{k}(:,ismember(round(data.time{k},4),round(time{k},4)));
end

    data.time = time;
    tcfg.toilim = cfg.toilim;
    tdata = ft_redefinetrial(tcfg,data);
    
end