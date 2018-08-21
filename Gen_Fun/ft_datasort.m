function data = ft_datasort(cfg,data)


if ~isfield(cfg,'trl')
switch cfg.method
       
    case 'var'
        [~,trl] = sort(max(cell2mat(cellfun(@(c) var(c,[],2),data.trial,'un',0))),'descend');
    case 'std'
        [~,trl] = sort(max(cell2mat(cellfun(@(c) std(c,[],2),data.trial,'un',0))),'descend');
    case 'kurtosis'
        [~,trl] = sort(max(cell2mat(cellfun(@(c) kurtosis(c,[],2),data.trial,'un',0))),'descend');
end
else
     trl = cfg.trl;
end

if isfield(data,'sampleinfo')
    data.sampleinfo = data.sampleinfo(trl,:);
end

if isfield(data,'trialinfo')
     data.trialinfo = data.trialinfo(trl,:);
end

if isfield(data,'trial')
    data.trial = data.trial(1,trl);
end

if isfield(data,'time')
    data.time = data.time(1,trl);
end

data.oldtrl = trl;