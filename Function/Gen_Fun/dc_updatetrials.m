function dc_updatetrials(dt_path,art_inpath,art_outpath)

% Current uses:
% updates events using locks as a references

%future uses:
% updates locks using events
% updates locks and events based on the loose assumption of trlnum
% switch labels between two arts

dt = CleanDir(dt_path);
art = CleanDir(art_inpath);

if ~all(ismember(dt,art))
   error('unequal dt and art lengths') 
end

for dts = 1:length(dt)

    cfg = [];
    load([dt_path filesep dt{dts}]);

    data = [];
    load([art_inpath filesep art{dts}]);
    
    if isfield(data,'condinfo')  %replacer starting 12/22, lets keep for ~4 months
        data.postart.locks = data.condinfo{1};
        data.postart.events = data.condinfo{2};
        data.postart.trlnum = data.condinfo{3};
        
        data.locks = data.preart{1};
        data.events = data.preart{2};
        data.trlnum = data.preart{3};
        
        
        data = rmfield(data,'condinfo');
        data = rmfield(data,'preart');
    end
    
    if ~isequaln(cfg.trl(:,4:end),data.locks)
       error('unequal locks beteen dt and art!') 
    end

    data.locks = cfg.trl(:,4:end);
    data.events = cfg.event;
    data.trlnum = cfg.trialnum;

    if ~isfolder(art_outpath)  
       mkdir(art_outpath) 
    end
    
    save([art_outpath filesep art{dts}],'data');
    data = rmfield(data,'trial');
    save([art_outpath filesep art{dts} '_info'],'data');
    
end