function dc_updatetrials(dt_path,art_inpath,art_outpath)

% Current uses:
% updates events using locks as a references

%future uses:
% updates locks using events
% updates locks and events based on the loose assumption of trlnum
% switch labels between two arts

dt = cellfun(@(c) strsplit(c,'.'),CleanDir(dt_path),'un',0);
dt = unique(cellfun(@(c) c{1},dt,'un',0));

fakeUI = figure;
fakeUI.UserData = dt;
fakeUI.Visible =  'off';
dc_select_labels(fakeUI,[],dt);
waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
dt = fakeUI.UserData;
close(fakeUI);


for dts = 1:length(dt)

    cfg = [];
    load([dt_path filesep dt{dts}]);

    data = [];
    load([art_inpath filesep dt{dts}]);
    
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

    data.postart.locks = data.locks(data.postart.trlnum,:);
    data.postart.events = data.events(data.postart.trlnum,:);
    

    if ~isfolder(art_outpath)  
       mkdir(art_outpath) 
    end
    
    info = rmfield(data,'trial');
    save([art_outpath filesep dt{dts}],'data','info','-v7.3');
end