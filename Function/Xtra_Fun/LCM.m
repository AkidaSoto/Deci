function [trl,trialinfo] = expfunor(cfg)

startstop ={cfg.DT.Starts{:} cfg.DT.Ends{:}};
startstop = cellfun(@num2str,startstop,'un',0);
sstime = cfg.DT.Toi;

cfg.DT = Exist(cfg.DT,'Beha',[]);
cfg.DT.Beha = Exist(cfg.DT.Beha,'Acc',[]);
cfg.DT.Beha = Exist(cfg.DT.Beha,'RT',[]);

cfg.DT = Exist(cfg.DT,'Conditions',1:length(cfg.DT.Markers));

if ~strcmp(cfg.file_ending,'.mat')
    event = ft_read_event(cfg.dataset);
    hdr   = ft_read_header(cfg.dataset);
else
    events = [];
    load(cfg.dataset,'events');
    event = events;
    hdr.Fs = 1;
end

trl = [];

event =  StandardizeEventMarkers(event);

if rem(length(find(ismember({event.value},startstop)')),2) ~= 0
%     base = strsplit(cfg.headerfile,'.');
%     movefile(strcat(base{1},'.*'),'C:\Users\User\Desktop\Kyle\TheLostandDamned');
%     x = strsplit(base{1},'\');
%     error(x{end})
    error('invalid number of start and end pairings, check data!');
end

startstopseg = reshape(find(ismember({event.value},startstop)'),[2 length(find(ismember({event.value},startstop)'))/2]);

if ~isempty(cfg.DT.Block)
    bstartstop = {cfg.DT.Block.Start{:} cfg.DT.Block.End{:}};
    bstartstop = cellfun(@num2str,bstartstop,'un',0);
    
    if rem(length(find(ismember({event.value},bstartstop)')),2)
        blockvec = find(ismember({event.value},bstartstop)');
        blockvec = blockvec(1:end-1);
        
        bstartstopseg = reshape(blockvec,[2 length(blockvec)/2]);
        
    else
        bstartstopseg = reshape(find(ismember({event.value},bstartstop)'),[2 length(find(ismember({event.value},bstartstop)'))/2]);
    end
    
    if ~isempty(cfg.DT.Block.Markers)
        bmarkers = cellfun(@num2str,cfg.DT.Block.Markers,'un',0);
        
        bmarkers = find(ismember({event.value},bmarkers)');
        bmarkers = cellfun(@str2num,{event(bmarkers).value});
    else
        
    end
end

for j = 1:length(startstopseg) 
    
    value = {event(startstopseg(1,j):startstopseg(2,j)).value};
    value = cellfun(@str2num,value);
    sample = {event(startstopseg(1,j):startstopseg(2,j)).sample};
    
    if ~ismember(399,value) && ~isempty(value) && sum(isnan(value))~=0
        value = [value(1:find(ismember(value,211:218))) 399 value(find(ismember(value,211:218))+1:end)];
        sample = [sample(1:find(ismember(value,211:218))) sample(find(ismember(value,211:218))) sample(find(ismember(value,211:218))+1:end)];
    end
    
    if ~ismember(499,value) && ~isempty(value) && sum(isnan(value))~=0 %empty and NAN trial segs.  
        value = [value 499];
        sample = [sample sample{length(sample)}]
    end
    
    if all(ismember(cfg.DT.Locks,value))
        
        begsample = sample{1} + sstime(1)*hdr.Fs;
        endsample = sample{end} + sstime(2)*hdr.Fs;
        
        offsets = begsample - [sample{ismember(value,cfg.DT.Locks)}];
        
        trl(end + 1,:) = [begsample endsample 0 offsets];
         
        for i = 1:length(cfg.DT.Markers)
            
            if any(ismember([cfg.DT.Markers{i}],value))
                
                trialinfo(size(trl,1),i) = cfg.DT.Markers{i}(ismember([cfg.DT.Markers{i}],value));
               
                
                
                
            end
        end
        %% Marker Make 2.0
        %Currently only works with markers that differ by exactly 8 between
        %their minimum and maximum code. More or less will result in
        %incorrect orientation. TODO - Change subratrction pattern to length to
        %allow variable comparisons
        mkval1 = [211,221;211,151;141,151;141,211;141,221];
        code=[450,550,650,700,750];
        for gw = 1:size(mkval1,1) 
            id1 = findelement(max(mkval1(gw,:)),cfg.DT.Markers);
            id2 = findelement(min(mkval1(gw,:)),cfg.DT.Markers);
            difofmk = (max(mkval1(gw,:)) - min(mkval1(gw,:))); 
            if any(ismember([cfg.DT.Markers{id1}],value)) && any(ismember([cfg.DT.Markers{id2}],value))
                new = abs([cfg.DT.Markers{id1}(ismember([cfg.DT.Markers{id1}],value)) - cfg.DT.Markers{id2}(ismember([cfg.DT.Markers{id2}],value))]-difofmk);
                if new > 4
                    new = [8 - new];
                end
                trialinfo(size(trl,1),length(cfg.DT.Markers)+gw) = new + code(gw);
            end
        end
        
        
%         %dif in primer orientations 2,3 - original code or debug purposes
%         if any(ismember([cfg.DT.Markers{2}],value)) && any(ismember([cfg.DT.Markers{3}],value))
%             
%             new = abs(diff([cfg.DT.Markers{2}(ismember([cfg.DT.Markers{2}],value)) cfg.DT.Markers{3}(ismember([cfg.DT.Markers{3}],value))])-10);
%             
%             if new > 4
%                 new = [8 - new];
%             end
%             
%             trialinfo(size(trl,1),length(cfg.DT.Markers)+1) = new + 450;
%         end
%         %dif in distractor to sarch field target 2,6
      
        
        if ~isempty(cfg.DT.Block)
            
            if find([event(startstopseg(1,j)).sample] >  [event(bstartstopseg(1,:)).sample],1,'last')
                k = 0;
            end
            
            if ~isempty(cfg.DT.Block.Markers)
                trialinfo(size(trl,1),size(trialinfo,2)+1) = bmarkers(find([event(startstopseg(1,j)).sample] >  [event(bstartstopseg(1,:)).sample],1,'last'));
            else
                trialinfo(size(trl,1),size(trialinfo,2)+1) = find([event(startstopseg(1,j)).sample] >  [event(bstartstopseg(1,:)).sample],1,'last');
            end
            
        end
        
    end
    
end

    function tarelement = findelement(marker,cells)
        for tridex = 1:numel(cells)
            if ismember(marker,cells{tridex})
               tarelement = tridex;
            end 
        end
    end

end


