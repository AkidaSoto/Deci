function [trl,trialinfo] = EEG_Polymerase(cfg)

%EEG_Polymerase will create
%trl: Relevant times, [TrialStart-Toi(1) TrialEnd+Toi(2) 0 trialnumber Lock1 Lock2 ....]
%trialinfo: Relevant Markers
%trialnumber will be transfered to it's own field during DefineTrialor

% Each Row represents 1 trial

startstop ={cfg.DT.Starts{:} cfg.DT.Ends{:}};
startstop = cellfun(@num2str,startstop,'un',0);
sstime = cfg.DT.Toi;

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
    error('invalid number of start and end pairings, check data!')
end

startstopseg = reshape(find(ismember({event.value},startstop)'),[2 length(find(ismember({event.value},startstop)'))/2]);

if ~isempty(cfg.DT.Block)
    
        bstartstop = [cfg.DT.Block.Markers{:}];
        bstartstop = arrayfun(@num2str,bstartstop,'un',0);
        bstartstop = [find(ismember({event.value},bstartstop))];
end

for j = 1:length(startstopseg)
    
    value = {event(startstopseg(1,j):startstopseg(2,j)).value};
    value = cellfun(@str2num,value);
    sample = {event(startstopseg(1,j):startstopseg(2,j)).sample};
    
    begsample = sample{1} + sstime(1)*hdr.Fs;
    endsample = sample{end} + sstime(2)*hdr.Fs;
    
    offsets = begsample - [sample{ismember(value,cfg.DT.Locks)}];
    
    trl(end + 1,1:4+length(cfg.DT.Locks)) = nan;
    trl(end ,1:4) = [begsample endsample 0 j];
    
    trl(end,logical([0 0 0 0 ismember(cfg.DT.Locks,value)])) = offsets;
    
    
    for i = 1:length(cfg.DT.Markers)
        
        if any(ismember([cfg.DT.Markers{i}],value))
            
            
            trialinfo(size(trl,1),i) = cfg.DT.Markers{i}(ismember([cfg.DT.Markers{i}],value));
        else
            trialinfo(size(trl,1),i) = nan;
        end
    end
    
    
    trialinfo(size(trl,1),length(cfg.DT.Markers)+1) = nan;
    if ~isempty(cfg.DT.Block)   
        
        
            trialinfo(size(trl,1),length(cfg.DT.Markers)+1) = [-1*find([event(startstopseg(1,j)).sample] >  [event(bstartstop).sample],1,'last')];
    else
    trialinfo(size(trl,1),length(cfg.DT.Markers)+1) = -1; 
    end
    
    
end

if ~isempty(cfg.DT.Block)
    if isfield(cfg.DT.Block,'Bisect')
        if cfg.DT.Block.Bisect
            bindex =  find(trialinfo(1,:) < 0);
            
            allbi = unique(trialinfo(:,bindex),'stable');
            
            for ab = allbi'
                
                onebi = find(ismember(trialinfo(:,bindex),ab));
                
                onebi = {onebi(1:floor(length(onebi)/2)) onebi(ceil(length(onebi)/2):end)};

                for eachsect = 1:2
                    
                    trialinfo(onebi{eachsect},bindex) = trialinfo(onebi{eachsect},bindex) -.5*(eachsect-1);
                    
                end
                
            end
            
        end
    end
end



if isfield(cfg.DT,'Displace')
    
    cfg.DT.Displace = Exist(cfg.DT.Displace,'Num',0);
    
    if cfg.DT.Displace.Num ~= 0
        
        blk = sort(unique(ceil(trialinfo(:,end))),'descend');
        
        blkpos = find(trialinfo(1,:) < 0);
        
        
        Displace = cfg.DT.Displace.Num;
        
        for dura = 1:cfg.DT.Displace.Duration
            
            for stat = 1:length(cfg.DT.Displace.Static)
                
                statictrialinfo = trialinfo(logical(sum(ismember(trialinfo,cfg.DT.Displace.Static(stat)),2)),:);
                statictrl = trl(logical(sum(ismember(trialinfo,cfg.DT.Displace.Static(stat)),2)),:);
                
                
                for dis = 1:length(blk)
                    
                    block{stat,dis,dura}  = statictrialinfo(find(ceil(statictrialinfo(:,blkpos)) == blk(dis)),:);
                    
                    if isempty(block{stat,dis,dura})
                        block{stat,dis,dura} = [];
                        continue;
                    end
                    
                    
                    shift{stat,dis,dura} = statictrl(find(ceil(statictrialinfo(:,blkpos)) == blk(dis)),:);
                    
                    if ~isempty(cfg.DT.Displace.Markers)
                        
                        for dmrk = 1:length(cfg.DT.Displace.Markers)
                            
                            dmrks{stat,dis,dura} = statictrialinfo(find(ceil(statictrialinfo(:,blkpos)) == blk(dis)),:);
                            
                            dmrks{stat,dis,dura} = dmrks{stat,dis,dura}(:,find(mean(ismember(dmrks{stat,dis,dura},cfg.DT.Displace.Markers{dmrk}))))+1000*dura;
                            
                            if cfg.DT.Displace.Num > 0
                                dmrks{stat,dis,dura} = [dmrks{stat,dis,dura}(1+[cfg.DT.Displace.Num+dura-1]:end);nan([cfg.DT.Displace.Num+dura-1 size(dmrks{stat,dis,dura},2)])];
                            else
                                dmrks{stat,dis,dura} = [nan([cfg.DT.Displace.Num+dura-1 size(dmrks{stat,dis,dura},2)]); dmrks{stat,dis,dura}(1:end-[cfg.DT.Displace.Num+dura-1])];
                            end
                            
                            block{stat,dis,dura} = [block{stat,dis,dura} dmrks{stat,dis,dura}];
                            
                        end
                    end
                    
                    
                end
                
            end
            
            Displace = Displace + 1;
            
            trl = cat(1,shift{:,:,dura});
            
            [strl,I] = sort(trl(:,1));
            
            trl = trl(I,:);
            
            trialinfo = cat(1,block{:,:,dura});
            trialinfo = trialinfo(I,:);
            
        end
        
        
        
    end
    
end


