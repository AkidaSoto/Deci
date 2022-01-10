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
    
    for blktype = 1:length(cfg.DT.Block.Markers)
        bmarks = [cfg.DT.Block.Markers{blktype}];
        bmarks = arrayfun(@num2str,bmarks,'un',0);
        bstartstop = [find(ismember({event.value},bmarks))];
        bsstype{blktype} = cellfun(@(c) find(ismember(bmarks,c)),{event.value},'un',0);
        bsstype{blktype} = bsstype{blktype}(~cellfun(@isempty,bsstype{blktype}));
        bss{blktype} = bstartstop;
       
    end
end

for j = 1:length(startstopseg)
    
    value = {event(startstopseg(1,j):startstopseg(2,j)).value};
    value = cellfun(@str2num,value);
    sample = {event(startstopseg(1,j):startstopseg(2,j)).sample};
    
    
    if ~all(ismember(cfg.DT.Locks,value))
        
        cfg.DT = Exist(cfg.DT,'NanLocks',false);
        
        if ~cfg.DT.NanLocks
        continue;
        end
    end
    
    minlockindex = ismember(cfg.DT.Locks,value([sample{:}] == min([sample{ismember(value,cfg.DT.Locks)}])));
    maxlockindex = ismember(cfg.DT.Locks,value([sample{:}] == max([sample{ismember(value,cfg.DT.Locks)}])));
    
    begsample = sample{ismember(value,cfg.DT.Locks(minlockindex))} + sstime(1)*hdr.Fs;
    endsample = sample{ismember(value,cfg.DT.Locks(maxlockindex))} + sstime(2)*hdr.Fs;

    offsets = begsample - [sample{ismember(value,cfg.DT.Locks)}];
    
    trl(end + 1,1:4+length(cfg.DT.Locks)) = nan;
    trl(end ,1:4) = [begsample endsample 0 size(trl,1)];
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
        
        for blktypes = 1:length(cfg.DT.Block.Markers)
            trialinfo(size(trl,1),length(cfg.DT.Markers)+blktypes) = power(1,blktypes)*[-1*find([event(startstopseg(1,j)).sample] >  [event(bss{blktype}).sample],1,'last')];
        end
    else
    trialinfo(size(trl,1),length(cfg.DT.Markers)+1) = -1; 
    end
    
    
end

if ~isempty(cfg.DT.Block)
    if isfield(cfg.DT.Block,'Bisect')
        
        cfg.DT.Block = Exist(cfg.DT.Block,'Num',2);
        
        if cfg.DT.Block.Bisect
            bindex =  find(trialinfo(1,:) < 0);
            
            allbi = unique(trialinfo(:,bindex),'stable');
            
            for ab = allbi'
                
                onebi = find(ismember(trialinfo(:,bindex),ab));
                

                minnum = min(onebi);
                minpairs = min(length(onebi)/cfg.DT.Block.Num);

                onebi = cat(1,reshape(onebi(1:end-rem(length(onebi),cfg.DT.Block.Num)), [[length(onebi)-rem(length(onebi),cfg.DT.Block.Num)]/cfg.DT.Block.Num cfg.DT.Block.Num]), ...
                    [nan([1 cfg.DT.Block.Num-rem(length(onebi),cfg.DT.Block.Num)]) onebi(end-rem(length(onebi),cfg.DT.Block.Num)+1:end)']);

                
                for eachsect = 1:cfg.DT.Block.Num
                    
                    trialinfo(floor(minpairs*[eachsect-1])+minnum:floor(minpairs*[eachsect])+minnum-1,bindex) = trialinfo(floor(minpairs*[eachsect-1])+minnum:floor(minpairs*[eachsect])+minnum-1,bindex) + [-1/cfg.DT.Block.Num]*(eachsect-1);
                    
                end
                
            end
            
        end
    end
end

if isfield(cfg.DT,'Displace')
    
    cfg.DT.Displace = Exist(cfg.DT.Displace,'Num',0);
    
    if cfg.DT.Displace.Num ~= 0
        
        
        blkpos = find(trialinfo(1,:) < 0);
        blk = sort(unique(cfg.DT.Displace.BlockType(trialinfo(:,blkpos))),'descend');
        
        Displace = cfg.DT.Displace.Num;
        
        if isempty(cfg.DT.Displace.Static)
            cfg.DT.Displace.Static = nan;
        end
        
        for dura = 1:cfg.DT.Displace.Duration
            
            for stat = 1:length(cfg.DT.Displace.Static)

                if ~isnan(cfg.DT.Displace.Static)
                    statictrialinfo = trialinfo(logical(sum(ismember(trialinfo,cfg.DT.Displace.Static(stat)),2)),:);
                    statictrl = trl(logical(sum(ismember(trialinfo,cfg.DT.Displace.Static(stat)),2)),:);
                else
                    statictrialinfo = trialinfo;
                    statictrl = trl;
                end
                
                for dis = 1:length(blk)
                    
                    block{stat,dis,dura}  = statictrialinfo(find(cfg.DT.Displace.BlockType(statictrialinfo(:,blkpos)) == blk(dis)),:);
                    
                    if isempty(block{stat,dis,dura})
                        block{stat,dis,dura} = [];
                        continue;
                    end
                    
                    
                    shift{stat,dis,dura} = statictrl(find(cfg.DT.Displace.BlockType(statictrialinfo(:,blkpos)) == blk(dis)),:);
                    
                    if ~isempty(cfg.DT.Displace.Markers)
                        
                        for dmrk = 1:length(cfg.DT.Displace.Markers)
                            
                            dmrks{stat,dis,dura} = statictrialinfo(find(cfg.DT.Displace.BlockType(statictrialinfo(:,blkpos)) == blk(dis)),:);
                            
                            dmrks{stat,dis,dura} = dmrks{stat,dis,dura}(:,find(mean(ismember(dmrks{stat,dis,dura},cfg.DT.Displace.Markers{dmrk}))))+1000*dura;
                            
                            if cfg.DT.Displace.Num > 0
                                dmrks{stat,dis,dura} = [dmrks{stat,dis,dura}(1+[cfg.DT.Displace.Num+dura-1]:end);nan([cfg.DT.Displace.Num+dura-1 size(dmrks{stat,dis,dura},2)])];
                            else
                                dmrks{stat,dis,dura} = [nan([-cfg.DT.Displace.Num+dura-1 size(dmrks{stat,dis,dura},2)]); dmrks{stat,dis,dura}(1:end+[cfg.DT.Displace.Num+dura-1])];
                            end
                            
                            block{stat,dis,dura} = [block{stat,dis,dura} dmrks{stat,dis,dura}];
                            
                        end
                        
                        if cfg.DT.Displace.AddLocks
                            
                             for locks = 1:size(shift{stat,dis,dura}(:,5:end),2)
                                 
                                 if cfg.DT.Displace.Num > 0
                                   change =  shift{stat,dis,dura}(1:end-[cfg.DT.Displace.Num+dura-1],1) - shift{stat,dis,dura}(1+[cfg.DT.Displace.Num+dura-1]:end,1);
                                   shift{stat,dis,dura}(:,end+1) = [shift{stat,dis,dura}(1+[cfg.DT.Displace.Num+dura-1]:end,4+locks) + change; nan([cfg.DT.Displace.Num+dura-1 1])];
                                 else
                                   change =  shift{stat,dis,dura}(1:end+[cfg.DT.Displace.Num+dura-1],1) - shift{stat,dis,dura}(1-[cfg.DT.Displace.Num+dura-1]:end,1); 
                                   shift{stat,dis,dura}(:,end+1) = [nan([-cfg.DT.Displace.Num+dura-1 1]) ; shift{stat,dis,dura}(1:end+[cfg.DT.Displace.Num+dura-1],4+locks)];
                                   
                                   shift{stat,dis,dura}(1-[cfg.DT.Displace.Num+dura-1]:end,4+locks) = shift{stat,dis,dura}(1-[cfg.DT.Displace.Num+dura-1]:end,4+locks) + change;
                                 end 
                             end
                             
                             if cfg.DT.Displace.Num > 0
                             shift{stat,dis,dura}(1:end-[cfg.DT.Displace.Num+dura-1],2) = shift{stat,dis,dura}(1+[cfg.DT.Displace.Num+dura-1]:end,2);
                             else
                             shift{stat,dis,dura}(1-[cfg.DT.Displace.Num+dura-1]:end,1) = shift{stat,dis,dura}(1:end+[cfg.DT.Displace.Num+dura-1],1);    
                             end
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
            %format short g
        end
        
        
        
    end
    
end


