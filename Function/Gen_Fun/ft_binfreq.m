function freq = ft_binfreq(cfg,freq)

if isfield(cfg,'latency')
    tmpcfg = [];
    tmpcfg.latency = cfg.latency;
    
    freq = ft_selectdata(cfg,freq);
end

if isfield(cfg,'frequency')
    tmpcfg = [];
    tmpcfg.frequency = cfg.frequency;
    
    freq = ft_selectdata(cfg,freq);
end

if isfield(cfg,'freqbin')
    
    if ischar(cfg.freqbin)
        
        if strcmp(cfg.freqbin,'waves')
            cfg.freqbin = ...
                [1 4 ; ...
                4 8 ; ...
                8 13; ...
                13 30; ...
                30 60;];
            
        else
            error('incorrect use of cfg.freqbin');
        end
        
    elseif isnumeric(cfg.freqbin)
        
        if size(cfg.freqbin,2) ~= 2
            error('nx2 numeric matrix for manual binning')
        end
        
    elseif iscell(cfg.freqbin)
        
        freqbin = [];
        
        if ismember('theta',cfg.freqbin)
            freqbin(end+1,:) = [4 8];
        end
        if ismember('alpha',cfg.freqbin)
             freqbin(end+1,:) = [8 13];
        end
        
        if ismember('beta',cfg.freqbin)
             freqbin(end+1,:) = [13 30];
        end
        
        if ismember('gamma',cfg.freqbin)
             freqbin(end+1,:) = [30 60];
        end
        
        
        cfg.freqbin = freqbin;
        
    end
    
    tmpcfg = [];
    tmpcfg.avgoverfreq = 'yes';
    for f = 1:size(cfg.freqbin,1)
        tmpcfg.frequency = cfg.freqbin(f,:);
        
        freqdata{f} = rmfield(ft_selectdata(tmpcfg,freq),'cfg');
        
    end
    
    tmpcfg.parameter = 'fourierspctrm';
    tmpcfg.appenddim = 'freq';
    freq = ft_appendfreq(tmpcfg,freqdata{:});
    
else
    error('need freqbin field');
end




if isfield(cfg,'timebin')
    
    if isnumeric(cfg.timebin)
        
        if length(cfg.timebin) ~= 1
            if size(cfg.timebin,2) ~= 2
                error('nx2 numeric matrix for manual binning')
            end
        else
            cfg.timebin = discretize(freq.time,linspace(min(freq.time),max(freq.time),cfg.timebin));
            for k = 1:length(unique(cfg.timebin))
                timebin(k,:) = [freq.time(find(cfg.timebin == k,1,'first')) freq.time(find(cfg.timebin == k,1,'last'))];
            end
            cfg.timebin = timebin;
        end
        
        tmpcfg = [];
        tmpcfg.avgovertime = 'yes';
        for f = 1:size(cfg.timebin,1)
            tmpcfg.latency = cfg.timebin(f,:);
            
            timedata{f} = rmfield(ft_selectdata(tmpcfg,freq),'cfg');
            
        end
        
        tmpcfg.parameter = 'fourierspctrm';
        tmpcfg.appenddim = 'time';
        freq = ft_appendfreq(tmpcfg,timedata{:});
        
    else
        error('numerics only accepted for timebins');
        
    end
    
    
else
    error('need timebin field');
    
end









end