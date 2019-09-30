function MakeMusic(File,Channel)



%% Load Data
data = [];
load(File);

bpm = 120;
bpsec = bpm/60;
secpb = 1/bpsec;

%% Segment Data

Chan = ismember(data.label,Channel);

for img = 1:length(data.trial)
 
    laststart = 1;
    
    for time = 1:size(data.trial{img},2)
    
       meanx = mean(data.trial{img}(Chan,laststart:time));
       xi =  data.trial{img}(Chan,time);
        
        if xi > 2*meanx || xi < 0
            
            seg{img}{end+1} = data.trial{img}(Chan,laststart:time);
            laststart = time;
        end
    end
end

%%

scfg.offset = data.condinfo{1};
scfg.toilim = [-secpb:secpb:20];
data = ft_datashift2(scfg,data);

tcfg.latency = [-secpb:secpb:20];
tcfg.keeptrials = 'yes';
tcfg.vartrllength = 2;
ERP = ft_timelockanalysis(tcfg,data);


fcfg.method        = 'wavelet';                                                  % Currently only uses 'wavelet' and 'hilbert'
fcfg.foi           = [6 12 35];                           % Frequency of Interest
fcfg.width         = 7 ;                                                         % Width
fcfg.gwidth        = 3;    
fcfg.pad       = 'nextpow2';

% Gwidth
fcfg.toi           = [-secpb:secpb:20];                                                  % Time Range to save
freq = ft_freqanalysis(fcfg,data);

end