function Plottor_CFC(Deci)

%% Deci Checks

if isfield(Deci.Plot.CFC,'freqhigh')
    if ischar(Deci.Plot.CFC.freqhigh)
        switch Deci.Plot.CFC.freqhigh
            case 'theta'
                Deci.Plot.CFC.freqhigh = [4 8];
            case 'beta'
                Deci.Plot.CFC.freqhigh = [12.5 30];
            case 'alpha'
                Deci.Plot.CFC.freqhigh = [8 12.5];
            case 'gamma'
                Deci.Plot.CFC.freqhigh = [30 50];
        end
    elseif isnumeric(Deci.Plot.CFC.freqhigh)
    else
        error(['cannot interrept freqhigh']);
    end
else
    error(['cannot interrept freqhigh']);
end

if isfield(Deci.Plot.CFC,'freqlow')
    if ischar(Deci.Plot.CFC.freqlow)
        switch Deci.Plot.CFC.freqlow
            case 'theta'
                Deci.Plot.CFC.freqlow = [4 8];
            case 'beta'
                Deci.Plot.CFC.freqlow = [12.5 30];
            case 'alpha'
                Deci.Plot.CFC.freqlow = [8 12.5];
            case 'gamma'
                Deci.Plot.CFC.freqlow = [30 50];
        end
    elseif isnumeric(Deci.Plot.CFC.freqlow)
    else
        error(['cannot interrept freqlow']);
    end
else
    error(['cannot interrept freqlow']);
end


if isequal(Deci.Plot.CFC.chanlow,'Reinhart-All')
    Deci.Plot.CFC.chanlow = [{'AF3'  } {'AF4'  } {'AF7'  } ...
        {'AF8'  } {'AFz'  } {'C1'   } {'C2'   } {'C3'   } {'C4'   } {'C5'   } ...
        {'C6'   } {'CP1'  } {'CP2'  } {'CP3'  } {'CP4'  } {'CP5'  } {'CP6'  } ...
        {'CPz'  } {'Cz'   } {'F1'   } {'F2'   } {'F3'   } {'F4'   } {'F5'   } ...
        {'F6'   } {'F7'   } {'F8'   } {'FC1'  } {'FC2'  } {'FC3'  } {'FC4'  } ...
        {'FC5'  } {'FC6'  } {'FCz'  } {'FT7'  } {'FT8'  } {'Fz'   } {'O1'   } ...
        {'O2'   } {'Oz'   } {'P1'   } {'P2'   } {'P3'   } {'P4'   } {'P5'   } ...
        {'P6'   } {'P7'   } {'P8'   } {'PO3'  } {'PO4'  } {'PO7'  } {'PO8'  } ...
        {'POz'  } {'Pz'   } {'T7'   } {'T8'   } {'TP10' } {'TP7'  } {'TP8'  } ...
        {'TP9'  } ] ;
end

if isequal(Deci.Plot.CFC.chanhigh,'Reinhart-All')
    Deci.Plot.CFC.chanhigh = [{'AF3'  } {'AF4'  } {'AF7'  } ...
        {'AF8'  } {'AFz'  } {'C1'   } {'C2'   } {'C3'   } {'C4'   } {'C5'   } ...
        {'C6'   } {'CP1'  } {'CP2'  } {'CP3'  } {'CP4'  } {'CP5'  } {'CP6'  } ...
        {'CPz'  } {'Cz'   } {'F1'   } {'F2'   } {'F3'   } {'F4'   } {'F5'   } ...
        {'F6'   } {'F7'   } {'F8'   } {'FC1'  } {'FC2'  } {'FC3'  } {'FC4'  } ...
        {'FC5'  } {'FC6'  } {'FCz'  } {'FT7'  } {'FT8'  } {'Fz'   } {'O1'   } ...
        {'O2'   } {'Oz'   } {'P1'   } {'P2'   } {'P3'   } {'P4'   } {'P5'   } ...
        {'P6'   } {'P7'   } {'P8'   } {'PO3'  } {'PO4'  } {'PO7'  } {'PO8'  } ...
        {'POz'  } {'Pz'   } {'T7'   } {'T8'   } {'TP10' } {'TP7'  } {'TP8'  } ...
        {'TP9'  } ] ;
end

%% Load

for  subject_list = 1:length(Deci.SubjectList)
    tic;
    for Channel = 1:length(Deci.Plot.CFC.chanhigh)
        
        display(['High Channel ' Deci.Plot.CFC.chanhigh{Channel} ' of ' num2str(length(Deci.Plot.CFC.chanhigh))])
        
        freq = [];
        load([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Lock filesep Deci.Plot.CFC.chanhigh{Channel}],'freq');
        
        foi = freq.freq >= round(Deci.Plot.CFC.freqhigh(1),4) & freq.freq <= round(Deci.Plot.CFC.freqhigh(2),4);
        toi = round(freq.time,4) >= Deci.Plot.CFC.latencyhigh(1) & round(freq.time,4) <= Deci.Plot.CFC.latencyhigh(2);
        
        HighChans{Channel} = freq;
        HighChans{Channel}.freq =  HighChans{Channel}.freq(foi);
        HighChans{Channel}.time =  HighChans{Channel}.time(toi);
        HighChans{Channel}.fourierspctrm  = HighChans{Channel}.fourierspctrm(:,:,foi,toi);
    end
    toc;
    
    acfg.parameter = 'fourierspctrm';
    acfg.appenddim = 'chan';
    HighChans = rmfield(ft_appendfreq(acfg,HighChans{:}),'cfg');
    
    tic;
    for Channel = 1:length(Deci.Plot.CFC.chanlow)
        
        display(['Low Channel ' Deci.Plot.CFC.chanlow{Channel} ' of ' num2str(length(Deci.Plot.CFC.chanlow))])
        
        freq = [];
        load([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Lock filesep Deci.Plot.CFC.chanlow{Channel}],'freq');
        
        foi = freq.freq >= round(Deci.Plot.CFC.freqlow(1),4) & freq.freq <= round(Deci.Plot.CFC.freqlow(2),4);
        toi = round(freq.time,4) >= Deci.Plot.CFC.latencylow(1) & round(freq.time,4) <= Deci.Plot.CFC.latencylow(2);
        
        LowChans{Channel} = freq;
        LowChans{Channel}.freq =  LowChans{Channel}.freq(foi);
        LowChans{Channel}.time =  LowChans{Channel}.time(toi);
        LowChans{Channel}.fourierspctrm  = LowChans{Channel}.fourierspctrm(:,:,foi,toi);
        
    end
    toc;
    
    acfg.parameter = 'fourierspctrm';
    acfg.appenddim = 'chan';
    LowChans = rmfield(ft_appendfreq(acfg,LowChans{:}),'cfg');
    
    for m = 1:length(Deci.Plot.CFC.methods)
        for Conditions = 1:length(Deci.Plot.Conditions)
            maxt = max(sum(ismember(freq.condinfo{2},Deci.Plot.Conditions{Conditions}),2));
            trl = sum(ismember(freq.condinfo{2},Deci.Plot.Conditions{Conditions}),2) == maxt;
            
            TrialCount{subject_list,Conditions} = length(find(trl));
            
            High = HighChans;
            High.fourierspctrm = HighChans.fourierspctrm(trl,:,:,:);
            Low = LowChans;
            Low.fourierspctrm = LowChans.fourierspctrm(trl,:,:,:);
            
            
            Deci.Plot.CFC.method = Deci.Plot.CFC.methods{m};
            CFCData{subject_list,Conditions,m} =  ft_singlecfc(Deci.Plot.CFC,Low,High);
            
        end
        
    end
    
end

end