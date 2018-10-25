function Plottor(Deci)


Freq    = [];
ERP     = [];

for subject_list = 1:length(Deci.SubjectList)
    if Deci.Plot.Freq.Plot ||  ~isempty(Deci.Plot.PRP)
        if ~isdir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}])
            error(['Freq Analysis not found for '  Deci.SubjectList{subject_list}])
        end
        
        if isempty(Deci.Plot.Freq.Toi)
            Deci.Plot.Freq.Toi = [-inf inf];
            warning('Cannot Find Toi for plot, setting as [-inf,inf]');
        end
        
        if isempty(Deci.Plot.Freq.Type)
            error('Cannot Find Type for Freq plot');
        end
        
        if isempty(Deci.Plot.Freq.Foi)
            Deci.Plot.Freq.Foi = [-inf inf];
            warning('Cannot Find Foi for Freq plot, setting as [-inf inf]');
        end
        
        if isempty(Deci.Plot.Freq.Roi)
            Deci.Plot.Freq.Roi =  [];
            warning('Cannot Find Roi for Freq plot, setting as maxmin');
        end
        
        if isempty(Deci.Plot.Freq.BslType)
            Deci.Plot.Freq.BslType =  'maxmin';
            warning('Cannot Find BslType for Freq plot, setting as relative');
        end
        
        if isempty(Deci.Plot.Freq.Bsl)
            error('Cannot have Bsl empty for plot');
        end
        
        if isequal(Deci.Plot.Freq.Channel,'Reinhart-All')
            Deci.Plot.Freq.Channel = [{'AF3'  } {'AF4'  } {'AF7'  } ...
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
        
        
        Freq.Conditions     = CleanDir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}]);
        Freq.Channels       = CleanDir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list} filesep Freq.Conditions{1}]);
        Freq.Channels       = cellfun(@(c) c(1:end-4),Freq.Channels,'un',0);
        
        if ~ischar(Deci.Plot.Freq.Channel)
            if ~any(ismember(Deci.Plot.Freq.Channel,Freq.Channels))
                error('Freq Channel Parameter has additional channels not found in Data')
            end
            
            Freq.Channels = Deci.Plot.Freq.Channel(ismember(Deci.Plot.Freq.Channel,Freq.Channels));
        end
    end
    
    if ~isempty(Deci.Plot.ERP) || ~isempty(Deci.Plot.PRP)
        if ~isdir([Deci.Folder.Analysis filesep 'Volt_ERP' filesep Deci.SubjectList{subject_list}])
            error(['ERP Analysis not found for '  Deci.SubjectList{subject_list}])
        end
        
        ERP.Conditions     = CleanDir([Deci.Folder.Analysis filesep 'Volt_ERP' filesep Deci.SubjectList{subject_list}]);
        
        if isempty(Deci.Plot.ERP.Toi)
            Deci.Plot.ERP.Toi = [-inf inf];
            warning('Cannot Find Toi for plot, setting as [-inf,inf]');
        end
        
    end
    
    if ~isempty(Deci.Plot.PRP)
        Deci.Plot.GA = 0;
        
        if isempty(Deci.Plot.PRP.label)
            Deci.Plot.Freq.BslType =  0;
            warning('Cannot Find label for PRP plot, setting as 0');
        end
        
        if isempty(Deci.Plot.PRP.ScatterToi)
            Deci.Plot.Freq.ScatterToi =  [-inf inf];
            warning('Cannot Find ScatterToi for PRPplot, setting as  [-inf inf]');
        end
        
        
        
    end
    
    if ~isempty(Deci.Plot.CFC)
        
        
        if ~isdir([Deci.Folder.Analysis filesep 'Cross_' Deci.Plot.CFC.method])
            error(['CFC for ' Deci.Plot.CFC.method ' not found']);
        end
        
        CFC.Conditions     = CleanDir([Deci.Folder.Analysis filesep 'Cross_' Deci.Plot.CFC.method filesep Deci.SubjectList{subject_list}]);
        
    end
    
    
end

%% Load

if Deci.Plot.Freq.Plot || ~isempty(Deci.Plot.PRP)
    for  Condition =  1:length(Freq.Conditions)
        for  subject_list = 1:length(Deci.SubjectList)
            for Channel = 1:length(Freq.Channels)
                
                display(['Channel ' Freq.Channels{Channel} ' of ' num2str(length(Freq.Channels))])
                
                freq = [];
                load([Deci.Folder.Analysis filesep 'Freq_' Deci.Plot.Freq.Type filesep Deci.SubjectList{subject_list} filesep Freq.Conditions{Condition} filesep Freq.Channels{Channel}],'freq');
                
                foi = freq.freq >= round(Deci.Plot.Freq.Foi(1),4) & freq.freq <= round(Deci.Plot.Freq.Foi(2),4);
                toi = round(freq.time,4) >= Deci.Plot.Freq.Toi(1) & round(freq.time,4) <= Deci.Plot.Freq.Toi(2);
                
                Chans{Channel} = freq;
                Chans{Channel}.freq =  Chans{Channel}.freq(foi);
                Chans{Channel}.time =  Chans{Channel}.time(toi);
                Chans{Channel}.powspctrm  =Chans{Channel}.powspctrm(:,foi,toi);
                Chans{Channel}.label = Freq.Channels(Channel);
                
            end
            
            acfg.parameter = 'powspctrm';
            acfg.appenddim = 'chan';
            Subjects{subject_list} = rmfield(ft_appendfreq(acfg,Chans{:}),'cfg');
            clear Chans;
            
            if  exist([Deci.Folder.Version  filesep 'Redefine' filesep 'BSL' filesep Deci.SubjectList{subject_list} filesep num2str(Condition) '.mat']) ~= 2 || ~isequal(Deci.Plot.Freq.Bsl,[0 0])
                
                acfg.latency = Deci.Plot.Freq.Bsl;
                acfg.avgovertime = 'yes';
                bsl = ft_selectdata(acfg,freq);
            else
                
                bsl = [];
                load([Deci.Folder.Version  filesep 'Redefine' filesep 'BSL' filesep Deci.SubjectList{subject_list} filesep num2str(Condition) '.mat']);
            end
            
            bsl = repmat(bsl.powspctrm(ismember(bsl.label,Freq.Channels),foi,:),[1 1 size(Subjects{subject_list}.powspctrm ,3)]);
            
            switch Deci.Plot.Freq.BslType
                case 'absolute'
                    Subjects{subject_list}.powspctrm = Subjects{subject_list}.powspctrm - bsl;
                case 'relative'
                    Subjects{subject_list}.powspctrm= Subjects{subject_list}.powspctrm ./ bsl;
                case 'relchange'
                    Subjects{subject_list}.powspctrm = (Subjects{subject_list} - bsl) ./ bsl;
                case 'db'
                    Subjects{subject_list}.powspctrm = 10*log10(Subjects{subject_list}.powspctrm ./ bsl);
            end
            
        end
        
        if Deci.Plot.GA
            facfg.parameter =  'powspctrm';
            FreqData{Condition,1} = rmfield(ft_freqgrandaverage(facfg,Subjects{:}),'cfg');
            
        else
            FreqData(Condition,:) = Subjects(:);
        end
        clear Subjects;
    end
    
    
    if ~isempty(Deci.Plot.PRP)
        if strcmp(Deci.Plot.PRP.Dim,'trials')
            for  Condition =  1:length(Freq.Conditions)
                for  subject_list = 1:length(Deci.SubjectList)
                    for Channel = 1:length(Freq.Channels)
                        
                        display(['Channel ' Freq.Channels{Channel} ' of ' num2str(length(Freq.Channels))])
                        
                        freq = [];
                        load([Deci.Folder.Analysis filesep 'Four_' Deci.Plot.Freq.Type filesep Deci.SubjectList{subject_list} filesep Freq.Conditions{Condition} filesep Freq.Channels{Channel}],'freq');
                        
                        foi = freq.freq >= round(Deci.Plot.Freq.Foi(1),4) & freq.freq <= round(Deci.Plot.Freq.Foi(2),4);
                        toi = round(freq.time,4) >= Deci.Plot.Freq.Toi(1) & round(freq.time,4) <= Deci.Plot.Freq.Toi(2);
                        
                        Chans{Channel} = freq;
                        Chans{Channel}.freq =  Chans{Channel}.freq(foi);
                        Chans{Channel}.time =  Chans{Channel}.time(toi);
                        Chans{Channel}.powspctrm  =Chans{Channel}.powspctrm(:,:,foi,toi);
                        Chans{Channel}.label = Freq.Channels(Channel);
                        
                    end
                    
                    acfg.parameter = 'powspctrm';
                    acfg.appenddim = 'chan';
                    Subjects{subject_list} = rmfield(ft_appendfreq(acfg,Chans{:}),'cfg');
                    clear Chans;
                    
                    if  exist([Deci.Folder.Version  filesep 'Redefine' filesep 'BSL' filesep Deci.SubjectList{subject_list} filesep num2str(Condition) '.mat']) ~= 2 || ~isequal(Deci.Plot.Freq.Bsl,[0 0])
                        
                        acfg.latency = Deci.Plot.Freq.Bsl;
                        acfg.avgovertime = 'yes';
                        bsl = ft_selectdata(acfg,freq);
                    else
                        
                        bsl = [];
                        load([Deci.Folder.Version  filesep 'Redefine' filesep 'BSL' filesep Deci.SubjectList{subject_list} filesep num2str(Condition) '.mat']);
                    end
                    
                    bsl = repmat(bsl.powspctrm(:,ismember(bsl.label,Freq.Channels),foi,:),[1 1 1 size(Subjects{subject_list}.powspctrm ,4)]);
                    
                    switch Deci.Plot.Freq.BslType
                        case 'absolute'
                            Subjects{subject_list}.powspctrm = Subjects{subject_list}.powspctrm - bsl;
                        case 'relative'
                            Subjects{subject_list}.powspctrm= Subjects{subject_list}.powspctrm ./ bsl;
                        case 'relchange'
                            Subjects{subject_list}.powspctrm = (Subjects{subject_list} - bsl) ./ bsl;
                        case 'db'
                            Subjects{subject_list}.powspctrm = 10*log10(Subjects{subject_list}.powspctrm ./ bsl);
                    end
                    
                end
                
                if Deci.Plot.GA
                    facfg.parameter =  'powspctrm';
                    FreqData{Condition,1} = rmfield(ft_freqgrandaverage(facfg,Subjects{:}),'cfg');
                    
                else
                    FreqData(Condition,:) = Subjects(:);
                end
                
            end
        end
    end
    
    
    
    if Deci.Plot.Math.Type ~= 0
        
        if isequal(Deci.Plot.Math.Form,'Hemifield')
            
            for subj = 1:size(FreqData,2)
                HemiData = FreqData(:,subj);
                Hemis = length(HemiData);
                
                for Hem = 1:Hemis
                    HemiData{Hem} = hemifieldflip(HemiData{Hem});
                end
                
                Hemidata = [FreqData(:,subj);HemiData];
                
                for Hem = 1:Hemis
                    scfg.operation = ['x' num2str(Hem) '-x' num2str(Hem+Hemis)];
                    
                    scfg.parameter = 'powspctrm';
                    
                    MathData{Hem,subj} = ft_math(scfg,Hemidata{:});
                    
                    scfg.channel =  MathData{Hem,subj}.label(cell2mat(cellfun(@(c) any(rem(c(isstrprop(c,'digit')),2)),  MathData{Hem,subj}.label,'un',0)));
                    MathData{Hem,subj} = ft_selectdata(scfg, MathData{Hem,subj});
                end
                
                
                
            end
            
            FreqData = MathData;
            
        else
            
            for cond = 1:length(Deci.Plot.Math.Form)
                for subj = 1:size(FreqData,2)
                    scfg.parameter = 'powspctrm';
                    scfg.operation = Deci.Plot.Math.Form{cond};
                    MathData{subj} = ft_math(scfg,FreqData{:,subj});
                end
                
                FreqData(length(Freq.Conditions)+cond,:) = MathData;
            end
            
            if Deci.Plot.Math.Type == 1
                FreqData = FreqData(length(Freq.Conditions)+1:end,:);
            end
        end
        
        
        
    end
    
end


if ~isempty(Deci.Plot.Freq.Wires)
    for subj = 1:size(FreqData,2)
        for cond = 1:size(FreqData,1)
            WireData{cond,subj} =  ft_selectdata(Deci.Plot.Freq.Wires,FreqData{cond,subj});
        end
    end
end


if ~isempty(Deci.Plot.ERP)|| ~isempty(Deci.Plot.PRP)
    for  Condition =  1:length(ERP.Conditions)
        for  subject_list = 1:length(Deci.SubjectList)
            
            time = [];
            load([Deci.Folder.Analysis filesep 'Volt_ERP' filesep Deci.SubjectList{subject_list} filesep ERP.Conditions{Condition}],'time');
            
            toi = round(time.time,4) >= Deci.Plot.ERP.Toi(1) & round(time.time,4) <= Deci.Plot.ERP.Toi(2);
            time.time =  time.time(toi);
            time.avg  =time.avg(:,toi);
            
            if ~ischar(Deci.Plot.ERP.Channel)
                if ~any(ismember(Deci.Plot.ERP.Channel,time.label))
                    error('Channel Parameter has additional channels not found in Data')
                end
                
                tcfg.channel = Deci.Plot.ERP.Channel;
                time = ft_selectdata(tcfg,time);
            end
            
            Subjects{subject_list} = time;
            clear time;
        end
        
        if Deci.Plot.GA
            facfg.parameter =  'avg';
            TimeData{Condition,1} = rmfield(ft_timelockgrandaverage(facfg,Subjects{:}),'cfg');
            
        else
            TimeData(Condition,:) = Subjects(:);
        end
        clear Subjects;
        
    end
    
    
    if Deci.Plot.Math.Type ~= 0
        
        if isequal(Deci.Plot.Math.Form,'Hemifield')
            
            for subj = 1:size(TimeData,2)
                HemiData = TimeData(:,subj);
                Hemis = length(HemiData);
                
                for Hem = 1:Hemis
                    HemiData{Hem} = hemifieldflip(HemiData{Hem});
                end
                
                Hemidata = [TimeData(:,subj);HemiData];
                
                for Hem = 1:Hemis
                    scfg.operation = ['x' num2str(Hem) '-x' num2str(Hem+Hemis)];
                    
                    scfg.parameter = 'avg';
                    
                    MathData{Hem,subj} = ft_math(scfg,Hemidata{:});
                    
                    scfg.channel =  MathData{Hem,subj}.label(cell2mat(cellfun(@(c) any(rem(c(isstrprop(c,'digit')),2)),  MathData{Hem,subj}.label,'un',0)));
                    MathData{Hem,subj} = ft_selectdata(scfg, MathData{Hem,subj});
                end
                
            end
            
            TimeData = MathData;
            
        else
            
            for cond = 1:length(Deci.Plot.Math.Form)
                for subj = 1:size(TimeData,2)
                    scfg.parameter = 'avg';
                    scfg.operation = Deci.Plot.Math.Form{cond};
                    MathData{subj} = ft_math(scfg,TimeData{:,subj});
                end
                
                TimeData(length(ERP.Conditions)+cond,:) = MathData;
            end
            
            if Deci.Plot.Math.Type == 1
                TimeData = TimeData(length(ERP.Conditions)+1:end,:);
            end
        end
        
        
        
    end
    
    
    
    
end

if ~isempty(Deci.Plot.CFC)
    for  Condition =  1:length(CFC.Conditions)
        for  subject_list = 1:length(Deci.SubjectList)
            
            
            cfc = [];
            load([Deci.Folder.Analysis filesep 'Cross_' Deci.Plot.CFC.method filesep Deci.SubjectList{subject_list} filesep CFC.Conditions{Condition}],'cfc');
            
            Subjects{subject_list} = cfc;
        end
        
        if Deci.Plot.GA
            
            CFCData{Condition,1} = Subjects{1};
            CFCData{Condition,1}.crsspctrm = mean(cell2mat(permute(cellfun(@(c) c.crsspctrm,Subjects,'un',0),[1 3 4 2])),4);
            
        else
            CFCData(Condition,:) = Subjects(:);
        end
        
        
        
    end
end

%% Plot


if Deci.Plot.GA
    Deci.SubjectList = {'Group 1'};
end

if ~isempty(Deci.Folder.Plot)
    mkdir(Deci.Folder.Plot);
end

cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';

if Deci.Plot.Freq.Plot
    
    for subj = 1:size(FreqData,2)
        
        if Deci.Plot.Freq.Square
            square(subj) = figure;
        end
        
        if Deci.Plot.Freq.Topo
            topo(subj)  = figure;
        end
        
        if ~isempty(Deci.Plot.Freq.Wires)
            wire(subj)  = figure;
        end
        
        for cond = 1:size(FreqData,1)
            
            
            
            
            
            if Deci.Plot.Freq.Topo
                if length(Freq.Channels) ~= 1
                    
                    
                    set(0, 'CurrentFigure', topo(subj) )
                    topo(subj).Visible = 'on';
                    cirky(subj,cond)    =  subplot(size(FreqData,1),1,cond);
                    ft_topoplotER(cfg, FreqData{cond,subj});
                    
                    
                    %                     hold on
                    %                     a = get(gca,'Children');
                    %                     conty = a(arrayfun(@(c) strcmp(c.Type,'contour'),a));
                    %                     surfy  = a(arrayfun(@(c) strcmp(c.Type,'surface'),a));
                    %
                    %                     spot = find(arrayfun(@(c) strcmp(c.Type,'contour'),a));
                    %
                    %                     actmin = min(min(surfy.CData(~isnan(surfy.CData))));
                    %                     actmax = max(max(surfy.CData(~isnan(surfy.CData))));
                    %
                    %                     h = figure;
                    %                     [~,cont] = contourf(gca,surfy.CData,linspace(actmin,actmax,6));
                    %                     a(spot).LevelList = cont.LevelList;
                    %                     a(spot).TextList = cont.TextList;
                    %                     a(spot).Fill = 'on';
                    %                     delete(h);
                    
                    %                     axis xy
                    title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' Cond'  num2str(cond)]);
                    colorbar('vert');
                    map = colormap('jet'); %'hot' 'gray'
                    colormap(map);
                    
                else
                    close topo
                end
            end
            
            
            if strcmpi(Deci.Plot.Scale,'log')
                FreqData{cond,subj}.freq = log(FreqData{cond,subj}.freq);
                
                if ~isempty(Deci.Plot.Freq.Wires)
                    WireData{cond,subj}.freq = log(WireData{cond,subj}.freq);
                end
            end
            
            if Deci.Plot.Freq.Square
                set(0, 'CurrentFigure', square(subj) )
                square(subj).Visible = 'on';
                subby(subj,cond) = subplot(size(FreqData,1),1,cond);
                ft_singleplotTFR(cfg,FreqData{cond,subj})
                
                title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' Cond'  num2str(cond)]);
                colorbar('vert');
                map = colormap('jet'); %'hot' 'gray'
                colormap(map);
            end
            
            
            
            
        end
        
        if ~isempty(Deci.Plot.Freq.Wires)
            
            set(0, 'CurrentFigure', wire(subj) )
            curvy(subj)    =  axes;
            
            for cond = 1:size(WireData,1)
                
                if strcmpi(Deci.Plot.Freq.Wires.Collapse,'Freq')
                    plot(WireData{cond,subj}.time,squeeze(WireData{cond,subj}.powspctrm));
                    
                    xlabel('Time');
                elseif  strcmpi(Deci.Plot.Freq.Wires.Collapse,'Time')
                    
                    if strcmpi(Deci.Plot.Scale,'log')
                        WireData{cond,subj}.freq = exp(WireData{cond,subj}.freq);
                    end
                    
                    plot(WireData{cond,subj}.freq,squeeze(WireData{cond,subj}.powspctrm));
                    xlabel('Freq')
                end
                ylabel('Power')
                hold on
            end
            legend(curvy(subj),strsplit(num2str(1:size(FreqData,1))))
            title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type]);
            
            if ~isempty(Deci.Folder.Plot)
                saveas(wire(subj),[Deci.Folder.Plot filesep Deci.SubjectList{subj} '_wire'],Deci.Plot.Save.Format);
            end
            
        end
        
    end
    
    
    for subj = 1:size(FreqData,2)
        if length(Freq.Channels) ~= 1
            
            if Deci.Plot.Freq.Topo
                set(0, 'CurrentFigure', topo(subj) )
                topo(subj).Visible = 'on';
                
                for r = 1:length(cirky(:))
                    
                    if length(Deci.Plot.Freq.Roi) == 2 && isnumeric(Deci.Plot.Freq.Roi)
                        cirky(r).CLim = Deci.Plot.Freq.Roi;
                    elseif strcmp(Deci.Plot.Freq.Roi,'maxmin')
                        cirky(r).CLim = [min([cirky.CLim]) max([cirky.CLim])];
                    elseif strcmp(Deci.Plot.Freq.Roi,'maxabs')
                        cirky(r).CLim = [-1*max(abs([cirky.CLim])) max(abs([cirky.CLim]))];
                    end
                    
                end
                
                if ~isempty(Deci.Folder.Plot)
                    saveas(topo(subj),[Deci.Folder.Plot filesep Deci.SubjectList{subj} '_topo'],Deci.Plot.Save.Format);
                end
                
            end
            
            
        end
        
        if Deci.Plot.Freq.Square
            set(0, 'CurrentFigure', square(subj))
            square(subj).Visible = 'on';
            
            for r = 1:length(subby(:))
                
                if length(Deci.Plot.Freq.Roi) == 2 && isnumeric(Deci.Plot.Freq.Roi)
                    subby(r).CLim = Deci.Plot.Freq.Roi;
                elseif strcmp(Deci.Plot.Freq.Roi,'maxmin')
                    subby(r).CLim = [min([subby.CLim]) max([subby.CLim])];
                elseif strcmp(Deci.Plot.Freq.Roi,'maxabs')
                    subby(r).CLim = [-1*max(abs([subby.CLim])) max(abs([subby.CLim]))];
                end
            end
            
            if strcmpi(Deci.Plot.Scale,'log')
                for r = 1:length(subby(:))
                    subby(r).YTickLabels = exp(subby(r).YTick);
                end
            end
            
            if ~isempty(Deci.Folder.Plot)
                saveas(square(subj),[Deci.Folder.Plot filesep Deci.SubjectList{subj} '_square'],Deci.Plot.Save.Format);
            end
            
        end
        
        
    end
end

if ~isempty(Deci.Plot.ERP)
    if Deci.Plot.ERP.Plot
        for subj = 1:size(TimeData,2)
            
            topoERP(subj)  = figure;
            
            
            for cond = 1:size(TimeData,1)
                
                set(0, 'CurrentFigure', topoERP)
                tippy(cond,subj)    = subplot(size(TimeData,1),1,cond);
                ft_topoplotER(cfg, TimeData{cond,subj});
                
                
            end
            
            wireERP(subj)  = figure;
            ft_singleplotER(cfg, TimeData{:,subj});
            
        end
        
        for r = 1:length(tippy(:))
            tippy(r).CLim = [min([tippy.CLim]) max([tippy.CLim])];
        end
    end
end

if ~isempty(Deci.Plot.PRP)
    
    
    PRPPlot = figure;
    colors = reshape(hsv(size(FreqData,1)),[3 size(FreqData,1) 1]);
    
    for subj = 1:size(FreqData,2)
        for cond = 1:size(FreqData,1)
            
            %Collapse Freq and Channel
            FreqData{cond,subj}.avg = permute(nanmean(nanmean(FreqData{cond,subj}.powspctrm,2),1),[1 3 2]);
            FreqData{cond,subj} = rmfield(FreqData{cond,subj},'powspctrm');
            FreqData{cond,subj}.label = {'mix'};
            FreqData{cond,subj}.dimord = 'chan_time';
            Hz = regexprep(num2str(num2str(minmax(FreqData{cond,subj}.freq)),'%.3g '),' +','-');
            FreqData{cond,subj} = rmfield(FreqData{cond,subj},'freq');
            
            %Collapse Channel
            TimeData{cond,subj}.avg = mean(TimeData{cond,subj}.avg,1);
            TimeData{cond,subj}.label = {'mix'};
            TimeData{cond,subj}.dimord = 'chan_time';
        end
    end
    
    
    
    for cond = 1:size(FreqData,1)
        
        clear SF ST
        
        for  subj = 1:size(FreqData,2)
            
            for k = 1:length(TimeData{cond,subj}.time)
                SF(subj,k) = FreqData{cond,subj}.avg(:,k);
                ST(subj,k) = TimeData{cond,subj}.avg(:,k);
            end
            
            ax4 = subplot(2,2,4);
            ax4.XAxisLocation = 'origin';
            ax4.YAxisLocation = 'origin';
            hold(ax4,'on');
            scattoi = TimeData{cond,subj}.time >= Deci.Plot.PRP.ScatterToi(1) &  TimeData{cond,subj}.time <= Deci.Plot.PRP.ScatterToi(2);
            scat(cond,subj,1) = nanmean(FreqData{cond,subj}.avg(:,scattoi),2);
            scat(cond,subj,2) = nanmean(TimeData{cond,subj}.avg(:,scattoi),2);
            scattitle = regexprep(num2str(num2str(minmax(TimeData{cond,subj}.time(scattoi))),'%.3g '),' +','-');
            
            scatter(ax4,scat(cond,subj,1),scat(cond,subj,2),'MarkerFaceColor',colors(:,cond))
            
            if Deci.Plot.PRP.label
                t =  text(ax4, scat(cond,subj,1),scat(cond,subj,2),Deci.SubjectList{subj});
                set(t, 'Clipping', 'on','Interpreter', 'none');
            end
            
            
        end
        
        for k = 1:length(TimeData{cond,subj}.time)
            RHO(k) = corr(SF(:,k),ST(:,k));
        end
        
        ax3 = subplot(2,2,3);
        ax3.XAxisLocation = 'origin';
        ax3.YAxisLocation = 'origin';
        RHOim = plot(ax3,TimeData{cond,subj}.time,RHO,'Color',colors(:,cond));
        hold(ax3,'on');
        
        tcfg.parameter = 'avg';
        SFreq = rmfield(ft_timelockgrandaverage(tcfg,FreqData{cond,:}),'cfg');
        STime = rmfield(ft_timelockgrandaverage(tcfg,TimeData{cond,:}),'cfg');
        
        ax1 = subplot(2,2,1);
        ax1.XAxisLocation = 'origin';
        ax1.YAxisLocation = 'origin';
        plot(ax1,FreqData{cond,subj}.time,SFreq.avg,'Color',colors(:,cond));
        hold(ax1,'on');
        
        ax2 = subplot(2,2,2);
        ax2.XAxisLocation = 'origin';
        ax2.YAxisLocation = 'origin';
        plot(ax2,TimeData{cond,subj}.time,STime.avg,'Color',colors(:,cond));
        hold(ax2,'on');
        set(ax2,'YDir','reverse');
        
        mxb = polyfit(scat(cond,:,1),scat(cond,:,2),1);
        slope(cond) = mxb(1);
        plot(ax4,ax4.XLim,polyval(mxb,ax4.XLim),'Color',colors(:,cond));
    end
    
    xlabel(ax1, 'Time');
    ylabel(ax1, 'Total Power ');
    
    xlabel(ax2, 'Time');
    ylabel(ax2, 'Voltage (uV)');
    
    xlabel(ax3, 'Rho value');
    ylabel(ax3, 'P value');
    title(ax3,{['Between Trial Correlation with p and rho']} );
    
    
    xlabel(ax3,'Time(secs)');
    ylabel(ax3, 'Rho' );
    title(ax3,{['Correlation by time']} )
    
    xlabel(ax4,'Total Power');
    ylabel(ax4, 'ERP Power (uV)' );
    
    title(ax4,{['Correlation Total Power-ERP by subject'] ['at times ' scattitle]})
    
    
end

if ~isempty(Deci.Plot.CFC)
    for subj = 1:size(CFCData,2)
        
        if Deci.Plot.CFC.Square
            CFCsquare(subj) = figure;
        end
        
        if Deci.Plot.CFC.Topo
            CFCtopo(subj)  = figure;
        end
        
        if Deci.Plot.CFC.Hist
            CFChist(subj)  = figure;
        end
        
        for cond = 1:size(CFCData,1)
            
            
            if Deci.Plot.CFC.Topo
                set(0, 'CurrentFigure', CFCtopo(subj) )
                ctopo(subj,cond)    =  subplot(size(CFCData,1),1,cond);
                
                Cross.labelcmb = CombVec(CFCData{1}.labellow',CFCData{1}.labelhigh')';
                Cross.freq = CFCData{1}.freqlow;
                Cross.cohspctrm = CFCData{cond,subj}.crsspctrm(:,:,1);
                Cross.dimord = CFCData{1}.dimord;
                
                Crosscfg =[];
                Crosscfg.foi = [-inf inf];
                Crosscfg.layout = Deci.Layout.Noeye;
                title([Deci.SubjectList{subj} ' ' Deci.Plot.CFC.method ' Cond '  num2str(cond)],'Interpreter', 'none');
                %               ft_topoplotCC(Crosscfg,Cross);
                
                ctopo(subj,cond).UserData = {Cross,Crosscfg,CFCData{cond,subj}.crsspctrm};
            end
            
            if Deci.Plot.CFC.Hist
                set(0, 'CurrentFigure', CFChist(subj) )
                CFChist(subj).Visible = 'on';
                chist(subj,cond)    =  subplot(size(CFCData,1),1,cond);
                chist(subj,cond).UserData = CFCData{cond,subj}.crsspctrm;
                bar(reshape(CFCData{cond,subj}.crsspctrm,[size(CFCData{cond,subj}.crsspctrm,1)*size(CFCData{cond,subj}.crsspctrm,2) size(CFCData{cond,subj}.crsspctrm,3)]));
                
                title([Deci.SubjectList{subj} ' ' Deci.Plot.CFC.method ' Cond '  num2str(cond)],'Interpreter', 'none');
                
            end
            
            
            if Deci.Plot.CFC.Square
                
                set(0, 'CurrentFigure', CFCsquare(subj) )
                csquare(subj,cond)    =  subplot(size(CFCData,1),1,cond);
                
                xdat = CFCData{1}.labellow;
                ydat = CFCData{1}.labelhigh;
                PLVim = imagesc(1:length(xdat),1:length(ydat),CFCData{cond,subj}.crsspctrm(:,:,1));
                xticks(1:length(xdat));
                xticklabels(xdat);
                yticks(1:length(ydat))
                yticklabels(ydat);
                title([Deci.SubjectList{subj} ' ' Deci.Plot.CFC.method ' Cond '  num2str(cond)],'Interpreter', 'none');
                
                csquare(subj,cond).UserData = CFCData{cond,subj}.crsspctrm;
                colorbar;
            end
        end
        
        if Deci.Plot.CFC.Hist
            set(0, 'CurrentFigure', CFChist(subj));
            UpdateAxes(chist(subj,:),Deci.Plot.CFC.Roi,'Y',0)
            
            for HistLim = 1:length(chist(subj,:))
                
                chist(subj,HistLim).YLim(1) = chist(subj,HistLim).YLim(1) *.98;
                chist(subj,HistLim).YLim(2) = chist(subj,HistLim).YLim(2) *1.02;
                
                xticklabels(chist(subj,HistLim), cell2mat(CombVec(CFCData{cond,subj}.labellow',CFCData{cond,subj}.labelhigh')'));
                xtickangle(chist(subj,HistLim),-20);
                
                legend(chist(subj,HistLim),[repmat('FreqLow Time ',[size(CFCData{cond,subj}.timelow',1) 1]) num2str([CFCData{cond,subj}.timelow']) repmat(' - FreqHigh Time ',[size(CFCData{cond,subj}.timelow',1) 1]) num2str([CFCData{cond,subj}.timelow'])])
            end
            
            
            
        end
        
        
        if Deci.Plot.CFC.Topo
            
            set(0, 'CurrentFigure', CFCtopo(subj) )
            uicontrol('style','text','position',[225 75 100 25],'String','Time of Interest');
            
            slide = uicontrol('style','slider','position',[75 10 400 20],...
                'min',1,'max',size(CFCData{1,subj}.crsspctrm,3),'callback',{@ChangeDimTopo,ctopo(subj,:),Deci.Plot.CFC.Roi}, ...
                'value',1,'SliderStep',[1/size(CFCData{1,subj}.crsspctrm,3) 1/size(CFCData{1,subj}.crsspctrm,3)]);
            
            ChangeDimTopo(slide,[],ctopo(subj,:),Deci.Plot.CFC.Roi);
            
            
            for tick = 1:size(CFCData{1,subj}.crsspctrm,3)
                uicontrol('style','text','position',[75+[[[slide.Position(3)]/5]*[tick-1]]+20 55 40 25],'String',num2str(round(CFCData{1,subj}.timelow(tick),2)));
                uicontrol('style','text','position',[75+[[[slide.Position(3)]/5]*[tick-1]]+20 30 40 25],'String',num2str(round(CFCData{1,subj}.timehigh(tick),2)));
            end
            uicontrol('style','text','position',[45 55 60 25],'String','FreqLow');
            uicontrol('style','text','position',[45 30 60 25],'String','FreqHigh');
        end
        
        
        if Deci.Plot.CFC.Square
            
            set(0, 'CurrentFigure', CFCsquare(subj) )
            uicontrol('style','text','position',[225 75 100 25],'String','Time of Interest');
            
            slide = uicontrol('style','slider','position',[75 10 400 20],...
                'min',1,'max',size(CFCData{1,subj}.crsspctrm,3),'callback',{@ChangeDim,csquare(subj,:),Deci.Plot.CFC.Roi,'C'}, ...
                'value',1,'SliderStep',[1/size(CFCData{1,subj}.crsspctrm,3) 1/size(CFCData{1,subj}.crsspctrm,3)]);
            
            UpdateAxes(csquare(subj,:),Deci.Plot.CFC.Roi,'C',1);
            
            for tick = 1:size(CFCData{1,subj}.crsspctrm,3)
                uicontrol('style','text','position',[75+[[[slide.Position(3)]/5]*[tick-1]]+20 55 40 25],'String',num2str(round(CFCData{1,subj}.timelow(tick),2)));
                uicontrol('style','text','position',[75+[[[slide.Position(3)]/5]*[tick-1]]+20 30 40 25],'String',num2str(round(CFCData{1,subj}.timehigh(tick),2)));
            end
            uicontrol('style','text','position',[45 55 60 25],'String','FreqLow');
            uicontrol('style','text','position',[45 30 60 25],'String','FreqHigh');
        end
        
    end
    
    
    
end

    function ChangeDim(popup,event,Axes,Roi,Lim)
        popup.Value = round(popup.Value);
        
        for i = 1:length(Axes)
            Axes(i).Children.CData = Axes(i).UserData(:,:,popup.Value);
        end
        
        UpdateAxes(Axes,Roi,Lim,1)
        
    end

    function UpdateAxes(Axes,Roi,Lim,Userdata)
        
        if Userdata == 1
            Dats = cell2mat(arrayfun(@(c) c.UserData,Axes,'un',0));
        else
            Dats = cell2mat(arrayfun(@(c) [c.Children.([Lim 'Data'])],Axes,'un',0));
        end
        
        for Axe = 1:length(Axes(:))
            if isequal(Roi,'maxmin')
                Axes(Axe).([Lim 'Lim']) = [min(Dats(:)) max(Dats(:))];
            elseif isequal(Roi,[0 1])
                Axes(Axe).([Lim 'Lim']) = [0 1];
            elseif length(Roi) == 2 && isnumeric(Roi)
                Axes(Axe).([Lim 'Lim']) = Roi;
            end
        end
    end

    function ChangeDimTopo(popup,event,Axes,Roi)
        popup.Value = round(popup.Value);
        
        if isequal(Roi,'maxmin')
            CLim = cell2mat(arrayfun(@(c) c.UserData{3},Axes,'un',0));
        elseif isequal(Roi,[0 1])
            CLim = [0 1];
        elseif length(Roi) == 2 && isnumeric(Roi)
            CLim = Roi;
        end
        
        for Axe =  1:length(Axes)
            set(Axes(Axe).Parent, 'currentaxes', Axes(Axe))
            
            NewCross = Axes(Axe).UserData{1};
            NewCross.cohspctrm =  Axes(Axe).UserData{3}(:,:,popup.Value);
            NewCrossCfg = Axes(Axe).UserData{2};
            NewCrossCfg.CLim = minmax(CLim(:)');
            ft_topoplotCC(NewCrossCfg ,NewCross);
            
        end
        
    end

end

