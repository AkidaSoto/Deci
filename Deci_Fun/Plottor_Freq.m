function Plottor_Freq(Deci)


cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';


%% Deci Checks
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

if length(Deci.Plot.Freq.Channel) == 1
    Deci.Plot.Freq.Topo = 0;
end

%% Load

for  subject_list = 1:length(Deci.SubjectList)
    
    tic;
    
    Freq.Channels       = CleanDir([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Lock]);
    Freq.Channels       = cellfun(@(c) c(1:end-4),Freq.Channels,'un',0);
    
    if ~ischar(Deci.Plot.Freq.Channel)
        if ~any(ismember(Deci.Plot.Freq.Channel,Freq.Channels))
            error('Freq Channel Parameter has additional channels not found in Data')
        end
        
        Freq.Channels = Deci.Plot.Freq.Channel(ismember(Deci.Plot.Freq.Channel,Freq.Channels));
    end
    
    for Channel = 1:length(Freq.Channels)
        
        display(['Channel ' Freq.Channels{Channel} ' of ' num2str(length(Freq.Channels))])
        
        freq = [];
        load([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Lock filesep Freq.Channels{Channel}],'freq');
        
        foi = freq.freq >= round(Deci.Plot.Freq.Foi(1),4) & freq.freq <= round(Deci.Plot.Freq.Foi(2),4);
        %toi = round(freq.time,4) >= Deci.Plot.Freq.Toi(1) & round(freq.time,4) <= Deci.Plot.Freq.Toi(2);
        
        Chans{Channel} = freq;
        Chans{Channel}.freq =  Chans{Channel}.freq(foi);
        %Chans{Channel}.time =  Chans{Channel}.time(toi);
        Chans{Channel}.powspctrm  =Chans{Channel}.fourierspctrm(:,:,foi,:);
        Chans{Channel}.label = Freq.Channels(Channel);
        
    end
    toc;
    
    acfg.parameter = 'powspctrm';
    acfg.appenddim = 'chan';
    Chans = rmfield(ft_appendfreq(acfg,Chans{:}),'cfg');
    
    acfg.latency = Deci.Plot.Freq.Bsl;
    acfg.avgovertime = 'yes';
    
    
    if Deci.Plot.Var
        a = figure;
        for Channel = 1:length(Chans.label)
            
            plot(reshape(squeeze(10*log10(mean(mean(abs(Chans.powspctrm(:,Channel,:,:)).^2,3),2)))',[1 size(Chans.powspctrm,1)*size(Chans.powspctrm,4)]));
            hold on
        end
        a.Visible = 'on';
        title([Deci.SubjectList{subject_list} ' Chan']);
        
        waitfor(a)
        
        b = figure;
        for Channel = 1:length(Chans.freq)
            
            plot(reshape(squeeze(10*log10(mean(mean(abs(Chans.powspctrm(:,:,Channel,:)).^2,3),2)))',[1 size(Chans.powspctrm,1)*size(Chans.powspctrm,4)]));
            hold on
        end
        b.Visible = 'on';
        title([Deci.SubjectList{subject_list} ' Freq'])
        waitfor(b)
        
    end
    
    
    
    for Conditions = 1:length(Deci.Plot.Conditions)
        
        maxt = max(sum(ismember(freq.condinfo{2},Deci.Plot.Conditions{Conditions}),2));
        
        trl = sum(ismember(freq.condinfo{2},Deci.Plot.Conditions{Conditions}),2) == maxt;
        
        toi = round(Chans.time,4) >= Deci.Plot.Freq.Toi(1) & round(Chans.time,4) <= Deci.Plot.Freq.Toi(2);
        
        
        
        Subjects{subject_list,Conditions} = Chans;
        Subjects{subject_list,Conditions}.powspctrm = Chans.powspctrm(trl,:,:,:);
        
        TrialCount{subject_list,Conditions} = length(find(trl));
        
        switch Deci.Plot.Freq.Type
            case 'TotalPower'
                Subjects{subject_list,Conditions}.powspctrm = permute(mean(abs(Subjects{subject_list,Conditions}.powspctrm).^2 ,1),[2 3 4 1]);
            case 'ITPC'
                Subjects{subject_list,Conditions}.powspctrm = permute(abs(mean(Subjects{subject_list,Conditions}.powspctrm./abs(Subjects{subject_list,Conditions}.powspctrm),1)),[2 3 4 1]) ;
            case 'TotalPower Std'
                Subjects{subject_list,Conditions}.powspctrm = permute(std(abs(Subjects{subject_list,Conditions}.powspctrm).^2 ,1),[2 3 4 1]);
            case 'TotalPower Mean/Std'
                Subjects{subject_list,Conditions}.powspctrm = permute(mean(abs(Subjects{subject_list,Conditions}.powspctrm).^2 ,1),[2 3 4 1])./permute(std(abs(Subjects{subject_list,Conditions}.powspctrm).^2 ,1),[2 3 4 1]);
                
        end

        bsl = ft_selectdata(acfg, Subjects{subject_list,Conditions});
        
        Subjects{subject_list,Conditions}.powspctrm =  Subjects{subject_list,Conditions}.powspctrm(:,:,toi);
        Subjects{subject_list,Conditions}.time = Subjects{subject_list,Conditions}.time(toi);
        bsl = repmat(bsl.powspctrm,[1 1 size(Subjects{subject_list,Conditions}.powspctrm ,3)]);
        
        
        switch Deci.Plot.Freq.BslType
            case 'absolute'
                Subjects{subject_list,Conditions}.powspctrm =  Subjects{subject_list,Conditions}.powspctrm - bsl;
            case 'relative'
                Subjects{subject_list,Conditions}.powspctrm=  Subjects{subject_list,Conditions}.powspctrm ./ bsl;
            case 'relchange'
                Subjects{subject_list,Conditions}.powspctrm = ( Subjects{subject_list,Conditions}.powspctrm - bsl) ./ bsl;
            case 'db'
                Subjects{subject_list,Conditions}.powspctrm = 10*log10( Subjects{subject_list,Conditions}.powspctrm ./ bsl);
        end
        
        
        Subjects{subject_list,Conditions}.dimord = 'chan_freq_time';
    end
    clear Chans;
end



if Deci.Plot.GA
    
    for conds = 1:size(Subjects,2)
        facfg.parameter =  'powspctrm';
        FreqData{conds} = rmfield(ft_freqgrandaverage(facfg,Subjects{:,conds}),'cfg');
        
        TotalCount{conds} = mean([TrialCount{:,conds}]);
        
    end
else
    FreqData = Subjects;
    TotalCount = TrialCount;
end
clear Subjects;


if ~isempty(Deci.Plot.PRP)
    if strcmp(Deci.Plot.PRP.Dim,'trials')
        for  Condition =  1:length(Deci.Plot.Conditions)
            for  subject_list = 1:length(Deci.SubjectList)
                for Channel = 1:length(Freq.Channels)
                    
                    display(['Channel ' Freq.Channels{Channel} ' of ' num2str(length(Freq.Channels))])
                    
                    freq = [];
                    load([Deci.Folder.Analysis filesep 'Four_' Deci.Plot.Freq.Type filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Conditions{Condition} filesep Freq.Channels{Channel}],'freq');
                    
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


if ~isempty(Deci.Plot.Math)
    
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
            for subj = 1:size(FreqData,1)
                scfg.parameter = 'powspctrm';
                scfg.operation = Deci.Plot.Math.Form{cond};
                MathData{subj} = ft_math(scfg,FreqData{subj,:});
            end
            
            FreqData(:,length(Deci.Plot.Conditions)+cond) = MathData;
            TotalCount(:,length(Deci.Plot.Conditions)+cond) = num2cell(nan(size(FreqData,1),1));
        end
        
        %             if Deci.Plot.Math.Type == 1
        %                 FreqData = FreqData(:,length(Deci.Plot.Conditions)+1:end);
        %                 TotalCount = TotalCount(:,length(Deci.Plot.Conditions)+1:end);
        %             end
        
    end
    
    
    
end


if strcmpi(Deci.Plot.Scale,'log')
    
    for i = 1:size(FreqData,1)
        for j = 1:size(FreqData,2)
            FreqData{i,j}.freq = log(FreqData{i,j}.freq);
        end
    end
end


if ~isempty(Deci.Plot.Freq.Wires)
    
    switch Deci.Plot.Freq.Wires.avg
        case 'freq'
            wcfg.avgoverfreq = 'yes';
        case 'time'
            wcfg.avgovertime = 'yes';
    end
    
    wcfg.avgoverchan = 'yes';
    
    for subj = 1:size(FreqData,1)
        for cond = 1:size(FreqData,2)
            WireData{subj,cond} =  ft_selectdata(wcfg,FreqData{subj,cond});
        end
    end
end


%% Plot


if Deci.Plot.GA
    Deci.SubjectList = {'Group 1'};
end

for cond = 1:length(Deci.Plot.Draw)
    for subj = 1:size(FreqData,1)
        
        if Deci.Plot.Freq.Square
            square(subj) = figure;
        end
        
        if Deci.Plot.Freq.Topo
            topo(subj)  = figure;
        end
        
        if ~isempty(Deci.Plot.Freq.Wires)
            wire(subj)  = figure;
        end
        
        for subcond = 1:length(Deci.Plot.Draw{cond})
            
            if Deci.Plot.Freq.Topo
                if length(Freq.Channels) ~= 1
                    
                    
                    set(0, 'CurrentFigure', topo(subj))
                    topo(subj).Visible = 'on';
                    cirky(subj,subcond)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                    ft_topoplotER(cfg, FreqData{subj,Deci.Plot.Draw{cond}(subcond)});
                    
                    title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' '  Deci.Plot.Subtitle{cond}{subcond} ' trial count: ' num2str(TotalCount{subj,Deci.Plot.Draw{cond}(subcond)})]);
                    colorbar('vert');
                    map = colormap('jet'); %'hot' 'gray'
                    colormap(map);
                    
                else
                    close topo
                end
                
                
                
            end
            
            if Deci.Plot.Freq.Square
                set(0, 'CurrentFigure', square(subj) )
                square(subj).Visible = 'on';
                subby(subj,subcond) = subplot(length(Deci.Plot.Draw{cond}),1,subcond );
                ft_singleplotTFR(cfg,FreqData{subj,Deci.Plot.Draw{cond}(subcond)})
                
                title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' ' Deci.Plot.Subtitle{cond}{subcond} ' trial count: ' num2str(TotalCount{subj,Deci.Plot.Draw{cond}(subcond)})]);
                colorbar('vert');
                map = colormap('jet'); %'hot' 'gray'
                colormap(map);
            end
            
           
            if ~isempty(Deci.Plot.Freq.Wires)
                
                set(0, 'CurrentFigure', wire(subj) )
                wire(subj).Visible = 'on';
                 
                if strcmpi(Deci.Plot.Freq.Wires.avg,'freq')
                    h = plot(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.time,squeeze(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm));
                    xlabel('Time');
                elseif  strcmpi(Deci.Plot.Freq.Wires.avg,'time')
                    
                    if strcmpi(Deci.Plot.Scale,'log')
                        WireData{subj,Deci.Plot.Draw{cond}(subcond)}.freq = exp(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.freq);
                    end
                    
                   h =  plot(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.freq,squeeze(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm));
                    xlabel('Freq')
                end
                ylabel('Power')
                hold on

                legend(h.Parent,Deci.Plot.Subtitle{cond})
               title(h.Parent,[Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' ' Deci.Plot.Title{cond} ' Wire'])
                
                if ~isempty(Deci.Folder.Plot)
                    saveas(wire(subj),[Deci.Folder.Plot filesep Deci.SubjectList{subj} '_wire'],Deci.Plot.Save.Format);
                end
                
            end
            
            
        end
        
    end
    
    for subj = 1:size(FreqData,1)
       
        
        if length(Freq.Channels) ~= 1
            
            if Deci.Plot.Freq.Topo
                set(0, 'CurrentFigure', topo(subj) )
                topo(subj).Visible = 'on';
                suptitle(Deci.Plot.Title{cond});
                
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
                    mkdir([Deci.Folder.Plot filesep Deci.Plot.Title{cond} filesep Deci.SubjectList{subj} '_topo']);
                    saveas(topo(subj),[Deci.Folder.Plot filesep Deci.Plot.Title{cond} filesep Deci.SubjectList{subj} '_topo'],Deci.Plot.Save.Format);
                end
                
            end
            
            
        end
        
        if Deci.Plot.Freq.Square
            set(0, 'CurrentFigure', square(subj))
            square(subj).Visible = 'on';
            suptitle(Deci.Plot.Title{cond});
            
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
                mkdir([Deci.Folder.Plot filesep Deci.Plot.Title{cond} filesep Deci.SubjectList{subj} '_square']);
                saveas(square(subj),[Deci.Folder.Plot filesep Deci.Plot.Title{cond} filesep Deci.SubjectList{subj} '_square'],Deci.Plot.Save.Format);
            end
            
        end
        
        
    end
    
end

end
