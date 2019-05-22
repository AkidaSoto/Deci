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
    
    c = CleanDir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Lock filesep]);
    
    Freq.Channels       = CleanDir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Lock filesep c{1} filesep]);
    Freq.Channels       = cellfun(@(c) c(1:end-4),Freq.Channels,'un',0);
    
    if ~ischar(Deci.Plot.Freq.Channel)
        if ~any(ismember(Deci.Plot.Freq.Channel,Freq.Channels))
            error('Freq Channel Parameter has additional channels not found in Data')
        end
        
        Freq.Channels = Deci.Plot.Freq.Channel(ismember(Deci.Plot.Freq.Channel,Freq.Channels));
    end
    
    
    for Conditions = 1:length(Deci.Analysis.Conditions)
        
        for Channel = 1:length(Freq.Channels)
            
            display(['Channel ' Freq.Channels{Channel} ' of ' num2str(length(Freq.Channels))])
            
            freq = [];
            
            
            switch Deci.Plot.Freq.Type
                case 'TotalPower'
                    load([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Analysis.CondTitle{Conditions} filesep Freq.Channels{Channel} '.mat'],'freq');
                case 'ITPC'
                    load([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Analysis.CondTitle{Conditions} filesep Freq.Channels{Channel} '.mat'],'freq');
                case 'TotalPower Mean/Var'
                    load([Deci.Folder.Analysis filesep 'Freq_TotalPowerVar' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Analysis.CondTitle{Conditions} filesep Freq.Channels{Channel} '.mat'],'freq');
            end
            
            foi = freq.freq >= round(Deci.Plot.Freq.Foi(1),4) & freq.freq <= round(Deci.Plot.Freq.Foi(2),4);

            Chans{Channel} = freq;
            Chans{Channel}.freq =  Chans{Channel}.freq(foi);
            %Chans{Channel}.time =  Chans{Channel}.time(toi);
            Chans{Channel}.powspctrm  =Chans{Channel}.powspctrm(:,foi,:);
            Chans{Channel}.label = Freq.Channels(Channel);
            
        end
        toc;
        
        acfg.parameter = 'powspctrm';
        acfg.appenddim = 'chan';
        Subjects{subject_list,Conditions} = rmfield(ft_appendfreq(acfg,Chans{:}),'cfg');
        
        acfg.latency = Deci.Plot.Freq.Bsl;
        acfg.avgovertime = 'yes';
        
        Subjects{subject_list,Conditions}.dimord = 'chan_freq_time';
        
        
    end
    clear Chans;
end


for Conditions = 1:size(Subjects,2)
    for subject_list = 1:size(Subjects,1)

        toi = round(Subjects{subject_list,Conditions}.time,4) >= Deci.Plot.Freq.Toi(1) & round(Subjects{subject_list,Conditions}.time,4) <= Deci.Plot.Freq.Toi(2);
        
        bsl = ft_selectdata(acfg, Subjects{subject_list,Conditions});
        
        Subjects{subject_list,Conditions}.powspctrm =  Subjects{subject_list,Conditions}.powspctrm(:,:,toi);
        Subjects{subject_list,Conditions}.time = Subjects{subject_list,Conditions}.time(toi);
        bsl = repmat(bsl.powspctrm,[1 1 size(Subjects{subject_list,Conditions}.powspctrm ,3)]);
        
        switch Deci.Plot.Freq.BslType
            case 'none'
                
            case 'absolute'
                Subjects{subject_list,Conditions}.powspctrm =  Subjects{subject_list,Conditions}.powspctrm - bsl;
            case 'relative'
                Subjects{subject_list,Conditions}.powspctrm=  Subjects{subject_list,Conditions}.powspctrm ./ bsl;
            case 'relchange'
                Subjects{subject_list,Conditions}.powspctrm = ( Subjects{subject_list,Conditions}.powspctrm - bsl) ./ bsl;
            case 'db'
                Subjects{subject_list,Conditions}.powspctrm = 10*log10( Subjects{subject_list,Conditions}.powspctrm ./ bsl);
        end
        
    end
end
        
if ~isempty(Deci.Plot.Math)
    for cond = 1:length(Deci.Plot.Math)
        for subj = 1:size(Subjects,1)
            scfg.parameter = 'powspctrm';
            scfg.operation = Deci.Plot.Math{cond};
            MathData{subj} = ft_math(scfg,Subjects{subj,:});
        end 
        Subjects(:,size(Subjects,2)+1) = MathData;
        %TrialCount(:,size(Subjects,2)+cond) = num2cell(nan(size(Subjects,1),1));
    end
end

if Deci.Plot.Freq.Wires
    
    for conds = 1:size(Subjects,2)
        for subj = 1:size(Subjects,1)
            bcfg.frequency = Deci.Plot.Freq.Foi;
            
            switch Deci.Plot.Freq.Ws.avg
                case 'freq'
                    bcfg.avgoverfreq = 'yes';
                case 'time'
                    bcfg.avgovertime = 'yes';
            end
            bcfg.nanmean = 'yes';
            WireSub{subj,conds} = ft_selectdata(bcfg,Subjects{subj,conds});
            
        end
    end
end


if Deci.Plot.GrandAverage
    
    if Deci.Plot.Freq.Wires
        save([Deci.Folder.Plot filesep 'preGrandAverageSubWire'],'WireSub')
    end
    
    for conds = 1:size(Subjects,2)
        facfg.parameter =  'powspctrm';
        facfg.type = 'mean';
        FreqData{conds} = rmfield(ft_freqgrandaverage(facfg,Subjects{:,conds}),'cfg');
        
        %TotalCount{conds} = mean([TrialCount{:,conds}]);
        
        if Deci.Plot.Freq.Wires
            facfg.type = 'mean';
            WireData{conds} = rmfield(ft_freqgrandaverage(facfg,WireSub{:,conds}),'cfg');
            
            facfg.type = 'sem';
            WireStd{conds} = rmfield(ft_freqgrandaverage(facfg,WireSub{:,conds}),'cfg');
        end
        
    end
    
    save([Deci.Folder.Plot filesep 'GrandAverageSubFreq'],'FreqData')
    
    if ~isempty(Deci.Plot.Freq.Wires)
        save([Deci.Folder.Plot filesep 'GrandAverageSubWire'],'WireData')
    end
else
    FreqData = Subjects;
    WireData = WireSub;
    %TotalCount = TrialCount;
    Deci.Plot.Freq.Ws.errorbars = 0;
end
clear Subjects;

if strcmpi(Deci.Plot.FreqYScale,'log')
    
    for i = 1:size(FreqData,1)
        for j = 1:size(FreqData,2)
            FreqData{i,j}.freq = log(FreqData{i,j}.freq);
        end
    end
end



%% Plot


if Deci.Plot.GrandAverage
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
        
        if Deci.Plot.Freq.Wires
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
                cfg.clim = 'maxmin';
                ft_singleplotTFR(cfg,FreqData{subj,Deci.Plot.Draw{cond}(subcond)})
                
                %title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' ' Deci.Plot.Subtitle{cond}{subcond} ' trial count: ' num2str(TotalCount{subj,Deci.Plot.Draw{cond}(subcond)})]);
                title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' ' Deci.Plot.Subtitle{cond}{subcond}]);

                colorbar('vert');
                map = colormap('jet'); %'hot' 'gray'
                colormap(map);
            end
            
            
            if Deci.Plot.Freq.Wires
                
                set(0, 'CurrentFigure', wire(subj) )
                wire(subj).Visible = 'on';
                
                if strcmpi(Deci.Plot.Freq.Ws.avg,'freq')
                    
                    if ~Deci.Plot.Freq.Ws.errorbars
                        
                        h = plot(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.time,squeeze(mean(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1)));
                    else
                        top = squeeze(nanmean(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1)) + squeeze(nanmean(WireStd{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1));
                        bot = squeeze(nanmean(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1)) - squeeze(nanmean(WireStd{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1));
                        
                        pgon = polyshape([WireData{subj,Deci.Plot.Draw{cond}(subcond)}.time fliplr(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.time)],[top' fliplr(bot')],'Simplify', false);
                        b = plot(pgon,'HandleVisibility','off');
                        hold on
                        b.EdgeAlpha = 0;
                        b.FaceAlpha = .15;
                        h = plot(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.time,squeeze(mean(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1)));
                        h.Color = b.FaceColor;
                        h.LineWidth = 1;
                        
                        
                    end
                    xlabel('Time');
                elseif  strcmpi(Deci.Plot.Freq.Ws.avg,'time')
                    
                    if strcmpi(Deci.Plot.Scale,'log')
                        WireData{subj,Deci.Plot.Draw{cond}(subcond)}.freq = exp(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.freq);
                    end
                    
                    h =  plot(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.freq,squeeze(mean(WireData{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1)));
                    xlabel('Freq')
                end
                ylabel('Power')
                hold on
                
                legend(h.Parent,Deci.Plot.Subtitle{cond})
                title(h.Parent,[Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' ' Deci.Plot.Title{cond} ' Wire'])
                
                if ~isempty(Deci.Folder.Plot)
                    mkdir([Deci.Folder.Plot filesep Deci.Plot.Title{cond}]);
                    saveas(wire(subj),[Deci.Folder.Plot filesep Deci.Plot.Title{cond} filesep Deci.SubjectList{subj} '_wire'],Deci.Plot.Save.Format);
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
                    mkdir([Deci.Folder.Plot filesep Deci.Plot.Title{cond}]);
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
            
            if strcmpi(Deci.Plot.FreqYScale,'log')
                for r = 1:length(subby(:))
                    subby(r).YTickLabels = exp(subby(r).YTick);
                end
            end
            
            if ~isempty(Deci.Folder.Plot)
                mkdir([Deci.Folder.Plot filesep Deci.Plot.Title{cond}]);
                saveas(square(subj),[Deci.Folder.Plot filesep Deci.Plot.Title{cond} filesep Deci.SubjectList{subj} '_square'],Deci.Plot.Save.Format);
            end
            
        end
        
        
    end
    
end

end
