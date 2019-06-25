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
    
    
    Freq.Channels       = CleanDir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Lock filesep Deci.Plot.IndexTitle{1} filesep]);
    Freq.Channels       = cellfun(@(c) c(1:end-4),Freq.Channels,'un',0);
    
    if ~ischar(Deci.Plot.Freq.Channel)
        if ~any(ismember(Deci.Plot.Freq.Channel,Freq.Channels))
            error('Freq Channel Parameter has additional channels not found in Data')
        end
        Freq.Channels = Deci.Plot.Freq.Channel(ismember(Deci.Plot.Freq.Channel,Freq.Channels));
    end
    
    
    for Conditions = 1:length(Deci.Plot.IndexTitle)
        
        for Channel = 1:length(Freq.Channels)
            
            display(['Channel ' Freq.Channels{Channel} ' of ' num2str(length(Freq.Channels))])
            
            freq = [];
            
            switch Deci.Plot.Freq.Type
                case 'TotalPower'
                    load([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.IndexTitle{Conditions} filesep Freq.Channels{Channel} '.mat'],'freq');
                case 'ITPC'
                    load([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.IndexTitle{Conditions} filesep Freq.Channels{Channel} '.mat'],'freq');
                case 'TotalPower Mean/Var'
                    load([Deci.Folder.Analysis filesep 'Freq_TotalPowerVar' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.IndexTitle{Conditions} filesep Freq.Channels{Channel} '.mat'],'freq');
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
        Subjects{subject_list,Conditions}.dimord = 'chan_freq_time';
    end
    clear Chans;
end

%% Baseline Correction
for Conditions = 1:size(Subjects,2)
    for subject_list = 1:size(Subjects,1)

        if ~strcmpi(Deci.Plot.BslRef,Deci.Plot.Lock)
            
            for Channel = 1:length(Freq.Channels)

                freq = [];
                    
                switch Deci.Plot.Freq.Type
                    case 'TotalPower'
                        load([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.BslRef filesep Deci.Plot.IndexTitle{Conditions} filesep Freq.Channels{Channel} '.mat'],'freq');
                    case 'ITPC'
                        load([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.BslRef filesep Deci.Plot.IndexTitle{Conditions} filesep Freq.Channels{Channel} '.mat'],'freq');
                    case 'TotalPower Mean/Var'
                        load([Deci.Folder.Analysis filesep 'Freq_TotalPowerVar' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.BslRef filesep Deci.Plot.IndexTitle{Conditions} filesep Freq.Channels{Channel} '.mat'],'freq');
                end
                
                foi = freq.freq >= round(Deci.Plot.Freq.Foi(1),4) & freq.freq <= round(Deci.Plot.Freq.Foi(2),4);
                Chans{Channel} = freq;
                Chans{Channel}.freq =  Chans{Channel}.freq(foi);
                %Chans{Channel}.time =  Chans{Channel}.time(toi);
                Chans{Channel}.powspctrm  =Chans{Channel}.powspctrm(:,foi,:);
                Chans{Channel}.label = Freq.Channels(Channel);
                
            end
            
            acfg.parameter = 'powspctrm';
            acfg.appenddim = 'chan';
            Bsl{subject_list,Conditions} = rmfield(ft_appendfreq(acfg,Chans{:}),'cfg');
            
            ccfg.latency = Deci.Plot.Freq.Bsl;
            ccfg.avgovertime = 'yes';
            
            toi = round(Bsl{subject_list,Conditions}.time,4) >= Deci.Plot.Freq.Toi(1) & round(Bsl{subject_list,Conditions}.time,4) <= Deci.Plot.Freq.Toi(2);
            bsl = ft_selectdata(ccfg, Bsl{subject_list,Conditions});
            
            Bsl{subject_list,Conditions}.powspctrm =  Bsl{subject_list,Conditions}.powspctrm(:,:,toi);
            Bsl{subject_list,Conditions}.time = Bsl{subject_list,Conditions}.time(toi);
            bsl = repmat(bsl.powspctrm,[1 1 size(Bsl{subject_list,Conditions}.powspctrm ,3)]);
            
        else
            
            ccfg.latency = Deci.Plot.Freq.Bsl;
            ccfg.avgovertime = 'yes';
            
            toi = round(Subjects{subject_list,Conditions}.time,4) >= Deci.Plot.Freq.Toi(1) & round(Subjects{subject_list,Conditions}.time,4) <= Deci.Plot.Freq.Toi(2);
            bsl = ft_selectdata(ccfg, Subjects{subject_list,Conditions});
            
            Subjects{subject_list,Conditions}.powspctrm =  Subjects{subject_list,Conditions}.powspctrm(:,:,toi);
            Subjects{subject_list,Conditions}.time = Subjects{subject_list,Conditions}.time(toi);
            bsl = repmat(bsl.powspctrm,[1 1 size(Subjects{subject_list,Conditions}.powspctrm ,3)]);
        end
        

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

%% Math
if ~isempty(Deci.Plot.Math)
    for conds = 1:length(Deci.Plot.Math)
        for subj = 1:size(Subjects,1)
            scfg.parameter = 'powspctrm';
            scfg.operation = Deci.Plot.Math{conds};
            MathData{subj} = ft_math(scfg,Subjects{subj,:});
        end 
        Subjects(:,size(Subjects,2)+1) = MathData;
    end
end

%% Data Management

if size(Subjects,1) == 1
    Deci.Plot.GrandAverage = false;  
end

for conds = 1:size(Subjects,2)
    
    if Deci.Plot.GrandAverage
        facfg.parameter =  'powspctrm';
        facfg.type = 'mean';
        facfg.keepindividual = 'yes';
        FreqData{1,conds} = rmfield(ft_freqgrandaverage(facfg,Subjects{:,conds}),'cfg');
    else
        Deci.Plot.Stat.do = false;
        FreqData(:,conds) = Subjects(:,conds);
    end
    
    for subj = 1:size(FreqData(:,conds),1)
        
        if Deci.Plot.Freq.Topo
            tcfg = [];
            tcfg.avgoverfreq = 'yes';
            tcfg.avgovertime = 'yes';
            tcfg.nanmean = 'yes';
            topodata{subj,conds} = ft_selectdata(tcfg,FreqData{subj,conds});
        end
        
        if Deci.Plot.Freq.Square
            tcfg = [];
            tcfg.avgoverchan = 'yes';
            tcfg.nanmean = 'yes';
            squaredata{subj,conds} = ft_selectdata(tcfg,FreqData{subj,conds});
        end
        
        if Deci.Plot.Freq.Wire
            tcfg = [];
            tcfg.avgoverchan = 'yes';
            tcfg.avgoverfreq = 'yes';
            tcfg.nanmean = 'yes';
            wiredata{subj,conds} = ft_selectdata(tcfg,FreqData{subj,conds});
        end
        
        if Deci.Plot.Freq.Bar
            tcfg = [];
            tcfg.avgoverchan = 'yes';
            tcfg.avgoverfreq = 'yes';
            tcfg.avgovertime = 'yes';
            tcfg.nanmean = 'yes';
            bardata{subj,conds} = ft_selectdata(tcfg,FreqData{subj,conds});
        end
    end
end

if Deci.Plot.Stat.do
    neighbours       = load('easycap_neighbours','neighbours');
    Deci.Plot.Stat.neighbours = neighbours.neighbours;
    Deci.Plot.Stat.ivar = 1;
    Deci.Plot.Stat.uvar = 2;
    Deci.Plot.Stat.tail = 1;
    Deci.Plot.Stat.statistic = 'depsamplesFmultivariate';
        
    for conds = 1:length(Deci.Plot.Draw)
        design = [];
        
        for subcond = 1:length(Deci.Plot.Draw{conds})
            for subj = 1:size(Subjects,1)
                design(1,subj+size(Subjects,1)*[subcond-1]) =  subcond;
                design(2,subj+size(Subjects,1)*[subcond-1]) = subj;
            end
        end
        
        Deci.Plot.Stat.design = design;
        
        if Deci.Plot.Freq.Topo
            [topostat{conds}] = ft_freqstatistics(Deci.Plot.Stat, topodata{Deci.Plot.Draw{conds}});
            if Deci.Plot.Stat.fdr
                topostat{conds}.prob =  fdr(topostat{conds}.prob,Deci.Plot.Stat.alpha);
            end
        end
        
        if Deci.Plot.Freq.Square
            [squarestat{conds}] = ft_freqstatistics(Deci.Plot.Stat, squaredata{Deci.Plot.Draw{conds}});
            if Deci.Plot.Stat.fdr
                squarestat{conds}.prob =  fdr(squarestat{conds}.prob,Deci.Plot.Stat.alpha);
            end
        end
        
        if Deci.Plot.Freq.Wire
            [wirestat{conds}] = ft_freqstatistics(Deci.Plot.Stat, wiredata{Deci.Plot.Draw{conds}});
            if Deci.Plot.Stat.fdr
                wirestat{conds}.prob =  fdr(wirestat{conds}.prob,Deci.Plot.Stat.alpha);
            end       
        end
        
        if Deci.Plot.Freq.Bar
            [barstat{conds}] = ft_freqstatistics(Deci.Plot.Stat, bardata{Deci.Plot.Draw{conds}});
            if Deci.Plot.Stat.fdr
                barstat{conds}.prob =  fdr(barstat{conds}.prob,Deci.Plot.Stat.alpha);
            end
        end
    end
end

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
    if Deci.Plot.Stat.do
        
        if Deci.Plot.Freq.Square
            tacfg = cfg;
            tacfg.parameter = 'stat';
            
            mask = double(squarestat{cond}.prob < Deci.Plot.Stat.alpha);
            mask(mask == 0) = .2;
            squarestat{cond}.mask = mask;
            
            tacfg.maskparameter = 'mask';
            
            if Deci.Plot.Stat.FPlots
                squaret(cond) = figure;
                squaret(cond).Visible = 'on';
                
                ft_singleplotTFR(tacfg,squarestat{cond})
                title([Deci.Plot.Stat.Type ' ' Deci.Plot.Title{cond} ' Square (alpha = ' num2str(Deci.Plot.Stat.alpha) ')']);
                colormap('jet'); %'hot' 'gray'
                ylabel('F score');
                xlabel('time');
            end
        end
        
        if Deci.Plot.Freq.Topo
            tcfg = cfg;
            tcfg.parameter = 'stat';
            
            mask = double(topostat{cond}.prob < Deci.Plot.Stat.alpha);
            mask(mask == 0) = .2;
            topostat{cond}.mask = mask;
            
            if Deci.Plot.Stat.FPlots
                topot(cond)  = figure;
                topot(cond).Visible = 'on';
                
                ft_topoplotER(tcfg, topostat{cond});
                title([Deci.Plot.Stat.Type ' ' Deci.Plot.Title{cond} ' Square (alpha = ' num2str(Deci.Plot.Stat.alpha) ')']);
                colormap('jet'); %'hot' 'gray'
            end
        end
        
        if Deci.Plot.Freq.Wire
            
            tcfg = cfg;
            tcfg.parameter = 'stat';
            
            mask = double(wirestat{cond}.prob < Deci.Plot.Stat.alpha);
            %mask(mask == 0) = .2;
            wirestat{cond}.mask = mask;
            
            if Deci.Plot.Stat.FPlots
                wiret(cond)  = figure;
                wiret(cond).Visible = 'on';
                
                plot(barstat{cond}.time,barstat{cond}.stat)
                title([Deci.Plot.Stat.Type ' ' Deci.Plot.Title{cond} ' Square (alpha = ' num2str(Deci.Plot.Stat.alpha) ')']);
            end
            
        end
        
        if Deci.Plot.Freq.Bar
            tcfg = cfg;
            tcfg.parameter = 'stat';
            
            mask = double(barstat{cond}.prob < Deci.Plot.Stat.alpha);
            mask(mask == 0) = .2;
            barstat{cond}.mask = mask;
            
            if Deci.Plot.Stat.FPlots
                bart(cond)  = figure;
                bart(cond).Visible = 'on';
                
                bar(barstat{cond}.stat)
                title([Deci.Plot.Stat.Type ' ' Deci.Plot.Title{cond} ' Square (alpha = ' num2str(Deci.Plot.Stat.alpha) ')']);
                colormap('jet'); %'hot' 'gray'
            end
        end
        
    end 

    for subj = 1:size(FreqData,1)
            
        if Deci.Plot.Freq.Square
            square(subj) = figure;
            
            if Deci.Plot.Stat.do
                ButtonH=uicontrol('Parent', square(subj),'Style','pushbutton','String','p Mask','Position',[10 10 100 25],'Visible','on','Callback',@pmask);
                ButtonH.UserData = @ones;
            end
        end
        
        if Deci.Plot.Freq.Topo
            topo(subj)  = figure;
            
            if Deci.Plot.Stat.do
                ButtonH=uicontrol('Parent', topo(subj),'Style','pushbutton','String','p Mask','Position',[10 10 100 25],'Visible','on','Callback',@pmask);
                ButtonH.UserData = @ones;
            end
        end
        
        if Deci.Plot.Freq.Wire
            wire(subj)  = figure;
        end

        %freq plots
        clear b;
        for subcond = 1:length(Deci.Plot.Draw{cond})
            
            if Deci.Plot.Freq.Topo
                if length(Freq.Channels) ~= 1
                    
                    
                    set(0, 'CurrentFigure', topo(subj))
                    topo(subj).Visible = 'on';
                    cirky(subj,subcond)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                    
                    pcfg = cfg;
                    if Deci.Plot.Stat.do
                        
                        pcfg.clim = 'maxmin';
                        pcfg.maskparameter ='mask';
                        FreqData{subj,Deci.Plot.Draw{cond}(subcond)}.mask = repmat(topostat{cond}.mask,[1 length(FreqData{subj,Deci.Plot.Draw{cond}(subcond)}.freq) length(FreqData{subj,Deci.Plot.Draw{cond}(subcond)}.time)]);
                    end
                    
                    ft_topoplotER(pcfg, FreqData{subj,Deci.Plot.Draw{cond}(subcond)});
                    
                    %title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' '  Deci.Plot.Subtitle{cond}{subcond} ' trial count: ' num2str(TotalCount{subj,Deci.Plot.Draw{cond}(subcond)})]);
                    title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' '  Deci.Plot.Subtitle{cond}{subcond}]);

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
                
                
                pcfg = cfg;
                if Deci.Plot.Stat.do
                    
                    pcfg.clim = 'maxmin';
                    pcfg.maskparameter ='mask';
                    FreqData{subj,Deci.Plot.Draw{cond}(subcond)}.mask = repmat(squarestat{cond}.mask,[length(FreqData{subj,Deci.Plot.Draw{cond}(subcond)}.label) 1 1]);
                end
                
                ft_singleplotTFR(pcfg,FreqData{subj,Deci.Plot.Draw{cond}(subcond)});
                title([Deci.Plot.Subtitle{cond}{subcond}]);

                colorbar('vert');
                map = colormap('jet'); %'hot' 'gray'
                colormap(map);
            end
            
            if Deci.Plot.Freq.Wire
                set(0, 'CurrentFigure', wire(subj) )
                wire(subj).Visible = 'on';
                
                top = squeeze(nanmean(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1)) + squeeze(nanstd(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1)/sqrt(size(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1)));
                bot = squeeze(nanmean(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1)) - squeeze(nanstd(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1)/sqrt(size(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1)));
                
                pgon = polyshape([wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.time fliplr(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.time)],[top' fliplr(bot')],'Simplify', false);
                b(subcond) = plot(pgon,'HandleVisibility','off');
                hold on
                b(subcond).EdgeAlpha = 0;
                b(subcond).FaceAlpha = .15;
                
                if Deci.Plot.Stat.do
                wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.mask = repmat(wirestat{cond}.mask,[length(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.label) 1 1]);
                end
            end
            
        end
        
        if Deci.Plot.Freq.Wire
            set(0, 'CurrentFigure', wire(subj) )
            wire(subj).Visible = 'on';
            
            
            pcfg = cfg;
            pcfg.clim = 'maxmin';
            
            if Deci.Plot.Stat.do
                pcfg.maskparameter ='mask';
                
            end
            pcfg.ylim = ylim;
            pcfg.graphcolor = lines;
            pcfg.linewidth = 1;
            ft_singleplotER(pcfg,wiredata{subj,Deci.Plot.Draw{cond}});
            
            arrayfun(@(c) uistack(c,'top'),b);
            
            axis tight
            hold on
            plot([wiredata{cond}.time(1), wiredata{cond}.time(end)], [0 0], 'k--') % hor. line
            plot([0 0], ylim, 'k--') % vert. l
            
            if Deci.Plot.Stat.do
                
                boxes = wire(subj).Children(2).Children.findobj('Type','Patch');
                
                if ~isempty(boxes)
                boxes.FaceAlpha = .35;
                uistack(boxes,'bottom')
                boxes.HandleVisibility = 'off';
                end
            end
            title([Deci.Plot.Subtitle{cond}{subcond}]);
            
            %h = plot(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.time,squeeze(mean(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1)));
            %h.Color = b.FaceColor;
            %h.LineWidth = 1;
            
            
            legend(Deci.Plot.Subtitle{cond})
            title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' ' Deci.Plot.Title{cond} ' Wire'])
            
            %                         if Deci.Plot.Freq.Wire &&subcond == 1
            %                          ax1 = axes;
            %                          hold on
            %                          imagesc(ax1,wirestatdraw{cond}.time,mean(ylim)/2,squeeze(wirestatdraw{cond}.prob <= .05),'alphadata',squeeze(wirestatdraw{cond}.prob <= .05))
            %                          ax1.Visible = 'off';
            %
            %                         end
            xlabel('Time');
        end
        
        
        if Deci.Plot.Freq.Bar
            
            barfig(subj)  = figure;
            
            %             if Deci.Plot.Stat.do
            %                 ButtonH=uicontrol('Parent', barfig(subj),'Style','pushbutton','String','p Mask','Position',[10 10 100 25],'Visible','on','Callback',@pmask);
            %                 ButtonH.UserData = 0;
            %             end
            
            set(0, 'CurrentFigure', barfig(subj) )
            barfig(subj).Visible = 'on';
            
            
            CleanBars(mean(cell2mat(arrayfun(@(c) c.powspctrm,[bardata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)),1), ...
                nanstd(cell2mat(arrayfun(@(c) c.powspctrm,[bardata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)),1) ...
                /sqrt(size(cell2mat(arrayfun(@(c) c.powspctrm,[bardata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)),1)));
            
            legend(Deci.Plot.Subtitle{cond})
            title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' ' Deci.Plot.Title{cond} ' Wire'])
            
            if Deci.Plot.Stat.do
                hold on
                
                if barstat{cond}.prob < Deci.Plot.Stat.alpha
                    plot([.75 1.25],[max(ylim) max(ylim)]*.90, '-k', 'LineWidth',2);
                    plot([1],[max(ylim)]*.95, '*k');
                end
            end
        end
        
%         if Deci.Plot.Stat.do
%             if Deci.Plot.Freq.Wire
%                 
%                 set(0, 'CurrentFigure', wire(subj) )
%                 wire(subj).Visible = 'on';
%                 
%                 wiret = repmat(squeeze(wirestat{cond}.prob)' <= .05 ,[ceil(diff(ylim)) 1]);
%                 
%                 h =imagesc(wirestat{cond}.time,min(ylim):max(ylim),wiret);
%                 h.AlphaData = squeeze(wirestat{cond}.prob)' <= .05;
%                 hold on
%                 
%                 set(gca,'YDir','normal');
%                 if Deci.Plot.Stat.do
%                     ButtonH=uicontrol('Parent', wire(subj),'Style','pushbutton','String','p Mask','Position',[10 10 100 25],'Visible','on','Callback',@pmask);
%                     ButtonH.UserData = @zeros;
%                 end
%             end
%         end
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
            suptitle([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' ' Deci.Plot.Title{cond}]);
            
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

    function pmask(PushButton, EventData)
        
        Axes = PushButton.Parent.Children.findobj('Type','Axes');
        Axes = Axes(arrayfun(@(c) ~isempty(c.String), [Axes.Title]));
        
        for a = 1:length(Axes)
            
            imag = Axes(a).Children.findobj('Type','Image');
            
            if isempty(imag)
               imag =  Axes(a).Children.findobj('Type','Surface');
            end
            
            if isempty(imag.UserData)
                imag.UserData = PushButton.UserData(size(imag.AlphaData));
            end
            
            placeholder = imag.UserData;
            imag.UserData = imag.AlphaData;
            imag.AlphaData = placeholder;
        end
    end

end
