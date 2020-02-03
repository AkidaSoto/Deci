function dc_plotsquare(Deci,Subjects,info)


cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';

%% Data Segmentation

for conds = 1:size(Subjects,2)
    
    if Deci.Plot.GrandAverage
        facfg.parameter =  info.parameter;
        facfg.type = 'mean';
        facfg.keepindividual = 'yes';
        
        if info.isfreq
            evalc('AvgData{1,conds} = ft_freqgrandaverage(facfg,Subjects{:,conds});');
        else
            evalc('AvgData{1,conds} = ft_timelockgrandaverage(facfg,Subjects{:,conds});');
        end
        AvgData{1,conds} = rmfield(AvgData{1,conds},'cfg');
        
    else
        Deci.Plot.Stat.do = false;
        AvgData(:,conds) = Subjects(:,conds);
    end
    
    for subj = 1:size(AvgData(:,conds),1)
        tcfg = [];
        tcfg.nanmean = Deci.Plot.nanmean;
        
        tcfg.latency = Deci.Plot.Square.Toi;
        tcfg.frequency = Deci.Plot.Square.Foi;
        tcfg.channel = Deci.Plot.Square.Channel;
        
        Segdata{subj,conds} = ft_selectdata(tcfg,AvgData{subj,conds});
        
        tcfg.avgoverchan = 'yes';
        SegStatdata{subj,conds} = ft_selectdata(tcfg,Segdata{subj,conds});
        
        if strcmpi(Deci.Plot.FreqYScale,'log')
            Segdata{subj,conds}.freq = log(Segdata{subj,conds}.freq);
            SegStatdata{subj,conds}.freq = log(SegStatdata{subj,conds}.freq);
        end
    end
end

%% Stat

if Deci.Plot.Stat.do
     StatData = dc_plotstat(Deci,SegStatdata,info);
end

%% Plot

if Deci.Plot.GrandAverage
    Deci.SubjectList = {'Group Average'};
end

for cond = 1:length(Deci.Plot.Draw)
    if Deci.Plot.Stat.do
        tacfg = cfg;
        tacfg.parameter = 'stat';
        StatData{cond}.mask = double(StatData{cond}.mask);
        StatData{cond}.mask(StatData{cond}.mask == 0) = .2;
        tacfg.maskparameter = 'mask';
        tacfg.colormap = Deci.Plot.ColorMap;
        
        if Deci.Plot.Stat.FPlots
            squaret(cond) = figure;
            squaret(cond).Visible = 'on';
            
            ft_singleplotTFR(tacfg,StatData{cond})
            title([Deci.Plot.Stat.Type ' ' Deci.Plot.Title{cond} ' Square (alpha = ' num2str(Deci.Plot.Stat.alpha) ')']);
            dc_pmask(squaret(cond))
            
            ylabel('F score');
            xlabel('time');
        end
    end
    
    for subj = 1:size(AvgData,1)
        
 
        square(subj) = figure;
        
        if Deci.Plot.Stat.do
            dc_pmask(square)
        end

        for subcond = 1:length(Deci.Plot.Draw{cond})
            
            set(0, 'CurrentFigure', square(subj) )
            square(subj).Visible = 'on';
            subby(subj,subcond) = subplot(length(Deci.Plot.Draw{cond}),1,subcond );
            
            pcfg = cfg;
            if Deci.Plot.Stat.do
                pcfg.clim = 'maxmin';
                pcfg.maskparameter ='mask';
                Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.mask = repmat(StatData{cond}.mask,[length(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.label) 1 1]);
            end
            
            pcfg.imagetype = 'contourf';
            pcfg.colormap = Deci.Plot.ColorMap;
            evalc('ft_singleplotTFR(pcfg,Segdata{subj,Deci.Plot.Draw{cond}(subcond)})');
            axis tight
            
            
            if Deci.Plot.Draw{cond}(subcond) <= size(info.lockers,2)
                xlims = xlim;
                ylims = ylim;
                
                for locks = 1:length([info.lockers(subj,Deci.Plot.Draw{cond}(subcond),:)])
                    hold on
                    
                    locktime = [info.lockers(subj,Deci.Plot.Draw{cond}(subcond),locks)/1000];
                    

                    if Deci.Plot.GrandAverage
                        if locktime > xlims(1) && locktime < xlims(2)
                            info.lockstd = [info.lockersstd(subj,Deci.Plot.Draw{cond}(subcond),locks)/1000];
                            plotlock = line([locktime locktime], ylims,'LineWidth',2,'Color','k','LineStyle','--','HandleVisibility','off');
                            
                            
                            if [locktime - info.lockstd] < xlims(1)
                                lockpstd(1) = xlims(1);
                            else
                                lockpstd(1) = [locktime - info.lockstd];
                            end
                            
                            if [locktime + info.lockstd] > xlims(2)
                                lockpstd(2) = xlims(2);
                            else
                                lockpstd(2) = [locktime + info.lockstd];
                            end
                            
                            lockpgon = polyshape([lockpstd fliplr(lockpstd)],sort([ylims ylims]),'Simplify', false);
                            lockb = plot(lockpgon,'HandleVisibility','off');
                            hold on
                            lockb.EdgeAlpha = 0;
                            lockb.FaceAlpha = .30;
                            lockb.FaceColor = plotlock.Color;
                            
                            arrayfun(@(c) uistack(c,'top'),lockb);
                            clear lockb
                            
                        end
                        
                    else
                        if locktime > xlims(1) && locktime < xlims(2)
                            plotlock = line([locktime locktime], ylims,'LineWidth',2,'Color','k','LineStyle','--','HandleVisibility','off');
                        end
                    end
                    
                end
                ylim(ylims)
                title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' '  Deci.Plot.Subtitle{cond}{subcond} ' (' num2str(info.trllen(subj,Deci.Plot.Draw{cond}(subcond))) ')']);
            else
                title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' '  Deci.Plot.Subtitle{cond}{subcond}]);
            end
        end
    end

    for subj = 1:size(AvgData,1)
        set(0, 'CurrentFigure', square(subj))
        suptitle([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' ' Deci.Plot.Title{cond}]);
        for r = 1:length(subby(:))
            if length(Deci.Plot.Roi) == 2 && isnumeric(Deci.Plot.Roi)
                subby(r).CLim = Deci.Plot.Roi;
            elseif strcmp(Deci.Plot.Roi,'maxmin')
                if ~ isempty(subby(r).Children.UserData)
                    subby(r).CLim = [min(arrayfun(@(c) min(c.Children.UserData(:)),subby(:))) max(arrayfun(@(c) max(c.Children.UserData(:)),subby(:)))];
                else
                    subby(r).CLim = [min(arrayfun(@(c) min(c.Children.ZData(:)),subby(:))) max(arrayfun(@(c) max(c.Children.ZData(:)),subby(:)))];
                end
                
                %subby(r).Children.LevelList = [linspace(subby(r).CLim(1),subby(r).CLim(2),12)];
                
            elseif strcmp(Deci.Plot.Roi,'maxabs')
                if ~ isempty(subby(r).Children.UserData)
                    subby(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.UserData(:))),subby(:))) max(arrayfun(@(c) max(abs(c.Children.UserData(:))),subby(:)))];
                else
                    subby(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.ZData(:))),subby(:))) max(arrayfun(@(c) max(abs(c.Children.ZData(:))),subby(:)))];
                    
                end
                %(r).Children.LevelListMode = 'auto';
                
                %subby(r).Children.LevelList = [linspace(subby(r).CLim(1),subby(r).CLim(2),12)];
                
            end
        end
        
        if strcmpi(Deci.Plot.FreqYScale,'log')
            for r = 1:length(subby(:))
                subby(r).YTickLabels = round(exp(subby(r).YTick),1);
            end
        end
        
        if ~isempty(Deci.Folder.Plot)
            %mkdir([Deci.Folder.Plot filesep Deci.Plot.Title{cond}]);
            %saveas(square(subj),[Deci.Folder.Plot filesep Deci.Plot.Title{cond} filesep Deci.SubjectList{subj} '_square'],Deci.Plot.Save.Format);
        end
    end
    
end
