function  dc_plotwire(Deci,Subjects,info)


cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';

%% DataSegmentation

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
        
        if Deci.Plot.GroupLevel
            AvgData2{2,conds} = AvgData{1,conds};
            AvgData2{1,conds} = AvgData{1,conds};
            AvgData2{2,conds}.powspctrm = AvgData{1,conds}.powspctrm(Deci.Plot.Groups{2},:,:,:);
            AvgData2{1,conds}.powspctrm = AvgData{1,conds}.powspctrm(Deci.Plot.Groups{1},:,:,:);
        end
        
    else
        Deci.Plot.Stat.do = false;
        AvgData(:,conds) = Subjects(:,conds);
    end
    
    for subj = 1:size(AvgData(:,conds),1)
        tcfg = [];
        tcfg.nanmean = Deci.Plot.nanmean;
        
        tcfg.latency = Deci.Plot.Wire.Toi;
        tcfg.frequency = Deci.Plot.Wire.Foi;
        tcfg.channel = Deci.Plot.Wire.Channel;
        
        Segdata{subj,conds} = ft_selectdata(tcfg,AvgData{subj,conds});
        
        tcfg.avgoverchan = 'yes';
        tcfg.avgoverfreq = 'yes';
        SegStatdata{subj,conds} = ft_selectdata(tcfg,Segdata{subj,conds});
        
        if strcmpi(Deci.Plot.FreqYScale,'log')
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
    
    if Deci.Plot.GroupLevel
    Deci.SubjectList = {'Group 1' 'Group 2'};    
    
     if Deci.Plot.GroupLevel
        AvgData =  AvgData2;

        for conds = 1:size(Subjects,2)
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
     end
    
    
    else
    Deci.SubjectList = {'Group Average'};
    end
end

for cond = 1:length(Deci.Plot.Draw)
    if Deci.Plot.Stat.do
        tcfg = cfg;
        tcfg.parameter = 'stat';
        StatData{cond}.mask = double(StatData{cond}.mask);
        if Deci.Plot.Stat.FPlots
            wiret(cond)  = figure;
            wiret(cond).Visible = 'on';
            plot(squeeze(Segdata{1}.time),squeeze(StatData{cond}.stat))
            title([Deci.Plot.Stat.Type ' ' Deci.Plot.Title{cond} ' Square (alpha = ' num2str(Deci.Plot.Stat.alpha) ')']);
        end
    end
    
    for subj = 1:size(AvgData,1)
        
        wire(subj)  = figure;
        
        for subcond = 1:length(Deci.Plot.Draw{cond})
            set(0, 'CurrentFigure', wire(subj) )
            wire(subj).Visible = 'on';
            
            top = squeeze(nanmean(nanmean(nanmean(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1),2),3)) + squeeze(nanstd(nanmean(nanmean(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,2),3),[],1))/sqrt(size(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1));
            bot = squeeze(nanmean(nanmean(nanmean(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1),2),3)) - squeeze(nanstd(nanmean(nanmean(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,2),3),[],1))/sqrt(size(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.powspctrm,1));
            
            if Deci.Plot.GrandAverage
                pgon = polyshape([Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.time fliplr(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.time)],[top' fliplr(bot')],'Simplify', false);
                b(subcond) = plot(pgon,'HandleVisibility','off');
                hold on
                b(subcond).EdgeAlpha = 0;
                b(subcond).FaceAlpha = .15;
            end
            
            if Deci.Plot.Stat.do
                Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.mask = squeeze(StatData{cond}.mask);
            end
        end
        
        set(0, 'CurrentFigure', wire(subj) )
        wire(subj).Visible = 'on';
        
        pcfg = cfg;
        pcfg.clim = 'maxmin';
        
        if Deci.Plot.Stat.do
            %pcfg.maskparameter ='mask'
            
            sigs  = Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.mask;
            sigs(sigs == 0) = nan;
            sigs = sigs .* Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.mask;
            
            colors = {'r' 'g' 'b'};
            
            mod = min(cell2mat(arrayfun(@(c) min(mean(c.powspctrm,1),[],'all'),[Segdata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)));
            modfunc = @minus;
            if max(cell2mat(arrayfun(@(c) max(abs(mean(c.powspctrm,1)),[],'all'),[Segdata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false))) < 0
                mod = max(cell2mat(arrayfun(@(c) max(mean(c.powspctrm,1),[],'all'),[Segdata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)));
            modfunc = @plus;
            end
            diffy = diff(ylim);
            hold on
            for z = 1:size(sigs,2)
            plot(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.time,sigs(:,z)*modfunc(mod,diffy*[.01*z]),'HandleVisibility','off','LineWidth',5,'Color',colors{z});
            end
        end
        
        %pcfg.ylim = ylim;
        pcfg.graphcolor = lines;
        pcfg.linewidth = 1;
        pcfg.maskstyle = 'box';
        ft_singleplotER(pcfg,Segdata{subj,Deci.Plot.Draw{cond}});
        
        if Deci.Plot.GrandAverage
            arrayfun(@(c) uistack(c,'top'),b);
            clear b
        end
        
        axis tight
        hold on
        plot([Segdata{cond}.time(1), Segdata{cond}.time(end)], [0 0], 'k--','HandleVisibility','off'); % hor. line
        
        if ~strcmpi(Deci.Plot.BslType,'none')
        plot([0 0], ylim, 'k--','HandleVisibility','off'); % vert. l
        end
        
%         if Deci.Plot.Stat.do
%             boxes = wire(subj).Children(2).Children.findobj('Type','Patch');
%             for bb = 1:length(boxes)
%                 if ~isempty(boxes)
%                     boxes(bb).FaceAlpha = .35;
%                     uistack(boxes(bb),'bottom')
%                     boxes(bb).HandleVisibility = 'off';
%                 end
%             end
%         end
        
        if max(Deci.Plot.Draw{cond}) <= size(info.trllen,2)
            legend(arrayfun(@(a,b) [ Deci.Plot.Freq.Type ' ' a{1} ' (' num2str(b) ')'] ,Deci.Plot.Subtitle{cond},info.trllen(subj,Deci.Plot.Draw{cond}),'UniformOutput',false));
        else
            legend([Deci.Plot.Freq.Type ' '  Deci.Plot.Subtitle{cond}]);
        end
        
        title([Deci.SubjectList{subj} ' ' Deci.Plot.Title{cond}], 'Interpreter', 'none');
        set(legend, 'Interpreter', 'none')
        xlim([Segdata{cond}.time(1) Segdata{cond}.time(end)])
        xlabel('Time');
        
        
        for subcond = 1:length(Deci.Plot.Draw{cond})
            if Deci.Plot.Draw{cond}(subcond) <= size(info.lockers,2) && ~Deci.Plot.GroupLevel
                xlims = xlim;
                ylims = ylim;
                
                for locks = 1:length([info.lockers(subj,Deci.Plot.Draw{cond}(subcond),:)])
                    hold on
                    
                    locktime = [info.lockers(subj,Deci.Plot.Draw{cond}(subcond),locks)/1000];
                    
                    
                    if locktime > xlims(1) && locktime < xlims(2)
                        plotlock = plot([locktime locktime], ylims,'LineWidth',2,'Color','k','LineStyle','--','HandleVisibility','off');
                        plotlock.Color(4) = .2;
                        
                        if Deci.Plot.GrandAverage
                            %this was breaking - JC 5/24/20
                            lockstd = [info.lockersstd(subj,Deci.Plot.Draw{cond}(subcond),locks)/1000];
                            if [locktime - lockstd] < xlims(1)
                                lockpstd(1) = xlims(1);
                            else
                                lockpstd(1) = [locktime - lockstd];
                            end
                            
                            if [locktime + lockstd] > xlims(2)
                                lockpstd(2) = xlims(2);
                            else
                                lockpstd(2) = [locktime + lockstd];
                            end
                            
                            lockpgon = polyshape([lockpstd fliplr(lockpstd)],sort([ylims ylims]),'Simplify', false);
                            lockb = plot(lockpgon,'HandleVisibility','off');
                            hold on
                            lockb.EdgeAlpha = 0;
                            lockb.FaceAlpha = .15;
                            lockb.FaceColor = plotlock.Color;
                            
                            arrayfun(@(c) uistack(c,'top'),lockb);
                            clear lockb
                            
                        end
                    end
                    ylim(ylims)
                end
            end
        end
        
    end
    
end
end