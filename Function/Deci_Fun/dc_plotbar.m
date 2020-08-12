function dc_plotbar(Deci,Subjects,info)


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
        
        tcfg.latency = Deci.Plot.Bar.Toi;
        tcfg.frequency = Deci.Plot.Bar.Foi;
        tcfg.channel = Deci.Plot.Bar.Channel;
        
        Segdata{subj,conds} = ft_selectdata(tcfg,AvgData{subj,conds});
        
        tcfg.avgoverchan = 'yes';
        tcfg.avgoverfreq = 'yes';
        tcfg.avgovertime = 'yes';
        
        SegStatdata{subj,conds} = ft_selectdata(tcfg,Segdata{subj,conds});
        
        if strcmpi(Deci.Plot.FreqYScale,'log')
            Segdata{subj,conds}.freq = log(Segdata{subj,conds}.freq);
            SegStatdata{subj,conds}.freq = log(SegStatdata{subj,conds}.freq);
        end
        
    end
    
end

%% Stats

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
        StatData{cond}.mask= double(StatData{cond}.mask);
        StatData{cond}.mask(StatData{cond}.mask == 0) = nan;
        
        if Deci.Plot.Stat.FPlots
            bart(cond)  = figure;
            bart(cond).Visible = 'on';
            bar(diag(squeeze(StatData{cond}.stat)),'stacked')
            hold on
            plot(1:length(squeeze(StatData{cond}.stat)),squeeze(StatData{cond}.mask).*squeeze(StatData{cond}.stat) + 1,'*k','HandleVisibility','off')
            title([Deci.Plot.Stat.Type ' ' Deci.Plot.Title{cond} ' Square (alpha = ' num2str(Deci.Plot.Stat.alpha) ')']);
        end
        
    end
    
    for subj = 1:size(AvgData,1)
        
        
        barfig(subj)  = figure;
        set(0, 'CurrentFigure', barfig(subj) )
        barfig(subj).Visible = 'on';
        
        CleanBars(mean(cell2mat(arrayfun(@(c) nanmean(nanmean(nanmean(c.powspctrm,2),3),4),[Segdata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)),1), ...
            nanstd(cell2mat(arrayfun(@(c) nanmean(nanmean(nanmean(c.powspctrm,2),3),4),[Segdata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)),[],1) ...
            /sqrt(size(cell2mat(arrayfun(@(c) nanmean(nanmean(nanmean(c.powspctrm,2),3),4),[Segdata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)),1)));
        
        l = legend(Deci.Plot.Subtitle{cond});
        title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' ' Deci.Plot.Title{cond}])
        set(l, 'Interpreter', 'none')
        if Deci.Plot.Stat.do
            hold on
            
            if isfield(StatData{cond},'prob')
                if StatData{cond}.prob < Deci.Plot.Stat.alpha
                    
                    if max(ylim) ~= 0
                    plot([.75 1.25],[max(ylim) max(ylim)]*.90, '-k', 'LineWidth',2,'HandleVisibility','off');
                     plot([1],[max(ylim)]*.95, '*k','HandleVisibility','off');
                    else
                    plot([.75 1.25],[min(ylim) min(ylim)]*.90, '-k', 'LineWidth',2,'HandleVisibility','off');
                     plot([1],[min(ylim)]*.95, '*k','HandleVisibility','off');
                    end
                   
                end
            end
        end
        
    end
    
end


end