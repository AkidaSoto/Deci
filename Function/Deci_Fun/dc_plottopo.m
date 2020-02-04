function dc_plottopo(Deci,Subjects,info)


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
        
        tcfg.latency = Deci.Plot.Topo.Toi;
        
        if info.isfreq
        tcfg.frequency = Deci.Plot.Topo.Foi;
        end
        
        tcfg.channel = Deci.Plot.Topo.Channel;
        
        Segdata{subj,conds} = ft_selectdata(tcfg,AvgData{subj,conds});
        
        if info.isfreq
        tcfg.avgoverfreq = 'yes';
        end
        
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
    Deci.SubjectList = {'Group Average'};
end

for cond = 1:length(Deci.Plot.Draw)
    if Deci.Plot.Stat.do
        tacfg = cfg;
        tacfg.parameter = 'stat';
        StatData{cond}.mask = double(StatData{cond}.mask);
        StatData{cond}.mask(StatData{cond}.mask == 0) = .2;
        tacfg.clim = 'maxmin';
        tacfg.maskparameter ='mask';
        tacfg.colormap = Deci.Plot.ColorMap;
        
        if Deci.Plot.Stat.FPlots
            topot(cond)  = figure;
            topot(cond).Visible = 'on';
            
            ft_topoplotTFR(tcfg, StatData{cond});
            title([Deci.Plot.Stat.Type ' ' Deci.Plot.Title{cond} ' Square (alpha = ' num2str(Deci.Plot.Stat.alpha) ')']);
            dc_pmask(topot(cond))
        end
    end
    
    clear topfig
    for subj = 1:size(AvgData,1)
        
        topo(subj)  = figure;
        
        if Deci.Plot.Stat.do
            dc_pmask(topo)
        end
        
        for subcond = 1:length(Deci.Plot.Draw{cond})
            
            set(0, 'CurrentFigure', topo(subj))
            topo(subj).Visible = 'on';
            topfig(subj,subcond)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
            
            pcfg = cfg;
            if Deci.Plot.Stat.do
                pcfg.clim = 'maxmin';
                pcfg.maskparameter ='mask';
                Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.mask = repmat(StatData{cond}.mask',[size(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.(info.parameter),1) 1 length(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.freq) length(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.time)]);
            end
            
            pcfg.imagetype = 'contourf';
            pcfg.comment = 'no';
            pcfg.style = 'fill';
            pcfg.markersymbol = '.';
            pcfg.colormap = Deci.Plot.ColorMap;
            pcfg.colorbar = 'yes';
            pcfg.contournum = 15;
            ft_topoplotER(pcfg, Segdata{subj,Deci.Plot.Draw{cond}(subcond)});
            
            if Deci.Plot.Draw{cond}(subcond) <= size(info.trllen,2)
                title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' '  Deci.Plot.Subtitle{cond}{subcond} ' (' num2str(info.trllen(subj,Deci.Plot.Draw{cond}(subcond))) ')']);
                
            else
                title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' '  Deci.Plot.Subtitle{cond}{subcond}]);
            end
            
        end
        
    end
    
    %% Normalize
    
    for subj = 1:size(AvgData,1)
        set(0, 'CurrentFigure', topo(subj) )
        suptitle(Deci.Plot.Title{cond});
        for r = 1:length(topfig(:))
            if length(Deci.Plot.Roi) == 2 && isnumeric(Deci.Plot.Roi)
                topfig(r).CLim = Deci.Plot.Roi;
            elseif strcmp(Deci.Plot.Roi,'maxmin')
                if ~ isempty(topfig(r).Children.UserData)
                    topfig(r).CLim = [min(arrayfun(@(c) min(c.Children.UserData(:)),topfig(:))) max(arrayfun(@(c) max(c.Children.UserData(:)),topfig(:)))];
                else
                    topfig(r).CLim= [min(arrayfun(@(c) min(c.Children.ZData(:)),topfig(:))) max(arrayfun(@(c) max(c.Children.ZData(:)),topfig(:)))];
                end
            elseif strcmp(Deci.Plot.Roi,'maxabs')
                if ~isempty(topfig(r).Children.findobj('Type','Contour').UserData)
                    topfig(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.findobj('Type','Contour').UserData(:))),topfig(:))) max(arrayfun(@(c) max(abs(c.Children.findobj('Type','Contour').UserData(:))),topfig(:)))];
                else
                    topfig(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.findobj('Type','Contour').ZData(:))),topfig(:))) max(arrayfun(@(c) max(abs(c.Children.findobj('Type','Contour').ZData(:))),topfig(:)))];
                end
            end
        end
    end
    
end


