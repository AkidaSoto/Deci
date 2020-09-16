function dc_plotsquare(Deci,Subjects,info)


cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';

%% Data Segmentation

if Deci.Plot.GrandAverage
    if Deci.Plot.GroupLevel
        
        Deci.Plot = Exist(Deci.Plot,'Groups',[]);
        
        if isempty(Deci.Plot.Groups)
            fakeUI = figure;
            fakeUI.UserData = Deci.SubjectList;
            fakeUI.Visible =  'off';
            dc_select_labels(fakeUI,[],Deci.SubjectList,{'Group 1','Group 2'});
            waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
            Deci.Plot.Groups{2} = ismember(Deci.SubjectList,fakeUI.UserData);
            Deci.Plot.Groups{1} = ~ismember(Deci.SubjectList,fakeUI.UserData);
            close(fakeUI);
        end
    end
end


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
     
     
         % Two-Way Data Management
    
    if Deci.Plot.Stat.twoway.do
        for cond = 1:length(Deci.Plot.Draw)
            
            if length(Deci.Plot.Draw{cond}) == 4 && ~strcmpi(Deci.Plot.Stat.Comp,'Bsl')
                
                anovadata = [Segdata{1,Deci.Plot.Draw{cond}}];
                Segdata{1,end+1} = Segdata{1};
                Segdata{1,end}.powspctrm = mean(cat(5,anovadata(3:4).powspctrm),5) - mean(cat(5,anovadata(1:2).powspctrm),5);
                Deci.Plot.Draw{end+1} = length(Segdata);
                StatData{end+1} = structfun(@(c) c(:,:,:,1),StatData{cond},'UniformOutput',false);
                
                Segdata{1,end+1} = Segdata{1};
                Segdata{1,end}.powspctrm = mean(cat(5,anovadata(2:4).powspctrm),5) - mean(cat(5,anovadata([1 3]).powspctrm),5);    
                Deci.Plot.Draw{end+1} = length(Segdata);
                StatData{end+1} = structfun(@(c) c(:,:,:,2),StatData{cond},'UniformOutput',false);

                Segdata{1,end+1} = Segdata{1};
                Segdata{1,end}.powspctrm =   [[anovadata(4).powspctrm] - [anovadata(3).powspctrm]]  -  [[anovadata(2).powspctrm] - [anovadata(1).powspctrm]];
                Deci.Plot.Draw{end+1} = length(Segdata);
                StatData{end+1} = structfun(@(c) c(:,:,:,3),StatData{cond},'UniformOutput',false);
                
                Deci.Plot.Title(end+1:end+3)        = Deci.Plot.Stat.twoway.Title;
                Deci.Plot.Subtitle(end+1:end+3)   = Deci.Plot.Stat.twoway.Subtitle;
                
            end
        end
        
        
    end
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
        tacfg = cfg;
        tacfg.parameter = 'stat';
        StatData{cond}.mask = double(StatData{cond}.mask);
        StatData{cond}.mask(StatData{cond}.mask == 0) = .2;
        StatData{cond}.freq = Segdata{1}.freq;
        StatData{cond}.time = Segdata{1}.time;
        StatData{cond}.label = {'all'}; 
        StatData{cond}.dimord = 'chan_freq_time';
        
        tacfg.colormap = Deci.Plot.ColorMap;
        
        if Deci.Plot.Stat.FPlots
            squaret(cond) = figure;
            squaret(cond).Position = Deci.Plot.Size;
            squaret(cond).Visible = 'on';
            
            tacfg.imagetype = Deci.Plot.ImageType;
            tacfg.colormap = Deci.Plot.ColorMap;
            tacfg.clim = 'maxmin';
            tacfg.colorbar = 'yes';
            tacfg.maskparameter = 'mask';
            
            
            for stat = 1:size(StatData{cond}.mask,4)
                
                statsub(stat)    =  subplot(size(StatData{cond}.mask,4),1,sub2ind([1 size(StatData{cond}.mask,4)],1,stat));
                tempdata = StatData{cond};
                tempdata.stat = tempdata.stat(:,:,:,stat);
                tempdata.mask = tempdata.mask(:,:,:,stat);
                tempdata.time = Segdata{1}.time;
                
                tacfg.colorbar = 'yes';
                
                statsub(stat).YAxis.Label.Visible = 'on';
                statsub(stat).YAxis.Label.String = ['Test #' num2str(stat)];
                statsub(stat).YAxis.Label.Rotation = 0;
                
                tacfg.contournum = 15;
                
                
                ft_singleplotTFR(tacfg,tempdata);
                %squaret(cond).SizeChangedFcn = {@(m,c) set(m,'Position',c.Position),m,c)
                
                ylabel('Freq');
                xlabel('time');
                title([Deci.Plot.Stat.Type ' ' Deci.Plot.Title{cond} ' Square (alpha = ' num2str(Deci.Plot.Stat.alpha) ')']);
            end

            dc_pmask(squaret(cond))
            
            childs =  squaret(cond).Children.findobj('Type','Axes');
            for r = 1:length(childs)
                childs(r).CLim = minmax([squaret(cond).Children.findobj('Type','Axes').CLim]);
                
                if strcmpi(Deci.Plot.FreqYScale,'log')
                childs(r).YTickLabels = round(exp(childs(r).YTick),1);
                end
            end

            
        end
    end
    
    for subj = 1:size(AvgData,1)
        
 
        square(subj) = figure;
        
        if Deci.Plot.Stat.do
            dc_pmask(square(subj))
        end

        for subcond = 1:length(Deci.Plot.Draw{cond})
            
            set(0, 'CurrentFigure', square(subj) )
            square(subj).Visible = 'on';
            subby(subj,subcond) = subplot(length(Deci.Plot.Draw{cond}),1,subcond );
            
            pcfg = cfg;
            if Deci.Plot.Stat.do
                pcfg.clim = 'maxmin';
                pcfg.maskparameter ='mask';
                Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.mask = StatData{cond}.mask; %repmat(,[length(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.label) 1 1]);
            end
            
            scfg.avgoverchan = 'yes';
            
            if Deci.Plot.GrandAverage
            scfg.avgoverrpt = 'yes';
            end
            
            tsquare = ft_selectdata(scfg,Segdata{subj,Deci.Plot.Draw{cond}(subcond)});
            
            pcfg.imagetype = Deci.Plot.ImageType;
            pcfg.colormap = Deci.Plot.ColorMap;
            evalc('ft_singleplotTFR(pcfg,tsquare)');
            
            
            
            if Deci.Plot.Draw{cond}(subcond) <= size(info.lockers,2) && ~Deci.Plot.GroupLevel
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
        
        childs = square(subj).Children.findobj('Type','Axes');
        
        k = [square.Children];
       
        for r = 1:length(childs)
            if length(Deci.Plot.Roi) == 2 && isnumeric(Deci.Plot.Roi)
                childs(r).CLim = Deci.Plot.Roi;
            elseif strcmp(Deci.Plot.Roi,'maxmin')
                    childs(r).CLim = minmax([k.findobj('Type','Axes').CLim]);
            elseif strcmp(Deci.Plot.Roi,'maxabs')
                    childs(r).CLim = [-1*max(abs(minmax([k.findobj('Type','Axes').CLim]))) max(abs(minmax([k.findobj('Type','Axes').CLim])))];
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
    
    for subj = 1:size(AvgData,1)
         set(0, 'CurrentFigure', square(subj))
       suptitle([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' ' Deci.Plot.Title{cond}]); 
    end
    
end
