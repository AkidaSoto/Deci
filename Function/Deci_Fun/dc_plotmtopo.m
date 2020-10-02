function dc_plottopo(Deci,Subjects,info)


warning('off','MATLAB:handle_graphics:exceptions:SceneNode');

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
        
        tcfg.latency = Deci.Plot.MTopo.Toi;
        
        if info.isfreq
        tcfg.frequency = Deci.Plot.MTopo.Foi;
        end
        
        tcfg.channel = Deci.Plot.MTopo.Channel;
        
        Segdata{subj,conds} = ft_selectdata(tcfg,AvgData{subj,conds});
        
%         if rem(numel(Segdata{subj,conds}.time),Deci.Plot.MTopo.ToiSegs) ~= 0
           Segdata{subj,conds}.powspctrm = permute(mean(reshape(Segdata{subj,conds}.powspctrm(:,:,:,1:end-rem(numel(Segdata{subj,conds}.time),Deci.Plot.MTopo.ToiSegs)),[size(Segdata{subj,conds}.powspctrm,1) size(Segdata{subj,conds}.powspctrm,2) size(Segdata{subj,conds}.powspctrm,3) numel(Segdata{subj,conds}.time(1:end-rem(numel(Segdata{subj,conds}.time),Deci.Plot.MTopo.ToiSegs)))/Deci.Plot.MTopo.ToiSegs Deci.Plot.MTopo.ToiSegs]),4),[1 2 3 5 4]);
           Segdata{subj,conds}.time = mean(reshape(Segdata{subj,conds}.time(1:end-rem(numel(Segdata{subj,conds}.time),Deci.Plot.MTopo.ToiSegs)),[numel(Segdata{subj,conds}.time(1:end-rem(numel(Segdata{subj,conds}.time),Deci.Plot.MTopo.ToiSegs)))/Deci.Plot.MTopo.ToiSegs Deci.Plot.MTopo.ToiSegs]),1);
           
        
%         else
%            Segdata{subj,conds}.powspctrm = permute(mean(reshape(Segdata{subj,conds}.powspctrm(:,:,:,1:end-rem(numel(Segdata{subj,conds}.time),Deci.Plot.MTopo.ToiSegs)),[size(Segdata{subj,conds}.powspctrm,1) size(Segdata{subj,conds}.powspctrm,2) size(Segdata{subj,conds}.powspctrm,3) numel(Segdata{subj,conds}.time(1:end-rem(numel(Segdata{subj,conds}.time),Deci.Plot.MTopo.ToiSegs)))/Deci.Plot.MTopo.ToiSegs Deci.Plot.MTopo.ToiSegs]),4),[1 2 3 5 4]);
%            Segdata{subj,conds}.time = mean(reshape(Segdata{subj,conds}.time,[numel(Segdata{subj,conds}.time)/Deci.Plot.MTopo.ToiSegs Deci.Plot.MTopo.ToiSegs]),2);
%         end
            
        
        if info.isfreq
        tcfg.avgoverfreq = 'yes';
        end

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
    
    
    
    % Two-Way Data Management
    
    if Deci.Plot.Stat.twoway.do
        for cond = 1:length(Deci.Plot.Draw)
            
            if length(Deci.Plot.Draw{cond}) == 4
                
                anovadata = [Segdata{1,Deci.Plot.Draw{cond}}];
                Segdata{1,end+1} = Segdata{1};
                Segdata{1,end}.powspctrm = mean(cat(5,anovadata([3 4]).powspctrm),5) - mean(cat(5,anovadata([1 2]).powspctrm),5);
                Deci.Plot.Draw{end+1} = length(Segdata);
                StatData{end+1} = structfun(@(c) c(:,1,:),StatData{cond},'UniformOutput',false);
                
                Segdata{1,end+1} = Segdata{1};
                Segdata{1,end}.powspctrm = mean(cat(5,anovadata([1 3]).powspctrm),5) - mean(cat(5,anovadata([2 4]).powspctrm),5);    
                Deci.Plot.Draw{end+1} = length(Segdata);
                StatData{end+1} = structfun(@(c) c(:,2,:),StatData{cond},'UniformOutput',false);

                Segdata{1,end+1} = Segdata{1};
                Segdata{1,end}.powspctrm =   [[anovadata(3).powspctrm] - [anovadata(4).powspctrm]]  -  [[anovadata(1).powspctrm] - [anovadata(2).powspctrm]];
                Deci.Plot.Draw{end+1} = length(Segdata);
                StatData{end+1} = structfun(@(c) c(:,3,:),StatData{cond},'UniformOutput',false);
                
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
        
        StatData{cond}.mask  = permute(StatData{cond}.mask,[1 4 3 2]);
        StatData{cond}.stat  = permute(StatData{cond}.stat,[1 4 3 2]);
        StatData{cond}.freq = mean(Segdata{1}.freq);
        StatData{cond}.time = Segdata{1}.time;
        StatData{cond}.label = Segdata{1}.label; 
        StatData{cond}.dimord = 'chan_freq_time';
        
        tacfg.clim = 'maxmin';
        tacfg.maskparameter ='mask';
        tacfg.colormap = Deci.Plot.ColorMap;
        
       
        
        if Deci.Plot.Stat.FPlots
            clear statsub
            
            topot(cond)  = figure;
            topot(cond).Visible = 'on';
     
            for stat = 1:size(StatData{cond}.mask,4)
                
                
                for time = 1:Deci.Plot.MTopo.ToiSegs
                    
                    statsub(stat,time)    =  subplot(size(StatData{cond}.mask,4),Deci.Plot.MTopo.ToiSegs,sub2ind([Deci.Plot.MTopo.ToiSegs size(StatData{cond}.mask,4)],time,stat));
                    tempdata = StatData{cond};
                    tempdata.stat = tempdata.stat(:,:,time,stat);
                    tempdata.mask = tempdata.mask(:,:,time,stat);
                    tempdata.time = Segdata{1}.time(time);
                    
                   
                    
                    tacfg.imagetype = Deci.Plot.ImageType;
                    tacfg.comment = 'no';
                    tacfg.style = 'fill';
                    tacfg.markersymbol = '.';
                    tacfg.colormap = Deci.Plot.ColorMap;
                    
                    if time == Deci.Plot.MTopo.ToiSegs
                    tacfg.colorbar = 'yes';
                    else
                    tacfg.colorbar = 'no';
                    end
                    
                    if time == 1
                        statsub(stat,time).YAxis.Label.Visible = 'on';
                        statsub(stat,time).YAxis.Label.String = ['Test #' num2str(stat)];
                        statsub(stat,time).YAxis.Label.Rotation = 0;
                    end
                    
                    tacfg.contournum = 15;
                    
                    
                    ft_topoplotTFR(tacfg, tempdata);
                    title(num2str(Segdata{1}.time(time)))
                    
                end
               
            end
            dc_pmask(topot(cond))
             

            childs = topot(cond).Children.findobj('Type','Axes');
            for r = 1:length(childs)
                childs(r).CLim = minmax([topot(cond).Children.findobj('Type','Axes').CLim]);
            end
            
           suptitle([Deci.Plot.Stat.Type ' ' Deci.Plot.Title{cond} ' Square (alpha = ' num2str(Deci.Plot.Stat.alpha) ')']);
               
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
            
            for time = 1:Deci.Plot.MTopo.ToiSegs
                
                topfig(subj,subcond,time)    =  subplot(length(Deci.Plot.Draw{cond}),Deci.Plot.MTopo.ToiSegs,sub2ind([Deci.Plot.MTopo.ToiSegs length(Deci.Plot.Draw{cond})],time,subcond));
                
                pcfg = cfg;
                if Deci.Plot.Stat.do
                    pcfg.clim = 'maxmin';
                    pcfg.maskparameter ='mask';
                    Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.mask = StatData{cond}.mask(:,:,time,:); %permute((:,:,time,:),[1 4 3 2]);
                end
                
                scfg.latency = [Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.time(time) Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.time(time)];
                scfg.avgoverfreq = 'yes';
                
                if Deci.Plot.GrandAverage
                    scfg.avgoverrpt = 'yes';
                end
                
                mtopo = ft_selectdata(scfg,Segdata{subj,Deci.Plot.Draw{cond}(subcond)});
                
                
                pcfg.imagetype = Deci.Plot.ImageType;
                pcfg.comment = 'no';
                pcfg.style = 'fill';
                pcfg.markersymbol = '.';
                pcfg.colormap = Deci.Plot.ColorMap;
                
                
                if time == Deci.Plot.MTopo.ToiSegs
                    pcfg.colorbar = 'yes';
                else
                    pcfg.colorbar = 'no';
                end
                
                
                pcfg.contournum = 15;
                ft_topoplotER(pcfg, mtopo);
                
                title([num2str(Segdata{subj,Deci.Plot.Draw{cond}(subcond)}.time(time)) 's']);
                
                if time  ==1
                    topfig(subj,subcond,time).YLabel.Visible = 'on';
                    topfig(subj,subcond,time).YLabel.String = Deci.Plot.Subtitle{cond}{subcond};
                end
                
            end
        end
        
    end
    
    %% Normalize
    
    for subj = 1:size(AvgData,1)
        set(0, 'CurrentFigure', topo(subj) )
        
         childs = topo(subj).Children.findobj('Type','Axes');
        
        for r = 1:length(childs)
            if length(Deci.Plot.Roi) == 2 && isnumeric(Deci.Plot.Roi)
                childs(r).CLim = Deci.Plot.Roi;
            elseif strcmp(Deci.Plot.Roi,'maxmin')
                if ~ isempty(childs(r).Children.UserData)
                    childs(r).CLim = [min(arrayfun(@(c) min(c.Children.UserData(:)),childs(:))) max(arrayfun(@(c) max(c.Children.UserData(:)),childs(:)))];
                else
                    childs(r).CLim= minmax([topo(subj).Children.findobj('Type','Axes').CLim]);
                end
            elseif strcmp(Deci.Plot.Roi,'maxabs')
                if ~isempty(childs(r).Children.findobj('Type','Contour').UserData)
                    childs(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.findobj('Type','Contour').UserData(:))),childs(:))) max(arrayfun(@(c) max(abs(c.Children.findobj('Type','Contour').UserData(:))),childs(:)))];
                else
                    childs(r).CLim = [-1*max(abs(minmax([topo(subj).Children.findobj('Type','Axes').CLim]))) max(abs(minmax([topo(subj).Children.findobj('Type','Axes').CLim])))];
                end
            end
        end
        
         suptitle([Deci.Plot.Title{cond} ' ' Deci.Plot.Freq.Type]);
    end
    
end



