function Plottor(Deci)

if isempty(Deci.Plot.Channel)
    Deci.Plot.Channel = 'all';
    warning('Cannot Find Channels for plot, setting as all');
end

if isempty(Deci.Plot.Toi)
    Deci.Plot.Toi = [-inf inf];
    warning('Cannot Find Toi for plot, setting as [-inf,inf]');
end

Freq    = [];
ERP     = [];

for subject_list = 1:length(Deci.SubjectList)
    if Deci.Plot.Freq.Plot ||  ~isempty(Deci.Plot.PRP)
        if ~isdir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}])
            error(['Freq Analysis not found for '  Deci.SubjectList{subject_list}])
        end
        
        if isempty(Deci.Plot.Freq.Type)
            error('Cannot Find Type for Freq plot');
        end
        
        if isempty(Deci.Plot.Freq.Foi)
            Deci.Plot.Freq.Foi = [-inf inf];
            warning('Cannot Find Foi for Freq plot, setting as [-inf inf]');
        end
        
        if isempty(Deci.Plot.Freq.Roi)
            Deci.Plot.Freq.Roi =  'maxmin';
            warning('Cannot Find Roi for Freq plot, setting as maxmin');
        end
        
        if isempty(Deci.Plot.Freq.BslType)
            Deci.Plot.Freq.BslType =  'maxmin';
            warning('Cannot Find BslType for Freq plot, setting as relative');
        end
        
        if isempty(Deci.Plot.Freq.Bsl)
            error('Cannot have Bsl empty for plot');
        end
        
        Freq.Conditions     = CleanDir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}]);
        Freq.Channels       = CleanDir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list} filesep Freq.Conditions{1}]);
        Freq.Channels       = cellfun(@(c) c(1:end-4),Freq.Channels,'un',0);
    end
    
    if Deci.Plot.ERP || ~isempty(Deci.Plot.PRP)
        if ~isdir([Deci.Folder.Analysis filesep 'Volt_ERP' filesep Deci.SubjectList{subject_list}])
            error(['ERP Analysis not found for '  Deci.SubjectList{subject_list}])
        end
        
        ERP.Conditions     = CleanDir([Deci.Folder.Analysis filesep 'Volt_ERP' filesep Deci.SubjectList{subject_list}]);
    end
    
    if ~isempty(Deci.Plot.PRP)
        
        if isempty(Deci.Plot.PRP.label)
            Deci.Plot.Freq.BslType =  0;
            warning('Cannot Find label for PRP plot, setting as 0');
        end
        
        if isempty(Deci.Plot.PRP.ScatterToi)
            Deci.Plot.Freq.ScatterToi =  [-inf inf];
            warning('Cannot Find ScatterToi for PRPplot, setting as  [-inf inf]');
        end
        
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
                toi = round(freq.time,4) >= Deci.Plot.Toi(1) & round(freq.time,4) <= Deci.Plot.Toi(2);
                
                freq.freq =  freq.freq(foi);
                freq.time =  freq.time(toi);
                freq.powspctrm  =freq.powspctrm(:,foi,toi);
                freq.label = Freq.Channels(Channel);
                
                Chans{Channel} = freq;
            end
            
            acfg.parameter = 'powspctrm';
            acfg.appenddim = 'chan';
            Subjects{subject_list} = rmfield(ft_appendfreq(acfg,Chans{:}),'cfg');
            clear Chans;
            
            acfg.baseline = Deci.Plot.Freq.Bsl;
            acfg.baselinetype = Deci.Plot.Freq.BslType;
            Subjects{subject_list} = rmfield(ft_freqbaseline(acfg,Subjects{subject_list}),'cfg');
            
        end
        
        if Deci.Plot.GA
            facfg.parameter =  'powspctrm';
            FreqData{Condition} = rmfield(ft_freqgrandaverage(facfg,Subjects{:}),'cfg');
            
        else
            FreqData(Condition,:) = Subjects(:);
        end
        clear Subjects;
    end
end

if Deci.Plot.ERP || ~isempty(Deci.Plot.PRP)
    for  Condition =  1:length(ERP.Conditions)
        for  subject_list = 1:length(Deci.SubjectList)
            
            time = [];
            load([Deci.Folder.Analysis filesep 'Volt_ERP' filesep Deci.SubjectList{subject_list} filesep ERP.Conditions{Condition}],'time');
            
            toi = round(time.time,4) >= Deci.Plot.Toi(1) & round(time.time,4) <= Deci.Plot.Toi(2);
            time.time =  time.time(toi);
            time.avg  =time.avg(:,toi);
            
            Subjects{subject_list} = time;
            clear time;
        end
        
        if Deci.Plot.GA
            facfg.parameter =  'avg';
            TimeData{Condition} = rmfield(ft_timelockgrandaverage(facfg,Subjects{:}),'cfg');
            
        else
            TimeData(Condition,:) = Subjects(:);
        end
        clear Subjects;
    end
end


%% Plot

cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';

if Deci.Plot.Freq.Plot
    
    for subj = 1:size(FreqData,2)
        square = figure;
        topo = figure;
        
        for cond = 1:size(FreqData,1)
            
            set(0, 'CurrentFigure', square)
            subby(subj,cond) = subplot(size(FreqData,1),1,cond);
            ft_singleplotTFR(cfg,FreqData{cond,subj})
            
            title(['Subj ' num2str(subj) ' ' Deci.Plot.Freq.Type ' Cond'  num2str(cond)]);
            colorbar('vert');
            map = colormap('jet'); %'hot' 'gray'
            colormap(map);
            
            if length(Freq.Channels) ~= 1
                
                set(0, 'CurrentFigure', topo)
                cirky(subj,cond)    =  subplot(size(FreqData,1),1,cond);
                ft_topoplotER(cfg, FreqData{cond,subj});
                
                hold on
                a = get(gca,'Children');
                conty = a(arrayfun(@(c) strcmp(c.Type,'contour'),a));
                surfy  = a(arrayfun(@(c) strcmp(c.Type,'surface'),a));
                
                spot = find(arrayfun(@(c) strcmp(c.Type,'contour'),a));
                
                actmin = min(min(surfy.CData(~isnan(surfy.CData))));
                actmax = max(max(surfy.CData(~isnan(surfy.CData))));
                
                h = figure;
                [~,cont] = contourf(gca,surfy.CData,linspace(actmin,actmax,6));
                a(spot).LevelList = cont.LevelList;
                a(spot).TextList = cont.TextList;
                a(spot).Fill = 'on';
                delete(h);
                
                axis xy
                title(['Subj ' num2str(subj) ' ' Deci.Plot.Freq.Type ' Cond'  num2str(cond)]);
                colorbar('vert');
                map = colormap('jet'); %'hot' 'gray'
                colormap(map);
            else
                close topo
            end
        end
    end
    
    if length(Freq.Channels) ~= 1
        set(0, 'CurrentFigure', topo)
        for r = 1:length(cirky(:))
            cirky(r).CLim = [min([cirky.CLim]) max([cirky.CLim])];
        end
    end
    
    set(0, 'CurrentFigure', square)
    for r = 1:length(subby(:))
        
        if length(Deci.Plot.Freq.Roi) == 2 && isnumeric(Deci.Plot.Freq.Roi)
            subby(r).CLim = Deci.Plot.Freq.Roi;
        elseif strcmp(Deci.Plot.Freq.Roi,'maxmin')
            subby(r).CLim = [min([subby.CLim]) max([subby.CLim])];
        elseif strcmp(Deci.Plot.Freq.Roi,'maxabs')
            subby(r).CLim = [-1*max(abs([subby.CLim])) max(abs([subby.CLim]))];
        end
    end
    
end

if Deci.Plot.ERP
    %Not yet available
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
                set(t, 'Clipping', 'on');
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
        
        scfg.parameter = 'avg';
        SFreq = rmfield(ft_timelockgrandaverage(scfg,FreqData{cond,:}),'cfg');
        STime = rmfield(ft_timelockgrandaverage(scfg,TimeData{cond,:}),'cfg');
        
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
    title(ax3,{['Between Trial Correlation with p and rho'] ['at times ' scattitle]} );
    
    
    xlabel(ax3,'Time(secs)');
    ylabel(ax3, 'Rho' );
    title(ax3,{['Correlation by time  only sigficiant by Fdr MC Corrected'] ['at times ' scattitle]} )
    
    xlabel(ax4,'Total Power');
    ylabel(ax4, 'ERP Power (uV)' );
    
    title(ax4,{['Correlation Total Power-ERP by subject'] ['at times ' scattitle]} )
    
    
end

end
