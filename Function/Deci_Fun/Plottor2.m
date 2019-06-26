function Plottor2(Deci)


%% File Checks

for subject_list = 1:length(Deci.SubjectList)
    
    if Deci.Run.Freq || Deci.Run.CFC
        if ~isempty(Deci.Plot.Freq) ||  ~isempty(Deci.Plot.PRP) ||  ~isempty(Deci.Plot.CFC)
            if ~isdir([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep Deci.SubjectList{subject_list}])
                error(['Freq Analysis not found for '  Deci.SubjectList{subject_list}])
            end
        end
    end
    
    if Deci.Run.ERP
        if ~isempty(Deci.Plot.ERP) || ~isempty(Deci.Plot.PRP)
            if ~isdir([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list}])
                error(['ERP Analysis not found for '  Deci.SubjectList{subject_list}])
            end
        end
    end
    
    if Deci.Run.Behavior
        if ~isempty(Deci.Plot.Behv)
            
            switch Deci.Plot.Behv.Source
                case 'Definition'
                    if exist([Deci.Folder.Version filesep 'Definition' filesep Deci.SubjectList{subject_list} '.mat'],'file') ~= 2
                        error(['Definition (Behv) Analysis not found for '  Deci.SubjectList{subject_list}])
                    end
                case 'Freq'
                    if ~isempty(Deci.Plot.Freq) ||  ~isempty(Deci.Plot.PRP) ||  ~isempty(Deci.Plot.CFC)
                        if ~isdir([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep Deci.SubjectList{subject_list}])
                            error(['Freq Analysis not found for '  Deci.SubjectList{subject_list}])
                        end
                    end
            end
            
        end
    end
end

if ~isempty(Deci.Folder.Plot)
    mkdir(Deci.Folder.Plot);
end


%% Split

if Deci.Run.Freq && ~isempty(Deci.Plot.Freq)
    Plottor_Freq2(Deci);
end

if Deci.Run.ERP && ~isempty(Deci.Plot.ERP)
    Plottor_ERP(Deci);
end

% if Deci.Run.PRP && ~isempty(Deci.Plot.PRP)
%
% end

if Deci.Run.CFC && ~isempty(Deci.Plot.CFC)
    Plottor_CFC2(Deci);
end

if Deci.Run.Behavior && ~isempty(Deci.Plot.Behv)
    Plottor_Behv(Deci);
end

if Deci.Run.Extra && ~isempty(Deci.Plot.Extra)
    for fun = find(Deci.Plot.Extra.List)
        
        Deci.Plot.Extra.Functions{fun}(Deci.Plot.Extra.Params{fun})
    end
end

%% Plot

if ~isempty(Deci.Plot.PRP)
    
    
    PRPPlot = figure;
    colors = reshape(hsv(size(FreqData,2)),[3 size(FreqData,2) 1]);
    
    for subj = 1:size(FreqData,1)
        for cond = 1:size(FreqData,2)
            
            %Collapse Freq and Channel
            FreqData{subj,cond}.avg = permute(nanmean(nanmean(FreqData{subj,cond}.powspctrm,2),1),[1 3 2]);
            FreqData{subj,cond} = rmfield(FreqData{subj,cond},'powspctrm');
            FreqData{subj,cond}.label = {'mix'};
            FreqData{subj,cond}.dimord = 'chan_time';
            Hz = regexprep(num2str(num2str(minmax(FreqData{subj,cond}.freq)),'%.3g '),' +','-');
            FreqData{subj,cond} = rmfield(FreqData{subj,cond},'freq');
            
            %Collapse Channel
            TimeData{subj,cond}.avg = mean(TimeData{subj,cond}.avg,1);
            TimeData{subj,cond}.label = {'mix'};
            TimeData{subj,cond}.dimord = 'chan_time';
        end
    end
    
    
    
    for cond = 1:size(FreqData,2)
        
        clear SF ST
        
        for  subj = 1:size(FreqData,1)
            
            for k = 1:length(TimeData{subj,cond}.time)
                SF(subj,k) = FreqData{subj,cond}.avg(:,k);
                ST(subj,k) = TimeData{subj,cond}.avg(:,k);
            end
            
            ax4 = subplot(2,2,4);
            ax4.XAxisLocation = 'origin';
            ax4.YAxisLocation = 'origin';
            hold(ax4,'on');
            scattoi = TimeData{subj,cond}.time >= Deci.Plot.PRP.ScatterToi(1) &  TimeData{subj,cond}.time <= Deci.Plot.PRP.ScatterToi(2);
            scat(subj,cond,1) = nanmean(FreqData{subj,cond}.avg(:,scattoi),2);
            scat(subj,cond,2) = nanmean(TimeData{subj,cond}.avg(:,scattoi),2);
            scattitle = regexprep(num2str(num2str(minmax(TimeData{subj,cond}.time(scattoi))),'%.3g '),' +','-');
            
            scatter(ax4,scat(subj,cond,1),scat(subj,cond,2),'MarkerFaceColor',colors(:,cond))
            
            if Deci.Plot.PRP.label
                t =  text(ax4, scat(subj,cond,1),scat(subj,cond,2),Deci.SubjectList{subj});
                set(t, 'Clipping', 'on','Interpreter', 'none');
            end
            
            
        end
        
        for k = 1:length(TimeData{subj,cond}.time)
            RHO(k) = corr(SF(:,k),ST(:,k));
        end
        
        ax3 = subplot(2,2,3);
        ax3.XAxisLocation = 'origin';
        ax3.YAxisLocation = 'origin';
        RHOim = plot(ax3,TimeData{subj,cond}.time,RHO,'Color',colors(:,cond));
        hold(ax3,'on');
        
        tcfg.parameter = 'avg';
        SFreq = rmfield(ft_timelockgrandaverage(tcfg,FreqData{cond,:}),'cfg');
        STime = rmfield(ft_timelockgrandaverage(tcfg,TimeData{cond,:}),'cfg');
        
        ax1 = subplot(2,2,1);
        ax1.XAxisLocation = 'origin';
        ax1.YAxisLocation = 'origin';
        plot(ax1,FreqData{subj,cond}.time,SFreq.avg,'Color',colors(:,cond));
        hold(ax1,'on');
        
        ax2 = subplot(2,2,2);
        ax2.XAxisLocation = 'origin';
        ax2.YAxisLocation = 'origin';
        plot(ax2,TimeData{subj,cond}.time,STime.avg,'Color',colors(:,cond));
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
    title(ax3,{['Between Trial Correlation with p and rho']} );
    
    
    xlabel(ax3,'Time(secs)');
    ylabel(ax3, 'Rho' );
    title(ax3,{['Correlation by time']} )
    
    xlabel(ax4,'Total Power');
    ylabel(ax4, 'ERP Power (uV)' );
    
    title(ax4,{['Correlation Total Power-ERP by subject'] ['at times ' scattitle]})
    
    
end

end

