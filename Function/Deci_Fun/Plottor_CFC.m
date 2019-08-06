function Plottor_CFC(Deci)

%% Load

if isfield(Deci.Plot.CFC,'freqhigh')
    if ischar(Deci.Plot.CFC.freqhigh)
        switch Deci.Plot.CFC.freqhigh
            case 'theta'
                Deci.Plot.CFC.freqhigh = [4 8];
            case 'beta'
                Deci.Plot.CFC.freqhigh = [12.5 30];
            case 'alpha'
                Deci.Plot.CFC.freqhigh = [8 12.5];
            case 'lowgamma'
                Deci.Plot.CFC.freqhigh = [30 55];
            case 'highgamma'
                Deci.Plot.CFC.freqhigh = [55 80];
        end
    elseif isnumeric(Deci.Plot.CFC.freqhigh)
    else
        error(['cannot interrept freqhigh']);
    end
else
    error(['cannot interrept freqhigh']);
end

if isfield(Deci.Plot.CFC,'freqlow')
    if ischar(Deci.Plot.CFC.freqlow)
        switch Deci.Plot.CFC.freqlow
            case 'theta'
                Deci.Plot.CFC.freqlow = [4 8];
            case 'beta'
                Deci.Plot.CFC.freqlow = [12.5 30];
            case 'alpha'
                Deci.Plot.CFC.freqlow = [8 12.5];
            case 'lowgamma'
                Deci.Plot.CFC.freqlow = [30 55];
            case 'highgamma'
                Deci.Plot.CFC.freqlow = [55 80];
        end
    elseif isnumeric(Deci.Plot.CFC.freqlow)
    else
        error(['cannot interrept freqlow']);
    end
else
    error(['cannot interrept freqlow']);
end


for  subject_list = 1:length(Deci.SubjectList)
    tic;
    for Conditions = 1:length(Deci.Plot.IndexTitle)
        cfc = [];
        
        
        for m = 1:length(Deci.Plot.CFC.methods)
            load([Deci.Folder.Analysis filesep 'CFC' filesep Deci.Plot.CFC.methods{m} filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.IndexTitle{Conditions}],'cfc');
            
            highfoi = cfc.freqhigh >= round(Deci.Plot.CFC.freqhigh(1),4) & cfc.freqhigh <= round(Deci.Plot.CFC.freqhigh(2),4);
            
            cfc.freqhigh = cfc.freqhigh(highfoi);
            cfc.cfcdata  = cfc.cfcdata(:,:,:,highfoi,:);
            
            lowfoi = cfc.freqlow >= round(Deci.Plot.CFC.freqlow(1),4) & cfc.freqlow <= round(Deci.Plot.CFC.freqlow(2),4);
            
            cfc.freqlow = cfc.freqlow(lowfoi);
            cfc.cfcdata  = cfc.cfcdata(:,:,lowfoi,:,:);
            
            
            Subjects{subject_list,Conditions,m} = cfc;
        end
        
    end
    
end


if ~isempty(Deci.Plot.Math)
    
    for m = 1:length(Deci.Plot.CFC.methods)
        for cond = 1:length(Deci.Plot.Math)
            for subj = 1:size(Subjects,1)
                operation = str2fun(['@(x)' regexprep(Deci.Plot.Math{cond},'x(\d*)','x{$1}')]);
                MathData{subj} = Subjects{subj,1,m};
                MathData{subj}.crsspctrm = [feval(operation, cellfun(@(c) c.crsspctrm,Subjects(subj,1:length(Deci.Plot.IndexTitle)+cond-1,m),'UniformOutput',false))];
            end
            Subjects(:,length(Deci.Plot.IndexTitle)+cond,m) = MathData;
            TotalCount(:,length(Deci.Plot.IndexTitle)+cond) = num2cell(nan(size(Subjects,1),1));
        end
    end
    
end

if Deci.Plot.GrandAverage
    
    for conds = 1:size(Subjects,2)
        for m = 1:length(Deci.Plot.CFC.methods)
            facfg.parameter =  'crsspctrm';
            facfg.type = 'mean';
            CFCData{1,conds,m} = Subjects{1,conds,m};
            CFCData{1,conds,m}.crsspctrm =  mean([cell2mat(cellfun(@(c) c.crsspctrm,Subjects(:,conds,m),'UniformOutput',false))],1);
            
            %TotalCount{1,conds} = mean([TrialCount{:,conds}]);
            
            if Deci.Plot.CFC.errorbars
                facfg.type = 'std';
                CFCStd{1,conds,m} = Subjects{1,conds,m};
                CFCStd{1,conds,m}.crsspsctrm = std([cell2mat(cellfun(@(c) c.crsspctrm,Subjects(:,conds,m),'UniformOutput',false))],1);
            end
            
            
        end
    end
    
    Deci.SubjectList = {'Group 1'};
else
    CFCData = Subjects;
    TotalCount = TrialCount;
    
    
    
    
end
clear Subjects;




for cond = 1:length(Deci.Plot.Draw)
    for met = 1:size(CFCData,3)
        for subj = 1:size(CFCData,1)
            
            % Chan x Chan plot, Time has a dial
            if Deci.Plot.CFC.Square
                CFCsquare(subj,met) = figure;
            end
            
            % Topo plot with lines for connections, Time has a Dial
            if Deci.Plot.CFC.Topo
                CFCtopo(subj,met)  = figure;
            end
            
            % Each Bar is a ChanxChan value, Y axis represent indpedent
            % unit. Bars are in group by timebins.
            if Deci.Plot.CFC.Hist
                CFChist(subj,met)  = figure;
            end
            
            % Same as Bar but represent as line plot instead of bar plots
            if Deci.Plot.CFC.Wire
                CFCwire(subj,met)  = figure;
            end
            
            % Phase of low freq vs power of high freq
            if Deci.Plot.CFC.Phase
                for chanh = 1:length(CFCData{subj,cond,met}.labelhigh)
                    CFCphase(subj,met,chanh)  = figure;
                end
            end
            
            % Same as previous plot but broken down by frequency
            
            for subcond = 1:length(Deci.Plot.Draw{cond})
                
                if Deci.Plot.CFC.Spectrum
                    for chanh = 1:length(CFCData{subj,subcond,met}.labelhigh)
                        CFCspectrum(subj,met,subcond,chanh)  = figure;
                    end
                end
                
                if Deci.Plot.CFC.Topo
                    set(0, 'CurrentFigure', CFCtopo(subj,met) )
                    CFCtopo(subj,met).Visible = 'on';
                    
                    ctopo(subj,subcond,met)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                    
                    
                    Cross.labelcmb = CombVec(CFCData{subj,subcond,1}.labellow',CFCData{subj,subcond,1}.labelhigh')';
                    Cross.freq = CFCData{subj,subcond,1}.freqlow;
                    Cross.cohspctrm = CFCData{subj,subcond,met}.crsspctrm;
                    Cross.dimord = CFCData{subj,subcond,1}.dimord;
                    
                    Crosscfg =[];
                    Crosscfg.foi = [min([cfc.freqhigh cfc.freqlow]) max([cfc.freqhigh cfc.freqlow])];
                    Crosscfg.layout = Deci.Layout.Noeye;
                    title([Deci.SubjectList{subj} ' ' Deci.Plot.CFC.methods{met} ' Cond '  num2str(cond)],'Interpreter', 'none');
                    ft_topoplotCC(Crosscfg,Cross);
                    
                    ctopo(subj,subcond,met).UserData = {Cross,Crosscfg,CFCData{subj,subcond,met}.crsspctrm};
                end
                
                if Deci.Plot.CFC.Hist
                    set(0, 'CurrentFigure', CFChist(subj,met) )
                    CFChist(subj,met).Visible = 'on';
                    
                    chist(subj,subcond,met)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                    chist(subj,subcond,met).UserData = CFCData{subj,subcond,met}.crsspctrm;
                    CleanBars(reshape(CFCData{subj,subcond,met}.crsspctrm,[size(CFCData{subj,subcond,met}.crsspctrm,1)*size(CFCData{subj,subcond,met}.crsspctrm,2) size(CFCData{subj,subcond,met}.crsspctrm,3)]));
                    
                    title([Deci.SubjectList{subj} ' ' Deci.Plot.CFC.methods{met} ' Cond '  num2str(cond)],'Interpreter', 'none');
                    
                end
                
                
                if Deci.Plot.CFC.Square
                    
                    set(0, 'CurrentFigure', CFCsquare(subj,met) )
                    CFCsquare(subj,met).Visible = 'on';
                    csquare(subj,subcond,met)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                    
                    xdat = CFCData{subj,subcond,1}.labellow;
                    ydat = CFCData{subj,subcond,1}.labelhigh;
                    PLVim = imagesc(1:length(xdat),1:length(ydat),CFCData{subj,subcond,met}.crsspctrm(:,:,1));
                    xticks(1:length(xdat));
                    xticklabels(xdat);
                    yticks(1:length(ydat))
                    yticklabels(ydat);
                    title([Deci.SubjectList{subj} ' ' Deci.Plot.CFC.methods{met} ' Cond '  num2str(cond)],'Interpreter', 'none');
                    
                    csquare(subj,subcond,met).UserData = CFCData{subj,subcond,met}.crsspctrm;
                    colorbar;
                end
                
                if Deci.Plot.CFC.Wire
                    set(0, 'CurrentFigure', CFCwire(subj,met) )
                    CFCwire(subj,met).Visible = 'on';
                    
                    
                    for chanh = 1:length(CFCData{subj,subcond,met}.labelhigh)
                        for chanl = 1:length(CFCData{subj,subcond,met}.labellow)
                            
                            cwire(subj,chanl,met)    =  subplot(length(CFCData{subj,subcond,met}.labellow),1,chanl);
                            
                            
                            if ~Deci.Plot.CFC.errorbars
                                plot(CFCData{subj,subcond,met}.timelow,squeeze(CFCData{subj,subcond,met}.crsspctrm(chanl,chanh,:)))
                                
                            else
                                top = squeeze(mean(CFCData{subj,Deci.Plot.Draw{cond}(subcond),met}.crsspctrm(chanl,chanh,:),1)) + squeeze(mean(CFCStd{subj,Deci.Plot.Draw{cond}(subcond),met}.crsspctrm(chanl,chanh,:),1));
                                bot = squeeze(mean(CFCData{subj,Deci.Plot.Draw{cond}(subcond),met}.crsspctrm(chanl,chanh,:),1)) - squeeze(mean(CFCStd{subj,Deci.Plot.Draw{cond}(subcond),met}.crsspctrm(chanl,chanh,:),1));
                                
                                pgon = polyshape([CFCData{subj,Deci.Plot.Draw{cond}(subcond)}.timelow fliplr(CFCData{subj,Deci.Plot.Draw{cond}(subcond)}.timelow)],[top' fliplr(bot')],'Simplify', false);
                                b = plot(pgon,'HandleVisibility','off');
                                hold on
                                b.EdgeAlpha = 0;
                                b.FaceAlpha = .15;
                                h =  plot(CFCData{subj,subcond,met}.timelow,squeeze(CFCData{subj,subcond,met}.crsspctrm(chanl,chanh,:)));
                                h.Color = b.FaceColor;
                                h.LineWidth = 1;
                                
                            end
                            hold on
                        end
                    end
                end
                
                if Deci.Plot.CFC.Phase
                    
                    
                    
                    for chanh = 1:length(CFCData{subj,subcond,met}.labelhigh)
                        set(0, 'CurrentFigure', CFCphase(subj,met,chanh) )
                        CFCphase(subj,met,chanh).Visible = 'on';
                        for chanl = 1:length(CFCData{subj,subcond,met}.labellow)
                            
                            cphase(subj,chanl,met)    =  subplot(length(CFCData{subj,subcond,met}.labellow),1,chanl);
                            
                            if ~Deci.Plot.CFC.errorbars
                                plot(linspace(-pi,pi,size(CFCData{subj,subcond,met}.cfcdata,5)),squeeze(mean(CFCData{subj,subcond,met}.cfcdata(chanl,chanh,:,:,:),[1 2 3 4])))
                                
                            else
                                top = squeeze(squeeze(mean(normalize(CFCData{subj,subcond,met}.cfcdata(chanl,chanh,:,:,:),5),[1 2 3 4])) + squeeze(std(normalize(CFCData{subj,subcond,met}.cfcdata(chanl,chanh,:,:,:),5),0,[1 2 3 4])));
                                bot = squeeze(squeeze(mean(normalize(CFCData{subj,subcond,met}.cfcdata(chanl,chanh,:,:,:),5),[1 2 3 4])) - squeeze(std(normalize(CFCData{subj,subcond,met}.cfcdata(chanl,chanh,:,:,:),5),0,[1 2 3 4])));
                                
                                pgon = polyshape([linspace(-pi,pi,size(CFCData{subj,subcond,met}.cfcdata,5)) fliplr(linspace(-pi,pi,size(CFCData{subj,subcond,met}.cfcdata,5)))],[top' fliplr(bot')],'Simplify', false);
                                b = plot(pgon,'HandleVisibility','off');
                                hold on
                                b.EdgeAlpha = 0;
                                b.FaceAlpha = .15;
                                h =  plot(linspace(-pi,pi,size(CFCData{subj,subcond,met}.cfcdata,5)),squeeze(mean(normalize(CFCData{subj,subcond,met}.cfcdata(chanl,chanh,:,:,:),5),[1 2 3 4])));
                                h.Color = b.FaceColor;
                                h.LineWidth = 1;
                                
                            end
                            hold on
                            xlabel('Phase (radians)');
                            ylabel('Amplitude (modulation index)');
                            legend([Deci.Plot.Subtitle{cond}]);
                            title([CFCData{subj,subcond,met}.labellow(chanl), CFCData{subj,subcond,met}.labelhigh(chanh)])
                            
                        end
                    end
                end
                
                
                if Deci.Plot.CFC.Spectrum
                    cols = length(CFCData{subj,subcond,met}.labellow);
                    rows = size(CFCData{subj,subcond,met}.cfcdata,3);
                    
                    for chanh = 1:length(CFCData{subj,subcond,met}.labelhigh)
                        set(0, 'CurrentFigure', CFCspectrum(subj,met,subcond,chanh) )
                        CFCspectrum(subj,met,subcond,chanh).Visible = 'on';
                        iterator = 0;
                        for chanl = 1:length(CFCData{subj,subcond,met}.labellow)
                            for row = 1:size(CFCData{subj,subcond,met}.cfcdata,3)
                                iterator = iterator + 1;
                                cspectrum(subj,chanl,met)    =  subplot(rows,cols,iterator);
                                imagesc(linspace(-pi,pi,size(CFCData{subj,subcond,met}.cfcdata,5)),CFCData{subj,subcond,met}.freqhigh,normalize(squeeze(CFCData{subj,subcond,met}.cfcdata(chanl, chanh, row, :,:)),2))
                                xlabel('Phase (radians)');
                                ylabel('High frequency');
                                title([CFCData{subj,subcond,met}.freqlow(row), CFCData{subj,subcond,met}.labellow(chanl), CFCData{subj,subcond,met}.labelhigh(chanh) Deci.Plot.Subtitle{cond}{subcond}]);
                            end
                            hold on
                        end
                    end
                end
                
            end
            
            
            
            if Deci.Plot.CFC.Hist
                set(0, 'CurrentFigure', CFChist(subj,met));
                UpdateAxes(chist(subj,:,met),Deci.Plot.CFC.Roi,'Y',0)
                suptitle(Deci.Plot.Title{cond});
                
                for HistLim = 1:length(chist(subj,:,met))
                    
                    chist(subj,HistLim,met).YLim(1) = chist(subj,HistLim,met).YLim(1) *.98;
                    chist(subj,HistLim,met).YLim(2) = chist(subj,HistLim,met).YLim(2) *1.02;
                    
                    Vect = CombVec(CFCData{subj,subcond,met}.labellow',CFCData{subj,subcond,met}.labelhigh')';
                    
                    xticklabels(chist(subj,HistLim,met), arrayfun(@(c) [Vect{c,1} Vect{c,2}],1:size(Vect,1),'un',0));
                    xtickangle(chist(subj,HistLim,met),-20);
                    
                    legend(chist(subj,HistLim,met),[repmat('FreqLow Time ',[size(CFCData{subj,subcond,met}.timelow',1) 1]) num2str([CFCData{subj,subcond,met}.timelow']) repmat(' - FreqHigh Time ',[size(CFCData{subj,subcond,met}.timelow',1) 1]) num2str([CFCData{subj,subcond,met}.timelow'])])
                end
            end
            
            if Deci.Plot.CFC.Wire
                for WireLim = 1:length(cwire(subj,:,met))
                    Vect = CombVec(CFCData{subj,subcond,met}.labellow',CFCData{subj,subcond,met}.labelhigh')';
                    legend(cwire(subj,WireLim,met),arrayfun(@(c) [Vect{c,1} Vect{c,2}],1:size(Vect,1),'un',0));
                    cwire(subj,WireLim,met).YLim = [min([cwire(subj,:,met).YLim]) max([cwire(subj,:,met).YLim])];
                end
            end
            
            %             if Deci.Plot.CFC.Phase
            %                 for PhaseLim = 1:length(cphase(subj,:,met))
            %                     Vect = CombVec(CFCData{subj,subcond,met}.labellow',CFCData{subj,subcond,met}.labelhigh')';
            %                     legend(cphase(subj,PhaseLim,met),arrayfun(@(c) [Vect{c,1} Vect{c,2}],1:size(Vect,1),'un',0));
            %                     cphase(subj,PhaseLim,met).YLim = [min([cphase(subj,:,met).YLim]) max([cphase(subj,:,met).YLim])];
            %                 end
            %             end
            
            if Deci.Plot.CFC.Topo
                
                set(0, 'CurrentFigure', CFCtopo(subj,met) )
                suptitle(Deci.Plot.Title{cond});
                uicontrol('style','text','position',[225 75 100 25],'String','Time of Interest');
                
                slide = uicontrol('style','slider','position',[75 10 400 20],...
                    'min',1,'max',size(CFCData{1,subj,met}.crsspctrm,3),'callback',{@ChangeDimTopo,ctopo(subj,:,met),Deci.Plot.CFC.Roi}, ...
                    'value',1,'SliderStep',[1/size(CFCData{1,subj,met}.crsspctrm,3) 1/size(CFCData{1,subj,met}.crsspctrm,3)]);
                
                ChangeDimTopo(slide,[],ctopo(subj,:,met),Deci.Plot.CFC.Roi);
                
                
                for tick = 1:size(CFCData{1,subj,met}.crsspctrm,3)
                    uicontrol('style','text','position',[75+[[[slide.Position(3)]/5]*[tick-1]]+20 55 40 25],'String',num2str(round(CFCData{1,subj,met}.timelow(tick),2)));
                    uicontrol('style','text','position',[75+[[[slide.Position(3)]/5]*[tick-1]]+20 30 40 25],'String',num2str(round(CFCData{1,subj,met}.timehigh(tick),2)));
                end
                uicontrol('style','text','position',[45 55 60 25],'String','FreqLow');
                uicontrol('style','text','position',[45 30 60 25],'String','FreqHigh');
            end
            
            
            if Deci.Plot.CFC.Square
                
                set(0, 'CurrentFigure', CFCsquare(subj,met) )
                suptitle(Deci.Plot.Title{cond});
                uicontrol('style','text','position',[225 75 100 25],'String','Time of Interest');
                
                slide = uicontrol('style','slider','position',[75 10 400 20],...
                    'min',1,'max',size(CFCData{1,subj,met}.crsspctrm,3),'callback',{@ChangeDim,csquare(subj,:,met),Deci.Plot.CFC.Roi,'C'}, ...
                    'value',1,'SliderStep',[1/size(CFCData{1,subj,met}.crsspctrm,3) 1/size(CFCData{1,subj,met}.crsspctrm,3)]);
                
                UpdateAxes(csquare(subj,:,met),Deci.Plot.CFC.Roi,'C',1);
                
                for tick = 1:size(CFCData{1,subj}.crsspctrm,3)
                    uicontrol('style','text','position',[75+[[[slide.Position(3)]/5]*[tick-1]]+20 55 40 25],'String',num2str(round(CFCData{1,subj,met}.timelow(tick),2)));
                    uicontrol('style','text','position',[75+[[[slide.Position(3)]/5]*[tick-1]]+20 30 40 25],'String',num2str(round(CFCData{1,subj,met}.timehigh(tick),2)));
                end
                uicontrol('style','text','position',[45 55 60 25],'String','FreqLow');
                uicontrol('style','text','position',[45 30 60 25],'String','FreqHigh');
            end
            
        end
        
    end
end

    function ChangeDim(popup,event,Axes,Roi,Lim)
        popup.Value = round(popup.Value);
        
        for i = 1:length(Axes)
            Axes(i).Children.CData = Axes(i).UserData(:,:,popup.Value);
        end
        
        UpdateAxes(Axes,Roi,Lim,1)
        
    end

    function UpdateAxes(Axes,Roi,Lim,Userdata)
        
        if Userdata == 1
            Dats = cell2mat(arrayfun(@(c) c.UserData,Axes,'un',0));
        else
            Dats = cell2mat(arrayfun(@(c) [c.Children.([Lim 'Data'])],Axes,'un',0));
        end
        
        for Axe = 1:length(Axes(:))
            if isequal(Roi,'maxmin')
                Axes(Axe).([Lim 'Lim']) = [min(Dats(:)) max(Dats(:))];
            elseif isequal(Roi,[0 1])
                Axes(Axe).([Lim 'Lim']) = [0 1];
            elseif length(Roi) == 2 && isnumeric(Roi)
                Axes(Axe).([Lim 'Lim']) = Roi;
            end
        end
    end

    function ChangeDimTopo(popup,event,Axes,Roi)
        popup.Value = round(popup.Value);
        
        if isequal(Roi,'maxmin')
            CLim = cell2mat(arrayfun(@(c) c.UserData{3},Axes,'un',0));
        elseif isequal(Roi,[0 1])
            CLim = [0 1];
        elseif length(Roi) == 2 && isnumeric(Roi)
            CLim = Roi;
        end
        
        for Axe =  1:length(Axes)
            set(Axes(Axe).Parent, 'currentaxes', Axes(Axe))
            
            NewCross = Axes(Axe).UserData{1};
            NewCross.cohspctrm =  Axes(Axe).UserData{3}(:,:,popup.Value);
            NewCrossCfg = Axes(Axe).UserData{2};
            NewCrossCfg.CLim = [min(CLim(:)') max(CLim(:)')];
            ft_topoplotCC(NewCrossCfg ,NewCross);
            
        end
        
    end

end