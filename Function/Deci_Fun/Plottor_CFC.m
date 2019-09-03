function Plottor_CFC(Deci)

%% Load

for  subject_list = 1:length(Deci.SubjectList)
    tic;
    for Conditions = 1:length(Deci.Plot.IndexTitle)
        cfc = [];
        
        
        for m = 1:length(Deci.Plot.CFC.methods)
            load([Deci.Folder.Analysis filesep 'CFC' filesep Deci.Plot.CFC.methods{m} filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.IndexTitle{Conditions}],'cfc');
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
            
            CFCSem{1,conds,m} = Subjects{1,conds,m};
            CFCSem{1,conds,m}.crsspctrm = std([cell2mat(cellfun(@(c) c.crsspctrm,Subjects(:,conds,m),'UniformOutput',false))],[],1)/sqrt(size(Subjects,1));
        end
    end
    
    Deci.SubjectList = {'Group 1'};
else
    CFCData = Subjects;
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
            %             if Deci.Plot.CFC.Hist
            %                 CFChist(subj,met)  = figure;
            %             end
            
            % Same as Bar but represent as line plot instead of bar plots
            if Deci.Plot.CFC.Wire
                CFCwire(subj,met)  = figure;
            end
            
            % Phase of low freq vs power of high freq
            %             if Deci.Plot.CFC.Phase
            %                 for chanh = 1:length(CFCData{subj,cond,met}.labelhigh)
            %                     CFCphase(subj,met,chanh)  = figure;
            %                 end
            %             end
            
            % Same as previous plot but broken down by frequency
            
            for subcond = 1:length(Deci.Plot.Draw{cond})
                
                %                 if Deci.Plot.CFC.Spectrum
                %                     for chanh = 1:length(CFCData{subj,subcond,met}.labelhigh)
                %                         CFCspectrum(subj,met,subcond,chanh)  = figure;
                %                     end
                %                 end
                
                if Deci.Plot.CFC.Topo
                    set(0, 'CurrentFigure', CFCtopo(subj,met) )
                    CFCtopo(subj,met).Visible = 'on';
                    
                    ctopo(subj,subcond,met)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                    
                    Crosscfg =[];
                    Crosscfg.layout = Deci.Layout.Noeye;
                    Crosscfg.toi = 1;
                    Crosscfg.foi1 = 1;
                    Crosscfg.foi2 = 1;
                    
                    title([Deci.SubjectList{subj} ' ' Deci.Plot.CFC.methods{met} ' Cond '  num2str(cond)],'Interpreter', 'none');
                    ft_topoplotCC_rob(Crosscfg,CFCData{subj,subcond,met});
                    
%                     Cross = [];
%                     Cross.cohspctrm = CFCData{subj,subcond,met}.crsspctrm(:,:,1,1,1);
%                     Cross.freq = 4;
%                     Cross.time = CFCData{subj,subcond,met}.timelow(1);
%                     Cross.dimord = 'chan_chan_freq';
%                     Cross.labellow = CFCData{subj,subcond,met}.labellow;
%                     Cross.labelhigh = CFCData{subj,subcond,met}.labelhigh;
%                     
%                     ft_connectivityplot(Crosscfg,Cross)
                    
                    ctopo(subj,subcond,met).UserData = {CFCData{subj,subcond,met},Crosscfg};
                end
                %
                %                 if Deci.Plot.CFC.Hist
                %                     set(0, 'CurrentFigure', CFChist(subj,met) )
                %                     CFChist(subj,met).Visible = 'on';
                %
                %                     chist(subj,subcond,met)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                %                     chist(subj,subcond,met).UserData = CFCData{subj,subcond,met}.crsspctrm;
                %                     CleanBars(reshape(CFCData{subj,subcond,met}.crsspctrm,[size(CFCData{subj,subcond,met}.crsspctrm,1)*size(CFCData{subj,subcond,met}.crsspctrm,2) size(CFCData{subj,subcond,met}.crsspctrm,3)]));
                %
                %                     title([Deci.SubjectList{subj} ' ' Deci.Plot.CFC.methods{met} ' Cond '  num2str(cond)],'Interpreter', 'none');
                %
                %                 end
                %
                
                if Deci.Plot.CFC.Square
                    
                    set(0, 'CurrentFigure', CFCsquare(subj,met) )
                    CFCsquare(subj,met).Visible = 'on';
                    csquare(subj,subcond,met)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                    
                    xdat = CFCData{subj,subcond,1}.labellow;
                    ydat = CFCData{subj,subcond,1}.labelhigh;
                    imagesc(1:length(xdat),1:length(ydat),CFCData{subj,subcond,met}.crsspctrm(:,:,1,1,1)');
                    xticks(1:length(xdat));
                    xticklabels(xdat);
                    yticks(1:length(ydat))
                    yticklabels(ydat);
                    title([Deci.SubjectList{subj} ' ' Deci.Plot.CFC.methods{met} ' Cond '  num2str(cond)],'Interpreter', 'none');
                    Crosscfg.toi = 1;
                    Crosscfg.foi1 = 1;
                    Crosscfg.foi2 = 1;
                    
                    
                    csquare(subj,subcond,met).UserData = {CFCData{subj,subcond,met},Crosscfg};
                    colorbar;
                end
                
                if Deci.Plot.CFC.Wire
                    set(0, 'CurrentFigure', CFCwire(subj,met) )
                    CFCwire(subj,met).Visible = 'on';
                    Crosscfg.toi = 1;
                    Crosscfg.foi1 = 1;
                    Crosscfg.foi2 = 1;
                    
                    cwire(subj,subcond,met)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                    
%                     top = squeeze(mean(CFCData{subj,Deci.Plot.Draw{cond}(subcond),met}.crsspctrm(:,:,:),1)) + squeeze(mean(CFCSem{subj,Deci.Plot.Draw{cond}(subcond),met}.crsspctrm(chanl,chanh,:),1));
%                     bot = squeeze(mean(CFCData{subj,Deci.Plot.Draw{cond}(subcond),met}.crsspctrm(chanl,chanh,:),1)) - squeeze(mean(CFCSem{subj,Deci.Plot.Draw{cond}(subcond),met}.crsspctrm(chanl,chanh,:),1));
%                     
%                     pgon = polyshape([CFCData{subj,Deci.Plot.Draw{cond}(subcond)}.timelow fliplr(CFCData{subj,Deci.Plot.Draw{cond}(subcond)}.timelow)],[top' fliplr(bot')],'Simplify', false);
%                     b = plot(pgon,'HandleVisibility','off');
%                     hold on
%                     b.EdgeAlpha = 0;
%                     b.FaceAlpha = .15;
                    
                    h =  plot(CFCData{subj,subcond,met}.timelow,squeeze(CFCData{subj,subcond,met}.crsspctrm(:,:,1,1,:)));
                    title([Deci.SubjectList{subj} ' ' Deci.Plot.CFC.methods{met} ' Cond '  num2str(cond)],'Interpreter', 'none');
%                     h.Color = b.FaceColor;
                    arrayfun(@(c) set(c,'lineWidth',1),h)
                    cwire(subj,subcond,met).UserData = {CFCData{subj,subcond,met},Crosscfg};
                    
                    hold on
                end
                
                %                 if Deci.Plot.CFC.Phase
                %
                %
                %
                %                     for chanh = 1:length(CFCData{subj,subcond,met}.labelhigh)
                %                         set(0, 'CurrentFigure', CFCphase(subj,met,chanh) )
                %                         CFCphase(subj,met,chanh).Visible = 'on';
                %                         for chanl = 1:length(CFCData{subj,subcond,met}.labellow)
                %
                %                             cphase(subj,chanl,met)    =  subplot(length(CFCData{subj,subcond,met}.labellow),1,chanl);
                %
                %                             if ~Deci.Plot.CFC.errorbars
                %                                 plot(linspace(-pi,pi,size(CFCData{subj,subcond,met}.cfcdata,5)),squeeze(mean(CFCData{subj,subcond,met}.cfcdata(chanl,chanh,:,:,:),[1 2 3 4])))
                %
                %                             else
                %                                 top = squeeze(squeeze(mean(normalize(CFCData{subj,subcond,met}.cfcdata(chanl,chanh,:,:,:),5),[1 2 3 4])) + squeeze(std(normalize(CFCData{subj,subcond,met}.cfcdata(chanl,chanh,:,:,:),5),0,[1 2 3 4])));
                %                                 bot = squeeze(squeeze(mean(normalize(CFCData{subj,subcond,met}.cfcdata(chanl,chanh,:,:,:),5),[1 2 3 4])) - squeeze(std(normalize(CFCData{subj,subcond,met}.cfcdata(chanl,chanh,:,:,:),5),0,[1 2 3 4])));
                %
                %                                 pgon = polyshape([linspace(-pi,pi,size(CFCData{subj,subcond,met}.cfcdata,5)) fliplr(linspace(-pi,pi,size(CFCData{subj,subcond,met}.cfcdata,5)))],[top' fliplr(bot')],'Simplify', false);
                %                                 b = plot(pgon,'HandleVisibility','off');
                %                                 hold on
                %                                 b.EdgeAlpha = 0;
                %                                 b.FaceAlpha = .15;
                %                                 h =  plot(linspace(-pi,pi,size(CFCData{subj,subcond,met}.cfcdata,5)),squeeze(mean(normalize(CFCData{subj,subcond,met}.cfcdata(chanl,chanh,:,:,:),5),[1 2 3 4])));
                %                                 h.Color = b.FaceColor;
                %                                 h.LineWidth = 1;
                %
                %                             end
                %                             hold on
                %                             xlabel('Phase (radians)');
                %                             ylabel('Amplitude (modulation index)');
                %                             legend([Deci.Plot.Subtitle{cond}]);
                %                             title([CFCData{subj,subcond,met}.labellow(chanl), CFCData{subj,subcond,met}.labelhigh(chanh)])
                %
                %                         end
                %                     end
                %                 end
                
                
                %                 if Deci.Plot.CFC.Spectrum
                %                     cols = length(CFCData{subj,subcond,met}.labellow);
                %                     rows = size(CFCData{subj,subcond,met}.cfcdata,3);
                %
                %                     for chanh = 1:length(CFCData{subj,subcond,met}.labelhigh)
                %                         set(0, 'CurrentFigure', CFCspectrum(subj,met,subcond,chanh) )
                %                         CFCspectrum(subj,met,subcond,chanh).Visible = 'on';
                %                         iterator = 0;
                %                         for chanl = 1:length(CFCData{subj,subcond,met}.labellow)
                %                             for row = 1:size(CFCData{subj,subcond,met}.cfcdata,3)
                %                                 iterator = iterator + 1;
                %                                 cspectrum(subj,chanl,met)    =  subplot(rows,cols,iterator);
                %                                 imagesc(linspace(-pi,pi,size(CFCData{subj,subcond,met}.cfcdata,5)),CFCData{subj,subcond,met}.freqhigh,normalize(squeeze(CFCData{subj,subcond,met}.cfcdata(chanl, chanh, row, :,:)),2))
                %                                 xlabel('Phase (radians)');
                %                                 ylabel('High frequency');
                %                                 title([CFCData{subj,subcond,met}.freqlow(row), CFCData{subj,subcond,met}.labellow(chanl), CFCData{subj,subcond,met}.labelhigh(chanh) Deci.Plot.Subtitle{cond}{subcond}]);
                %                             end
                %                             hold on
                %                         end
                %                     end
                %                 end
                
            end
            
            
            
            %             if Deci.Plot.CFC.Hist
            %                 set(0, 'CurrentFigure', CFChist(subj,met));
            %                 UpdateAxes(chist(subj,:,met),Deci.Plot.CFC.Roi,'Y',0)
            %                 suptitle(Deci.Plot.Title{cond});
            %
            %                 for HistLim = 1:length(chist(subj,:,met))
            %
            %                     chist(subj,HistLim,met).YLim(1) = chist(subj,HistLim,met).YLim(1) *.98;
            %                     chist(subj,HistLim,met).YLim(2) = chist(subj,HistLim,met).YLim(2) *1.02;
            %
            %                     Vect = CombVec(CFCData{subj,subcond,met}.labellow',CFCData{subj,subcond,met}.labelhigh')';
            %
            %                     xticklabels(chist(subj,HistLim,met), arrayfun(@(c) [Vect{c,1} Vect{c,2}],1:size(Vect,1),'un',0));
            %                     xtickangle(chist(subj,HistLim,met),-20);
            %
            %                     legend(chist(subj,HistLim,met),[repmat('FreqLow Time ',[size(CFCData{subj,subcond,met}.timelow',1) 1]) num2str([CFCData{subj,subcond,met}.timelow']) repmat(' - FreqHigh Time ',[size(CFCData{subj,subcond,met}.timelow',1) 1]) num2str([CFCData{subj,subcond,met}.timelow'])])
            %                 end
            %             end
            

            
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
                
                
                text = uicontrol('style','text','position',[225 75 100 25],'String',['Time of Interest: '  num2str(CFCData{subj,subcond,met}.timelow(1)) ' x ' num2str(CFCData{subj,subcond,met}.timehigh(1)) ]);
                uicontrol('style','pushbutton','position',[175 75 45 25],'String','<','callback',{@ChangeDim,ctopo(subj,:,met),-1,'toi',text,'Time of Interest: ','Topo'})
                uicontrol('style','pushbutton','position',[325 75 45 25],'String','>','callback',{@ChangeDim,ctopo(subj,:,met),+1,'toi',text,'Time of Interest: ','Topo'})
                
                
                text = uicontrol('style','text','position',[225+100+25+50+50 75 100 25],'String',['FreqLow of Interest: '  CFCData{subj,subcond,met}.freqlow{1}]);
                uicontrol('style','pushbutton','position',[175+100+25+50+50 75 45 25],'String','<','callback',{@ChangeDim,ctopo(subj,:,met),-1,'foi1',text,'FreqLow of Interest: ','Topo'})
                uicontrol('style','pushbutton','position',[325+100+25+50+50 75 45 25],'String','>','callback',{@ChangeDim,ctopo(subj,:,met),+1,'foi1',text,'FreqLow of Interest: ','Topo'})
                
                 ChangeDim([],[],ctopo(subj,:,met),0,'foi1',text,'FreqLow of Interest: ','Topo');
                
                text = uicontrol('style','text','position',[225+100+25+50+100+25+50+50+50 75 100 25],'String',['FreqHigh Interest: '   CFCData{subj,subcond,met}.freqhigh{1} ]);
                uicontrol('style','pushbutton','position',[175+100+25+50+100+25+50+50+50 75 45 25],'String','<','callback',{@ChangeDim,ctopo(subj,:,met),-1,'foi2',text,'FreqHigh Interest: ','Topo'})
                uicontrol('style','pushbutton','position',[325+100+25+50+100+25+50+50+50 75 45 25],'String','>','callback',{@ChangeDim,ctopo(subj,:,met),+1,'foi2',text,'FreqHigh Interest: ','Topo'})
                
            end
            
            
            if Deci.Plot.CFC.Square
                
                set(0, 'CurrentFigure', CFCsquare(subj,met) )
                suptitle(Deci.Plot.Title{cond});
                text = uicontrol('style','text','position',[225 75 100 25],'String',['Time of Interest: '  num2str(CFCData{subj,subcond,met}.timelow(1)) ' x ' num2str(CFCData{subj,subcond,met}.timehigh(1)) ]);
                
                uicontrol('style','pushbutton','position',[175 75 45 25],'String','<','callback',{@ChangeDim,csquare(subj,:,met),-1,'toi',text,'Time of Interest: ','Square'})
                uicontrol('style','pushbutton','position',[325 75 45 25],'String','>','callback',{@ChangeDim,csquare(subj,:,met),+1,'toi',text,'Time of Interest: ','Square'})
               
                
                text = uicontrol('style','text','position',[225+100+25+50+50 75 100 25],'String',['FreqLow of Interest: '  CFCData{subj,subcond,met}.freqlow{1}]);
                uicontrol('style','pushbutton','position',[175+100+25+50+50 75 45 25],'String','<','callback',{@ChangeDim,csquare(subj,:,met),-1,'foi1',text,'FreqLow of Interest: ','Square'})
                uicontrol('style','pushbutton','position',[325+100+25+50+50 75 45 25],'String','>','callback',{@ChangeDim,csquare(subj,:,met),+1,'foi1',text,'FreqLow of Interest: ','Square'})
                
                ChangeDim([],[],csquare(subj,:,met),0,'foi1',text,'FreqLow of Interest: ','Square');
                
                text = uicontrol('style','text','position',[225+100+25+50+100+25+50+50+50 75 100 25],'String',['FreqHigh Interest: '   CFCData{subj,subcond,met}.freqhigh{1} ]);
                uicontrol('style','pushbutton','position',[175+100+25+50+100+25+50+50+50 75 45 25],'String','<','callback',{@ChangeDim,csquare(subj,:,met),-1,'foi2',text,'FreqHigh Interest: ','Square'})
                uicontrol('style','pushbutton','position',[325+100+25+50+100+25+50+50+50 75 45 25],'String','>','callback',{@ChangeDim,csquare(subj,:,met),+1,'foi2',text,'FreqHigh Interest: ','Square'})
            end
            
            if Deci.Plot.CFC.Wire
                
                set(0, 'CurrentFigure', CFCwire(subj,met) )
                suptitle(Deci.Plot.Title{cond});
              
                text = uicontrol('style','text','position',[225+100+25+50+50 75 100 25],'String',['FreqLow of Interest: '  CFCData{subj,subcond,met}.freqlow{1}]);
                uicontrol('style','pushbutton','position',[175+100+25+50+50 75 45 25],'String','<','callback',{@ChangeDim,cwire(subj,:,met),-1,'foi1',text,'FreqLow of Interest: ','Wire'})
                uicontrol('style','pushbutton','position',[325+100+25+50+50 75 45 25],'String','>','callback',{@ChangeDim,cwire(subj,:,met),+1,'foi1',text,'FreqLow of Interest: ','Wire'})
                
                
                text = uicontrol('style','text','position',[225+100+25+50+100+25+50+50+50 75 100 25],'String',['FreqHigh Interest: '   CFCData{subj,subcond,met}.freqhigh{1} ]);
                uicontrol('style','pushbutton','position',[175+100+25+50+100+25+50+50+50 75 45 25],'String','<','callback',{@ChangeDim,cwire(subj,:,met),-1,'foi2',text,'FreqHigh Interest: ','Wire'})
                uicontrol('style','pushbutton','position',[325+100+25+50+100+25+50+50+50 75 45 25],'String','>','callback',{@ChangeDim,cwire(subj,:,met),+1,'foi2',text,'FreqHigh Interest: ','Wire'})
                
                
                
                for WireLim = 1:length(cwire(subj,:,met))
                    Vect = CombVec(CFCData{subj,subcond,met}.labellow',CFCData{subj,subcond,met}.labelhigh')';
                    legend(cwire(subj,WireLim,met),arrayfun(@(c) [Vect{c,1} Vect{c,2}],1:size(Vect,1),'un',0));
                    cwire(subj,WireLim,met).YLim = [min([cwire(subj,:,met).YLim]) max([cwire(subj,:,met).YLim])];
                end
            end
            
            
        end
        
    end
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

    function ChangeDim(popup,event,Axes,Change,Dim,TextUI,Text,Type)
        
        
        switch Dim
            
            case 'toi'
                for Axe =  1:length(Axes)
                    set(Axes(Axe).Parent, 'currentaxes', Axes(Axe))
                    
                    NewCross = Axes(Axe).UserData{1};
                    
                    NewCrossCfg = Axes(Axe).UserData{2};
                    NewCrossCfg.toi = NewCrossCfg.toi + Change;
                    
                    if NewCrossCfg.toi < 1 || NewCrossCfg.toi > length(Axes(Axe).UserData{1}.timelow)
                        break
                    end
                    
                    switch Type
                        case 'Topo'
                            ft_topoplotCC_rob(NewCrossCfg ,NewCross);
                            %A = findobj(Axes(Axe).Children,'Type','Patch');
                            Axes(Axe).CLim = [minmax([Axes(:).CLim])];
                            
                        case 'Square'
                            Axes(Axe).Children.CData = Axes(Axe).UserData{1}.crsspctrm(:,:,NewCrossCfg.foi1,NewCrossCfg.foi2,NewCrossCfg.toi);
                            Axes(Axe).CLim = [minmax([Axes(:).CLim])];
                    end
                    
                    Axes(Axe).UserData{2}.toi = Axes(Axe).UserData{2}.toi + Change;
                end
                
                if NewCrossCfg.toi >= 1 && NewCrossCfg.toi <= length(Axes(Axe).UserData{1}.timelow)
                    TextUI.String = [Text num2str(Axes(Axe).UserData{1}.timelow(NewCrossCfg.toi) ) ' x ' num2str(Axes(Axe).UserData{1}.timehigh(NewCrossCfg.toi))];
                end
            case 'foi1'
                for Axe =  1:length(Axes)
                    set(Axes(Axe).Parent, 'currentaxes', Axes(Axe))
                    
                    NewCross = Axes(Axe).UserData{1};
                    
                    NewCrossCfg = Axes(Axe).UserData{2};
                    NewCrossCfg.foi1 = NewCrossCfg.foi1 + Change;
                    
                    if NewCrossCfg.foi1 < 1 || NewCrossCfg.foi1 > length(Axes(Axe).UserData{1}.freqlow)
                        break
                    end
                    
                    switch Type
                        case 'Topo'
                            ft_topoplotCC_rob(NewCrossCfg ,NewCross);
                             Axes(Axe).CLim = [minmax([Axes(:).CLim])];
                        case 'Square'
                            Axes(Axe).Children.CData = Axes(Axe).UserData{1}.crsspctrm(:,:,NewCrossCfg.foi1,NewCrossCfg.foi2,NewCrossCfg.toi);
                            Axes(Axe).CLim = [minmax([Axes(:).CLim])];
                        case 'Wire'
                            
                            for w = 1:length(Axes(Axe).Children)
                                Axes(Axe).Children(w).YData = Axes(Axe).UserData{1}.crsspctrm(:,:,NewCrossCfg.foi1,NewCrossCfg.foi2,NewCrossCfg.toi)
                            end
                            
                    end
                    
                    Axes(Axe).UserData{2}.foi1 = Axes(Axe).UserData{2}.foi1 + Change;
                end
                
                if NewCrossCfg.foi1 >= 1 && NewCrossCfg.foi1 <= length(Axes(Axe).UserData{1}.freqlow)
                    TextUI.String = [Text Axes(Axe).UserData{1}.freqlow{NewCrossCfg.foi1}];
                end
            case 'foi2'
                for Axe =  1:length(Axes)
                    set(Axes(Axe).Parent, 'currentaxes', Axes(Axe))
                    
                    NewCross = Axes(Axe).UserData{1};
                    
                    NewCrossCfg = Axes(Axe).UserData{2};
                    NewCrossCfg.foi2 = NewCrossCfg.foi2 + Change;
                    
                    if NewCrossCfg.foi2 < 1 || NewCrossCfg.foi2 > length(Axes(Axe).UserData{1}.freqhigh)
                        break
                    end
                    
                    switch Type
                        
                        case 'Topo'
                            ft_topoplotCC_rob(NewCrossCfg ,NewCross);
                             Axes(Axe).CLim = [minmax([Axes(:).CLim])];
                        case 'Square'
                            Axes(Axe).Children.CData = Axes(Axe).UserData{1}.crsspctrm(:,:,NewCrossCfg.foi1,NewCrossCfg.foi2,NewCrossCfg.toi);
                            Axes(Axe).CLim = [minmax([Axes(:).CLim])];
                    end
                    
                    Axes(Axe).UserData{2}.foi2 = Axes(Axe).UserData{2}.foi2 + Change;
                end
                
                if NewCrossCfg.foi2 >= 1 && NewCrossCfg.foi2 <= length(Axes(Axe).UserData{1}.freqhigh)
                    TextUI.String = [Text Axes(Axe).UserData{1}.freqhigh{NewCrossCfg.foi2}];
                end
        end
    end

end