function Plottor_CFC(Deci)

%% Deci Checks

if isfield(Deci.Plot.CFC,'freqhigh')
    if ischar(Deci.Plot.CFC.freqhigh)
        switch Deci.Plot.CFC.freqhigh
            case 'theta'
                Deci.Plot.CFC.freqhigh = [4 8];
            case 'beta'
                Deci.Plot.CFC.freqhigh = [12.5 30];
            case 'alpha'
                Deci.Plot.CFC.freqhigh = [8 12.5];
            case 'gamma'
                Deci.Plot.CFC.freqhigh = [30 50];
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
            case 'gamma'
                Deci.Plot.CFC.freqlow = [30 50];
        end
    elseif isnumeric(Deci.Plot.CFC.freqlow)
    else
        error(['cannot interrept freqlow']);
    end
else
    error(['cannot interrept freqlow']);
end


if isequal(Deci.Plot.CFC.chanlow,'Reinhart-All')
    Deci.Plot.CFC.chanlow = [{'AF3'  } {'AF4'  } {'AF7'  } ...
        {'AF8'  } {'AFz'  } {'C1'   } {'C2'   } {'C3'   } {'C4'   } {'C5'   } ...
        {'C6'   } {'CP1'  } {'CP2'  } {'CP3'  } {'CP4'  } {'CP5'  } {'CP6'  } ...
        {'CPz'  } {'Cz'   } {'F1'   } {'F2'   } {'F3'   } {'F4'   } {'F5'   } ...
        {'F6'   } {'F7'   } {'F8'   } {'FC1'  } {'FC2'  } {'FC3'  } {'FC4'  } ...
        {'FC5'  } {'FC6'  } {'FCz'  } {'FT7'  } {'FT8'  } {'Fz'   } {'O1'   } ...
        {'O2'   } {'Oz'   } {'P1'   } {'P2'   } {'P3'   } {'P4'   } {'P5'   } ...
        {'P6'   } {'P7'   } {'P8'   } {'PO3'  } {'PO4'  } {'PO7'  } {'PO8'  } ...
        {'POz'  } {'Pz'   } {'T7'   } {'T8'   } {'TP10' } {'TP7'  } {'TP8'  } ...
        {'TP9'  } ] ;
end

if isequal(Deci.Plot.CFC.chanhigh,'Reinhart-All')
    Deci.Plot.CFC.chanhigh = [{'AF3'  } {'AF4'  } {'AF7'  } ...
        {'AF8'  } {'AFz'  } {'C1'   } {'C2'   } {'C3'   } {'C4'   } {'C5'   } ...
        {'C6'   } {'CP1'  } {'CP2'  } {'CP3'  } {'CP4'  } {'CP5'  } {'CP6'  } ...
        {'CPz'  } {'Cz'   } {'F1'   } {'F2'   } {'F3'   } {'F4'   } {'F5'   } ...
        {'F6'   } {'F7'   } {'F8'   } {'FC1'  } {'FC2'  } {'FC3'  } {'FC4'  } ...
        {'FC5'  } {'FC6'  } {'FCz'  } {'FT7'  } {'FT8'  } {'Fz'   } {'O1'   } ...
        {'O2'   } {'Oz'   } {'P1'   } {'P2'   } {'P3'   } {'P4'   } {'P5'   } ...
        {'P6'   } {'P7'   } {'P8'   } {'PO3'  } {'PO4'  } {'PO7'  } {'PO8'  } ...
        {'POz'  } {'Pz'   } {'T7'   } {'T8'   } {'TP10' } {'TP7'  } {'TP8'  } ...
        {'TP9'  } ] ;
end

%% Load

for  subject_list = 1:length(Deci.SubjectList)
    tic;
    for Channel = 1:length(Deci.Plot.CFC.chanhigh)
        
        display(['High Channel ' Deci.Plot.CFC.chanhigh{Channel} ' of ' num2str(length(Deci.Plot.CFC.chanhigh))])
        
        freq = [];
        load([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Lock filesep Deci.Plot.CFC.chanhigh{Channel}],'freq');
        
        foi = freq.freq >= round(Deci.Plot.CFC.freqhigh(1),4) & freq.freq <= round(Deci.Plot.CFC.freqhigh(2),4);
        toi = round(freq.time,4) >= Deci.Plot.CFC.latencyhigh(1) & round(freq.time,4) <= Deci.Plot.CFC.latencyhigh(2);
        
        HighChans{Channel} = freq;
        HighChans{Channel}.freq =  HighChans{Channel}.freq(foi);
        HighChans{Channel}.time =  HighChans{Channel}.time(toi);
        HighChans{Channel}.fourierspctrm  = HighChans{Channel}.fourierspctrm(:,:,foi,toi);
    end
    toc;
    
    acfg.parameter = 'fourierspctrm';
    acfg.appenddim = 'chan';
    HighChans = rmfield(ft_appendfreq(acfg,HighChans{:}),'cfg');
    
    tic;
    for Channel = 1:length(Deci.Plot.CFC.chanlow)
        
        display(['Low Channel ' Deci.Plot.CFC.chanlow{Channel} ' of ' num2str(length(Deci.Plot.CFC.chanlow))])
        
        freq = [];
        load([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Lock filesep Deci.Plot.CFC.chanlow{Channel}],'freq');
        
        foi = freq.freq >= round(Deci.Plot.CFC.freqlow(1),4) & freq.freq <= round(Deci.Plot.CFC.freqlow(2),4);
        toi = round(freq.time,4) >= Deci.Plot.CFC.latencylow(1) & round(freq.time,4) <= Deci.Plot.CFC.latencylow(2);
        
        LowChans{Channel} = freq;
        LowChans{Channel}.freq =  LowChans{Channel}.freq(foi);
        LowChans{Channel}.time =  LowChans{Channel}.time(toi);
        LowChans{Channel}.fourierspctrm  = LowChans{Channel}.fourierspctrm(:,:,foi,toi);
        
    end
    toc;
    
    acfg.parameter = 'fourierspctrm';
    acfg.appenddim = 'chan';
    LowChans = rmfield(ft_appendfreq(acfg,LowChans{:}),'cfg');
    
    for m = 1:length(Deci.Plot.CFC.methods)
        for Conditions = 1:length(Deci.Plot.Conditions)
            maxt = max(sum(ismember(freq.condinfo{2},Deci.Plot.Conditions{Conditions}),2));
            trl = sum(ismember(freq.condinfo{2},Deci.Plot.Conditions{Conditions}),2) == maxt;
            
            TrialCount{subject_list,Conditions} = length(find(trl));
            
            High = HighChans;
            High.fourierspctrm = HighChans.fourierspctrm(trl,:,:,:);
            Low = LowChans;
            Low.fourierspctrm = LowChans.fourierspctrm(trl,:,:,:);
            
            
            Deci.Plot.CFC.method = Deci.Plot.CFC.methods{m};
            Subjects{subject_list,Conditions,m} =  ft_singlecfc(Deci.Plot.CFC,Low,High);
            
        end
        
    end
    clear HighChans
    clear LowChans
end


if Deci.Plot.GA
    
    for conds = 1:size(Subjects,2)
        for m = 1:length(Deci.Plot.CFC.methods)
            facfg.parameter =  'crsspctrm';
            
            CFCData{1,conds,m} = Subjects{1,conds,m};
            CFCData{1,conds,m}.crsspctrm =  mean([cell2mat(cellfun(@(c) c.crsspctrm,Subjects(:,conds,m),'UniformOutput',false))],1);
            
            TotalCount{1,conds} = mean([TrialCount{:,conds}]);
        end
    end
    
     Deci.SubjectList = {'Group 1'};
else
    CFCData = Subjects;
    TotalCount = TrialCount;
    
   

    
end
clear Subjects;

if ~isempty(Deci.Plot.Math)
    
    for m = 1:length(Deci.Plot.CFC.methods)
        for cond = 1:length(Deci.Plot.Math.Form)
            
            
            for subj = 1:size(CFCData,1)
                
                
                operation = str2fun(['@(x)' regexprep(Deci.Plot.Math.Form{cond},'x(\d*)','x{$1}')]);
                
                MathData(subj) = CFCData(subj,1,m);
                MathData{subj}.crsspctrm = [feval(operation, cellfun(@(c) c.crsspctrm,CFCData(subj,1:length(Deci.Plot.Conditions)+cond-1,m),'UniformOutput',false))];
                
                
            end
            CFCData(:,length(Deci.Plot.Conditions)+cond,m) = MathData;
            TotalCount(:,length(Deci.Plot.Conditions)+cond) = num2cell(nan(size(CFCData,1),1));
        end
        
        
    end
    
    
end


for cond = 1:length(Deci.Plot.Draw)
    for met = 1:size(CFCData,3)
        for subj = 1:size(CFCData,1)
            
            if Deci.Plot.CFC.Square
                CFCsquare(subj,met) = figure;
            end
            
            if Deci.Plot.CFC.Topo
                CFCtopo(subj,met)  = figure;
            end
            
            if Deci.Plot.CFC.Hist
                CFChist(subj,met)  = figure;
            end
            
            for subcond = 1:length(Deci.Plot.Draw{cond})
                
                if Deci.Plot.CFC.Topo
                    set(0, 'CurrentFigure', CFCtopo(subj,met) )
                    CFCtopo(subj,met).Visible = 'on';
                    
                    ctopo(subj,subcond,met)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                    
                    
                    Cross.labelcmb = CombVec(CFCData{subj,subcond,1}.labellow',CFCData{subj,subcond,1}.labelhigh')';
                    Cross.freq = CFCData{subj,subcond,1}.freqlow;
                    Cross.cohspctrm = CFCData{subj,subcond,met}.crsspctrm(:,:,2);
                    Cross.dimord = CFCData{subj,subcond,1}.dimord;
                    
                    Crosscfg =[];
                    Crosscfg.foi = [min([Deci.Plot.CFC.freqhigh Deci.Plot.CFC.freqlow]) max([Deci.Plot.CFC.freqhigh Deci.Plot.CFC.freqlow])];
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
            NewCrossCfg.CLim = minmax(CLim(:)');
            ft_topoplotCC(NewCrossCfg ,NewCross);
            
        end
        
    end

end