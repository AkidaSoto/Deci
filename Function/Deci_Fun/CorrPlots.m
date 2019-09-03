function CorrPlots(Deci,params)

for freq = 1:length(params.Freq)
    
    for param = 1:length(params.Variable)
        
        
        for Cond = 1:length(Deci.Analysis.CondTitle)
            for subject_list = 1:length(Deci.SubjectList)
                
                
                
                Channel = CleanDir([Deci.Folder.Analysis filesep 'Extra' filesep 'Corr' filesep Deci.Plot.Extra.Corr.Freq{freq} '_' Deci.Plot.Extra.Corr.Variable{param} filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Analysis.CondTitle{Cond}]);
                
                for Choi = 1:length(Channel)
                    
                    R = []; P = []; Beta = [];
                    
                    load([Deci.Folder.Analysis filesep 'Extra' filesep 'Corr' filesep Deci.Plot.Extra.Corr.Freq{freq} '_' Deci.Plot.Extra.Corr.Variable{param} filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Analysis.CondTitle{Cond} filesep Channel{Choi}],'R','P');
                    
                    RChan{Choi} = R;
                    PChan{Choi} = P;
                    
                    if any(any(R.powspctrm > 1))
                      k = 0;  
                    end
                    
                end
                
                
                
                RSub{subject_list} = rmfield(ft_appendfreq(struct('parameter','powspctrm','appenddim','chan'),RChan{:}),'cfg');
                RSub{subject_list}.trialinfo = RChan{1}.trialinfo;
                PSub{subject_list} = rmfield(ft_appendfreq(struct('parameter','powspctrm','appenddim','chan'),PChan{:}),'cfg');
                PSub{subject_list}.trialinfo = PChan{1}.trialinfo;
                
                
                
                if strcmpi(RSub{subject_list}.dimord,'rpt_chan_freq_time')
                    
                    if isfield(RSub{subject_list},'trialinfo') && ismember(RSub{subject_list}.trialinfo,params.Varnum)
                        RSub{subject_list} = rmfield(ft_selectdata(struct('trials',find(params.Varnum == RSub{subject_list}.trialinfo),'avgoverrpt','yes'),RSub{subject_list}),'cfg');
                        PSub{subject_list} = rmfield(ft_selectdata(struct('trials',find(params.Varnum == P.trialinfo),'avgoverrpt','yes'),PSub{subject_list}),'cfg');
                    else
                        RSub{subject_list} = rmfield(ft_selectdata(struct('avgoverrpt','yes'),RSub{subject_list}),'cfg');
                        PSub{subject_list} = rmfield(ft_selectdata(struct('avgoverrpt','yes'),PSub{subject_list}),'cfg');
                    end
                end
                
            end
            
            RCond{Cond} = rmfield(ft_freqgrandaverage(struct('parameter','powspctrm','keepindividual','no'),RSub{:}),'cfg');
            PCond{Cond} = rmfield(ft_freqgrandaverage(struct('parameter','powspctrm','keepindividual','no'),PSub{:}),'cfg');
            
        end
        
        %% Math
        if ~isempty(Deci.Plot.Math)
            for conds = 1:length(Deci.Plot.Math)
                scfg.parameter = 'powspctrm';
                scfg.operation = Deci.Plot.Math{conds};
                MathData = ft_math(scfg,RCond{:});
                RCond{size(RCond,1)+1} = MathData;
            end
        end
        
        %%
        
        for cond = 1:length(Deci.Plot.Draw)
            
            for subj = 1:length(Deci.SubjectList)
                
                if Deci.Plot.Extra.Corr.Square
                    Corrsquare(subj) = figure;
                    suptitle([Deci.Plot.Title{cond} ' '  Deci.Plot.Extra.Corr.Freq{freq} '_' Deci.Plot.Extra.Corr.Variable{param}]);
                end
                
                % Topo plot with lines for connections, Time has a Dial
                if Deci.Plot.Extra.Corr.Topo
                    Corrtopo(subj)  = figure;
                    suptitle([Deci.Plot.Title{cond} ' '  Deci.Plot.Extra.Corr.Freq{freq} '_' Deci.Plot.Extra.Corr.Variable{param}]);
                    
                end
                
                if Deci.Plot.Extra.Corr.Wire
                    Corrwire(subj)  = figure;
                end
                
                for subcond = 1:length(Deci.Plot.Draw{cond})
                    
                    if Deci.Plot.Extra.Corr.Topo
                        set(0, 'CurrentFigure', Corrtopo(subj) )
                        Corrtopo(subj).Visible = 'on';
                        
                        ctopo(subj,subcond)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                        
                        Crosscfg =[];
                        Crosscfg.layout = Deci.Layout.Noeye;
                        Crosscfg.toi = 1;
                        Crosscfg.foi = 1;
                        
                        title([Deci.SubjectList{subj} ' Cond '  num2str(cond)],'Interpreter', 'none');
                        
                        %RCond{subj,subcond}.mask = ft_selectdata(struct('latency',RCond{subj,subcond}.time(1),'frequency',RCond{subj,subcond}.freq(1)),PCond{subj,subcond});
                        %RCond{subj,subcond}.mask = RCond{subj,subcond}.mask.powspctrm;
                        
                        ft_topoplotER(Crosscfg,ft_selectdata(struct('latency',RCond{subj,subcond}.time(1),'frequency',RCond{subj,subcond}.freq(1)),RCond{subj,subcond}));
                        ctopo(subj,subcond).UserData = {RCond{subj,subcond},Crosscfg};
                    end
                    
                    if Deci.Plot.Extra.Corr.Square
                        
                        set(0, 'CurrentFigure', Corrsquare(subj) )
                        Corrsquare(subj).Visible = 'on';
                        csquare(subj,subcond)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                        
                        Crosscfg =[];
                        Crosscfg.layout = Deci.Layout.Noeye;
                        Crosscfg.choi = 1;
                        
                        title([Deci.SubjectList{subj} ' Cond '  num2str(cond)],'Interpreter', 'none');
                        
                        ft_singleplotTFR(Crosscfg,ft_selectdata(struct('channel',RCond{subj,subcond}.label(1)),RCond{subj,subcond}));
                        
                        csquare(subj,subcond).UserData = {RCond{subj,subcond},Crosscfg};
                        colorbar;
                    end
                    
                    if Deci.Plot.Extra.Corr.Wire
                        set(0, 'CurrentFigure', CFCwire(subj) )
                        CFCwire(subj).Visible = 'on';
                        
                        Crosscfg.chan = 1;
                        Crosscfg.foi = 1;
                        
                        cwire(subj,subcond)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                        
                        top = squeeze(mean(CFCData{subj,Deci.Plot.Draw{cond}(subcond)}.crsspctrm(chanl,chanh,:),1)) + squeeze(mean(CFCSem{subj,Deci.Plot.Draw{cond}(subcond)}.crsspctrm(chanl,chanh,:),1));
                        bot = squeeze(mean(CFCData{subj,Deci.Plot.Draw{cond}(subcond)}.crsspctrm(chanl,chanh,:),1)) - squeeze(mean(CFCSem{subj,Deci.Plot.Draw{cond}(subcond)}.crsspctrm(chanl,chanh,:),1));
                        
                        pgon = polyshape([CFCData{subj,Deci.Plot.Draw{cond}(subcond)}.timelow fliplr(CFCData{subj,Deci.Plot.Draw{cond}(subcond)}.timelow)],[top' fliplr(bot')],'Simplify', false);
                        b = plot(pgon,'HandleVisibility','off');
                        hold on
                        b.EdgeAlpha = 0;
                        b.FaceAlpha = .15;
                        
                        h =  plot(CFCData{subj,subcond}.timelow,squeeze(CFCData{subj,subcond}.crsspctrm(chanl,chanh,:)));
                        h.Color = b.FaceColor;
                        h.LineWidth = 1;
                        
                        hold on
                    end
                end
                
                
                if Deci.Plot.Extra.Corr.Wire
                    for WireLim = 1:length(cwire(subj,:))
                        Vect = CombVec(CFCData{subj,subcond}.labellow',CFCData{subj,subcond}.labelhigh')';
                        legend(cwire(subj,WireLim),arrayfun(@(c) [Vect{c,1} Vect{c,2}],1:size(Vect,1),'un',0));
                        cwire(subj,WireLim).YLim = [min([cwire(subj,:).YLim]) max([cwire(subj,:).YLim])];
                    end
                end
                
                if Deci.Plot.Extra.Corr.Topo
                    
                    set(0, 'CurrentFigure', Corrtopo(subj) )
                    
                    text = uicontrol('style','text','position',[225 75 100 25],'String',['Time of Interest: '  num2str(RCond{subj,subcond}.time(1))]);
                    uicontrol('style','pushbutton','position',[175 75 45 25],'String','<','callback',{@ChangeDim,ctopo(subj,:),-1,'toi',text,'Time of Interest: ','Topo'})
                    uicontrol('style','pushbutton','position',[325 75 45 25],'String','>','callback',{@ChangeDim,ctopo(subj,:),+1,'toi',text,'Time of Interest: ','Topo'})
                    
                    
                    text = uicontrol('style','text','position',[225+100+25+50+50 75 100 25],'String',['Freq of Interest: '  num2str(RCond{subj,subcond}.freq(1))]);
                    uicontrol('style','pushbutton','position',[175+100+25+50+50 75 45 25],'String','<','callback',{@ChangeDim,ctopo(subj,:),-1,'foi',text,'Freq of Interest: ','Topo'})
                    uicontrol('style','pushbutton','position',[325+100+25+50+50 75 45 25],'String','>','callback',{@ChangeDim,ctopo(subj,:),+1,'foi',text,'Freq of Interest: ','Topo'})
                    
                    ChangeDim([],[],ctopo(subj,:),-1,'toi',text,'Time of Interest: ','Topo')
                end
                
                
                if Deci.Plot.Extra.Corr.Square
                    
                    set(0, 'CurrentFigure', Corrsquare(subj) )
                    text = uicontrol('style','pushbutton','position',[225 75 125 25],'String',['Channel of Interest: '  RCond{subj,subcond}.label{1}],'callback',{@ChangeDim,csquare(subj,:),0,'choi',text,'Channel of Interest: ','Square'});
                    
                    ChangeDim([],[],csquare(subj,:),0,'choi',text,'Channel of Interest: ','Square')
                end
                
            end
            
        end
    end
    
    
    
end



    function ChangeDim(popup,event,Axes,Change,Dim,TextUI,Text,Type)
        
        if strcmpi(Dim,'choi')
            if ~isempty(popup)
                fakeUI = figure;
                fakeUI.UserData = Axes(1).UserData{1}.label(Axes(1).UserData{2}.choi);
                fakeUI.Visible =  'off';
                select_labels(fakeUI,[],Axes(1).UserData{1}.label);
                waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
                choi = find(ismember(Axes(1).UserData{1}.label,fakeUI.UserData));
                close(fakeUI);
            else
                choi = 1;
            end
            
        end
        
        for Axe =  1:length(Axes)
            set(0, 'CurrentFigure', Axes(Axe).Parent )
            set(Axes(Axe).Parent, 'currentaxes', Axes(Axe))
            
            NewCrossCfg = Axes(Axe).UserData{2};
            
            switch Dim
                
                case 'toi'
                    NewCrossCfg.toi = NewCrossCfg.toi + Change;
                    if NewCrossCfg.toi < 1 || NewCrossCfg.toi > length(Axes(Axe).UserData{1}.time)
                        break
                    end
                    Axes(Axe).UserData{2}.toi = Axes(Axe).UserData{2}.toi + Change;
                    TextUI.String = [Text num2str(Axes(Axe).UserData{1}.time(NewCrossCfg.toi) )];
                    
                    
                case 'foi'
                    NewCrossCfg.foi = NewCrossCfg.foi + Change;
                    if NewCrossCfg.foi < 1 || NewCrossCfg.foi > length(Axes(Axe).UserData{1}.freq)
                        break
                    end
                    Axes(Axe).UserData{2}.foi = Axes(Axe).UserData{2}.foi + Change;
                    TextUI.String = [Text num2str(Axes(Axe).UserData{1}.freq(NewCrossCfg.foi) )];
                case 'choi'
                    NewCrossCfg.choi = choi;
                    Axes(Axe).UserData{2}.choi = NewCrossCfg.choi;
                    popup.String = [Text Axes(Axe).UserData{1}.label{NewCrossCfg.choi}];
            end
            
            switch Type
                case 'Topo'
                    ft_topoplotTFR(NewCrossCfg,ft_selectdata(struct('latency',Axes(Axe).UserData{1}.time(NewCrossCfg.toi),'frequency',Axes(Axe).UserData{1}.freq(NewCrossCfg.foi)),Axes(Axe).UserData{1}));
                case 'Square'
                    ft_singleplotTFR(NewCrossCfg,ft_selectdata(struct('channel',Axes(Axe).UserData{1}.label(NewCrossCfg.choi)),Axes(Axe).UserData{1}));
            end
        end
        
        for Axe =  1:length(Axes)
            
            switch Type
                case 'Topo'
                    Axes(Axe).CLim = [-1*max(arrayfun(@(c) nanmax(nanmax(abs(c.Children(end).CData))),Axes,'UniformOutput',1)) max(arrayfun(@(c) nanmax(nanmax(abs(c.Children(end).CData))),Axes,'UniformOutput',1))];
                case 'Square'
                    Axes(Axe).CLim = [-1*max(arrayfun(@(c) nanmax(nanmax(abs(c.Children.CData))),Axes,'UniformOutput',1)) max(arrayfun(@(c) nanmax(nanmax(abs(c.Children.CData))),Axes,'UniformOutput',1)) ];
            end
            
            
        end
        
    end
end


