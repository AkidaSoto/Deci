function Plottor_Conn(Deci,Params)

%% Load
cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';

disp('----------------------');
display(' ')

display(['Plotting Conn']);

for ConnList = 1:length(Params.List)
    clear Subjects
    
    for  subject_list = 1:length(Deci.SubjectList)
        
        display(['Loading Plottor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
        
        for Conditions = 1:length(Deci.Plot.CondTitle)
            
            load([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.CondTitle{Conditions} filesep strjoin(Params.List{ConnList},'_') '.mat'],'conn');
            
            %toi = round(freq.time,4) >= Tois(1) & round(freq.time,4) <= Tois(2);
            
            if isfield(conn,'trllen')
                trllen(subject_list,Conditions) = conn.trllen;
            else
                trllen(subject_list,Conditions) = nan;
            end
            
            if isfield(conn,'lockers')
                lockers(subject_list,Conditions,:) = conn.lockers;
            else
                lockers(subject_list,Conditions,:) = [];
            end
            
            conn.dimord = 'freq_freq_time';
            %$conn.freq = conn.freqlow;
            conn.label = conn.chanlow;
            
            Subjects{subject_list,Conditions} = conn;
        end
        
    end
    
    %% Baseline Correction
    display(' ');
    display(['Using Lock: ' Deci.Plot.Lock]);
    display(['Using Ref: ' Deci.Plot.BslRef ' at times ' strrep(regexprep(num2str(Deci.Plot.Freq.Bsl),' +',' '),' ','-')]);
    
    if Params.Bsl
        
        for Conditions = 1:size(Subjects,2)
            for subject_list = 1:size(Subjects,1)
                
                if ~strcmpi(Deci.Plot.BslRef,Deci.Plot.Lock)
                    
                    load([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.BslRef filesep Deci.Plot.CondTitle{Conditions} filesep strjoin(Params.List{ConnList},'_') '.mat'],'conn');
                    
                    ccfg.latency = Deci.Plot.Freq.Bsl;
                    ccfg.avgovertime = 'yes';
                    
                    Bsl{subject_list,Conditions} = conn;
                    
                    toi1 = round(Bsl{subject_list,Conditions}.time,4) >= Deci.Plot.Freq.Bsl(1) & round(Bsl{subject_list,Conditions}.time,4) <= Deci.Plot.Freq.Bsl(2);
                    
                    Bsl{subject_list,Conditions}.param =  Bsl{subject_list,Conditions}.param(:,:,toi1);
                    Bsl{subject_list,Conditions}.time = Bsl{subject_list,Conditions}.time(toi1);
                    Bsl{subject_list,Conditions} = ft_selectdata(ccfg, Bsl{subject_list,Conditions});
                    Bsl{subject_list,Conditions}.param = repmat(Bsl{subject_list,Conditions}.param,[1 1 size(Subjects{subject_list,Conditions}.param ,3)]);
                    
                else
                    ccfg.latency = Deci.Plot.Freq.Bsl;
                    ccfg.avgovertime = 'yes';
                    
                    toi1 = round(Subjects{subject_list,Conditions}.time,4) >= Deci.Plot.Freq.Bsl(1) & round(Subjects{subject_list,Conditions}.time,4) <= Deci.Plot.Freq.Bsl(2);
                    Bsl{subject_list,Conditions} =Subjects{subject_list,Conditions};
                    Bsl{subject_list,Conditions}.param =  Bsl{subject_list,Conditions}.param(:,:,toi1);
                    Bsl{subject_list,Conditions}.time = Bsl{subject_list,Conditions}.time(toi1);
                    Bsl{subject_list,Conditions} = ft_selectdata(ccfg, Bsl{subject_list,Conditions});
                    Bsl{subject_list,Conditions}.param = repmat(Bsl{subject_list,Conditions}.param,[1 1 size(Subjects{subject_list,Conditions}.param ,3)]);
                end
                
                switch Deci.Plot.Freq.BslType
                    case 'none'
                    case 'absolute'
                        Subjects{subject_list,Conditions}.param =  Subjects{subject_list,Conditions}.param - Bsl{subject_list,Conditions}.param;
                    case 'relative'
                        Subjects{subject_list,Conditions}.param=  Subjects{subject_list,Conditions}.param ./ Bsl{subject_list,Conditions}.param;
                    case 'relchange'
                        Subjects{subject_list,Conditions}.param = ( Subjects{subject_list,Conditions}.param - Bsl{subject_list,Conditions}.param) ./ Bsl{subject_list,Conditions}.param;
                    case 'db'
                        Subjects{subject_list,Conditions}.param = 10*log10( Subjects{subject_list,Conditions}.param ./ Bsl{subject_list,Conditions}.param);
                end
                
                %Subjects{subject_list,Conditions}.time = Subjects{subject_list,Conditions}.time(toi);
                %Subjects{subject_list,Conditions}.param = Subjects{subject_list,Conditions}.param(:,:,toi);
                
            end
        end
    end
    
    %% Math
    if ~isempty(Deci.Plot.Math)
        
        display(' ')
        display(['Doing ' num2str(length(Deci.Plot.Math)) ' Maths'] )
        
        for conds = 1:length(Deci.Plot.Math)
            for subj = 1:size(Subjects,1)
                scfg.parameter = 'param';
                scfg.operation = Deci.Plot.Math{conds};
                evalc('MathData{subj} = ft_math(scfg,Subjects{subj,:})');
            end
            
            Subjects(:,size(Subjects,2)+1) = MathData;
        end
    end
    
    %% Data Management
    if size(Subjects,1) == 1
        Deci.Plot.GrandAverage = false;
    end
    
    if any(~isnan(trllen))
        trlstd = nanstd(trllen,[],1);
        trllen = nanmean(trllen,1);
    end
    
    if any(~isnan(lockers))
        lockersstd = nanstd(lockers,[],1);
        lockers = nanmean(lockers,1);
    end
    
    for conds = 1:size(Subjects,2)
        
        if Deci.Plot.GrandAverage
            facfg.parameter =  'param';
            facfg.type = 'mean';
            facfg.keepindividual = 'yes';
            evalc('ConnData{1,conds} = ft_timelockgrandaverage(facfg,Subjects{:,conds});');
            ConnData{1,conds}.param = ConnData{1,conds}.individual;
            ConnData{1,conds} = rmfield(ConnData{1,conds},'cfg');
            
            WiretData{1,conds} = ConnData{1,conds};
            
            WiretData{1,conds}.avg = squeeze(WiretData{1,conds}.individual);
            WiretData{1,conds}.dimord = 'subj_time';
            
        else
            Params.Stat = false;
            ConnData(:,conds) = Subjects(:,conds);
        end
        
                    
        if Params.Wire
            facfg.parameter =  'param';
            facfg.type = 'mean';
            facfg.keepindividual = 'no';
            evalc('wiredata{1,conds} = ft_timelockgrandaverage(facfg,Subjects{:,conds});');
            wiredata{1,conds} = rmfield(wiredata{1,conds},'cfg');
        end
        
        
        for subj = 1:size(ConnData(:,conds),1)
   
            if Params.Bar.do
                tcfg = [];
                tcfg.nanmean = Deci.Plot.nanmean;
                tcfg.latency = Deci.Plot.Freq.Bar.Toi;
                
                bardata{subj,conds} = ft_selectdata(tcfg,ConnData{subj,conds});
                
                tcfg.avgovertime = 'yes';
                
                bartdata{subj,conds} = ft_selectdata(tcfg,bardata{subj,conds});
            end
            
        end
    end
    
    %% Stat
    
    if Params.Stat
        
        display(' ')
        display('Calculating Statistics')
        
        neighbours       = load('easycap_neighbours','neighbours');
        Deci.Plot.Stat.neighbours = neighbours.neighbours;
        Deci.Plot.Stat.ivar = 1;
        
        Deci.Plot.Stat.parameter = 'avg';
        
        for conds = 1:length(Deci.Plot.Draw)
            design = [];
            
            switch Deci.Plot.Stat.Comp
                case 'Conditions'
                    if length(Deci.Plot.Draw{conds}) ~= 1
                        Deci.Plot.Stat.uvar = 2;
                        
                        for subcond = 1:length(Deci.Plot.Draw{conds})
                            for subj = 1:size(Subjects,1)
                                design(1,subj+size(Subjects,1)*[subcond-1]) =  subcond;
                                design(2,subj+size(Subjects,1)*[subcond-1]) = subj;
                            end
                        end
                        
                        Deci.Plot.Stat.design = design;
                        
                        if length(Deci.Plot.Draw{conds}) > 2
                            Deci.Plot.Stat.tail = 1;
                            Deci.Plot.Stat.statistic = 'depsamplesFmultivariate';
                            Deci.Plot.Stat.clustertail      = 1;
                        else
                            Deci.Plot.Stat.statistic = 'depsamplesT';
                            Deci.Plot.Stat.tail = 0;
                            Deci.Plot.Stat.clustertail      = 0;
                        end
                        if  Params.Wire
                            
                            [wirestat{conds}] = ft_timelockstatistics(Deci.Plot.Stat, WiretData{:,Deci.Plot.Draw{conds}});
                        end
                        
                        if Params.Bar.do
                            [barstat{conds}] = ft_freqstatistics(Deci.Plot.Stat, bartdata{:,Deci.Plot.Draw{conds}});
                        end
                    else
                        
                        if  Params.Wire
                            [wirestat{conds}.mask] = ones(size(ConnData{:,Deci.Plot.Draw{conds}}.param(1,:)));
                        end
                        
                        if Params.do
                            [barstat{conds}.mask] = ones(size(bartdata{Deci.Plot.Draw{conds}}.param(1,:)));
                        end
                        
                    end
                case 'Bsl'
                    Deci.Plot.Stat.tail = 0;
                    Deci.Plot.Stat.statistic = 'indepsamplesT';
                    
                    Deci.Plot.Stat.design = ones(size(Subjects,1));
                    
                    if  Params.Wire
                        
                        allwiretdata = ft_timelockgrandaverage(struct('parameter','param'),ConnData{:,Deci.Plot.Draw{conds}});
                        allwiretdata.avg = squeeze(allwiretdata.avg);
                        allwiretdata.dimord = 'subj_time';
                        [wirestat{conds}] = ft_timelockstatistics(Deci.Plot.Stat, allwiretdata);
                    end
                    
                    if Params.Bar.do
                        allbartdata = ft_freqgrandaverage(struct('parameter','param'),bartdata{:,Deci.Plot.Draw{conds}});
                        allbartdata.dimord = 'rpt_freq_freq_time';
                        [barstat{conds}] = ft_freqstatistics(Deci.Plot.Stat, allbartdata);
                    end
                    
            end
            
        end
    end
    
    %% Plot
    
    if Deci.Plot.GrandAverage
        SubjectList = {'Group Average'};
        
    end
    
    
    
    for cond = 1:length(Deci.Plot.Draw)
        clear cirky subby
        for subj = 1:size(ConnData,1)
            
            if Deci.Plot.Extra.Conn.Wire
                wire(subj)  = figure;
            end
            
            for subcond = 1:length(Deci.Plot.Draw{cond})
                
                
                if Deci.Plot.Extra.Conn.Wire
                    set(0, 'CurrentFigure', wire(subj) )
                    wire(subj).Visible = 'on';
                    
                    top = squeeze(nanmean(nanmean(nanmean(ConnData{subj,Deci.Plot.Draw{cond}(subcond)}.param,1),2),3)) + squeeze(nanstd(nanmean(nanmean(ConnData{subj,Deci.Plot.Draw{cond}(subcond)}.param,2),3),[],1))/sqrt(size(ConnData{subj,Deci.Plot.Draw{cond}(subcond)}.param,1));
                    bot = squeeze(nanmean(nanmean(nanmean(ConnData{subj,Deci.Plot.Draw{cond}(subcond)}.param,1),2),3)) - squeeze(nanstd(nanmean(nanmean(ConnData{subj,Deci.Plot.Draw{cond}(subcond)}.param,2),3),[],1))/sqrt(size(ConnData{subj,Deci.Plot.Draw{cond}(subcond)}.param,1));
                    
                    if Deci.Plot.GrandAverage
                        pgon = polyshape([ConnData{subj,Deci.Plot.Draw{cond}(subcond)}.time fliplr(ConnData{subj,Deci.Plot.Draw{cond}(subcond)}.time)],[top' fliplr(bot')],'Simplify', false);
                        b(subcond) = plot(pgon,'HandleVisibility','off');
                        hold on
                        b(subcond).EdgeAlpha = 0;
                        b(subcond).FaceAlpha = .15;
                    end
                    
                    if Params.Stat
                       wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.mask = permute(wirestat{cond}.mask,[2 3 1]);
                    end
                    
                    
                end
                
            end
            
            if Deci.Plot.Extra.Conn.Wire
                set(0, 'CurrentFigure', wire(subj) )
                wire(subj).Visible = 'on';
                
                pcfg = cfg;
                pcfg.clim = 'maxmin';
                
                if Params.Stat
                    pcfg.maskparameter ='mask';
                    
                end
                
                pcfg.ylim = ylim;
                pcfg.graphcolor = lines;
                pcfg.linewidth = 1;
                pcfg.parameter = 'avg';
                
                %ConnData(subj,Deci.Plot.Draw{cond})= cellfun(@(c) rmfield(rmfield(c,'chanlow'),'chanhigh'), ConnData(subj,Deci.Plot.Draw{cond}),'un',0);
                %ConnData(subj,Deci.Plot.Draw{cond})= cellfun(@(c) rmfield(rmfield(c,'freqlow'),'freqhigh'), ConnData(subj,Deci.Plot.Draw{cond}),'un',0);
                
                pcfg.maskstyle = 'box';
                %ConnData(subj,Deci.Plot.Draw{cond}) = cellfun(@(c) setfield(c,'dimord','chan_freq_time'), ConnData(subj,Deci.Plot.Draw{cond}),'un',0);
                ft_singleplotER(pcfg,wiredata{subj,Deci.Plot.Draw{cond}});
                
                if Deci.Plot.GrandAverage
                    arrayfun(@(c) uistack(c,'top'),b);
                    clear b
                end
                
                hold on
                axis tight
                plot([ConnData{cond}.time(1), ConnData{cond}.time(end)], [0 0], 'k--') % hor. line
                plot([0 0], ylim, 'k--') % vert. l
                
                if Params.Stat
                    boxes = wire(subj).Children(2).Children.findobj('Type','Patch');
                    for bb = 1:length(boxes)
                        if ~isempty(boxes)
                            boxes(bb).FaceAlpha = .35;
                            uistack(boxes(bb),'bottom')
                            boxes(bb).HandleVisibility = 'off';
                        end
                    end
                end
                
            if max(Deci.Plot.Draw{cond}) <= size(trllen,2)
                legend(arrayfun(@(a,b) [ a{1} ' (' num2str(b) ')'] ,Deci.Plot.Subtitle{cond},trllen(subj,Deci.Plot.Draw{cond}),'UniformOutput',false));
            else
                legend([Deci.Plot.Subtitle{cond}{subcond}]);
            end
                
                title([Deci.SubjectList{subj} ' ' strjoin(Params.List{ConnList},' ')], 'Interpreter', 'none');
                set(legend, 'Interpreter', 'none')
                xlim([ConnData{cond}.time(1) ConnData{cond}.time(end)])
                xlabel('Time');
                
                
                for subcond = 1:length(Deci.Plot.Draw{cond})
                    if Deci.Plot.Draw{cond}(subcond) <= size(lockers,2)
                        xlims = xlim;
                        ylims = ylim;
                        
                        for locks = 1:length([lockers(subj,Deci.Plot.Draw{cond}(subcond),:)])
                            hold on
                            
                            locktime = [lockers(subj,Deci.Plot.Draw{cond}(subcond),locks)/1000];
                            lockstd = [lockersstd(subj,Deci.Plot.Draw{cond}(subcond),locks)/1000];
                            
                            if locktime > xlims(1) && locktime < xlims(2)
                                plotlock = plot([locktime locktime], ylims,'LineWidth',2,'Color','k','LineStyle','--','HandleVisibility','off');
                                plotlock.Color(4) = .2;
                                
                                if Deci.Plot.GrandAverage
                                    
                                    if [locktime - lockstd] < xlims(1)
                                        lockpstd(1) = xlims(1);
                                    else
                                        lockpstd(1) = [locktime - lockstd];
                                    end
                                    
                                    if [locktime + lockstd] > xlims(2)
                                        lockpstd(2) = xlims(2);
                                    else
                                        lockpstd(2) = [locktime + lockstd];
                                    end
                                    
                                    lockpgon = polyshape([lockpstd fliplr(lockpstd)],sort([ylims ylims]),'Simplify', false);
                                    lockb = plot(lockpgon,'HandleVisibility','off');
                                    hold on
                                    lockb.EdgeAlpha = 0;
                                    lockb.FaceAlpha = .15;
                                    lockb.FaceColor = plotlock.Color;
                                    
                                    arrayfun(@(c) uistack(c,'top'),lockb);
                                    clear lockb
                                    
                                end
                            end
                            ylim(ylims)
                        end
                    end
                end
                
            end
            
            if Deci.Plot.Extra.Conn.Bar.do
                barfig(subj)  = figure;
                set(0, 'CurrentFigure', barfig(subj) )
                barfig(subj).Visible = 'on';
                
                CleanBars(mean(cell2mat(arrayfun(@(c) nanmean(nanmean(nanmean(c.param,2),3),4),[bardata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)),1), ...
                    nanstd(cell2mat(arrayfun(@(c) nanmean(nanmean(nanmean(c.param,2),3),4),[bardata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)),[],1) ...
                    /sqrt(size(cell2mat(arrayfun(@(c) nanmean(nanmean(nanmean(c.param,2),3),4),[bardata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)),1)));
                
                l = legend(Deci.Plot.Subtitle{cond});
                title([Deci.SubjectList{subj} ' ' Deci.Plot.Freq.Type ' ' Deci.Plot.Title{cond} ' Wire'])
                set(l, 'Interpreter', 'none')
                if Params.Stat
                    hold on
                    
                    if barstat{cond}.prob < Deci.Plot.Stat.alpha
                        plot([.75 1.25],[max(ylim) max(ylim)]*.90, '-k', 'LineWidth',2);
                        plot([1],[max(ylim)]*.95, '*k');
                    end
                end
            end
            
        end
    end
    
end