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
    
    
    Current = Params.List{ConnList};
    chanl = Current{1};
    chanh = Current{2};
    freqlow = Current{3};
    freqhigh = Current{4};
    conntype = Current{5};
    
    if ischar(conntype)
        conntype = {conntype};
    end
    
    %% do plots for 1 conoi at a time
    for conoi = 1:length(conntype)
        clear Subjects
        
        
        if ischar(chanl)
            chanl = {chanl};
        end
        
        if ischar(chanh)
            chanh = {chanh};
        end
        
        if ischar(freqlow)
            freqlow = {freqlow};
        end
        
        if ischar(freqhigh)
            freqhigh = {freqhigh};
        end
        
        %% load 1 full conoi
        
        for  subject_list = 1:length(Deci.SubjectList)
            
            display(['Loading Plottor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
            for Conditions = 1:length(Deci.Plot.CondTitle)
                
                for choil = 1:length(chanl)
                    
                    for choih = 1:length(chanh)
                        
                        for foil = 1:length(freqlow)
                            
                            for foih = 1:length(freqhigh)
                                
                                connfile = strjoin([chanl(choil) chanh(choih) freqlow(foil) freqhigh(foih) conntype(conoi)],'_');
                                
                                load([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.CondTitle{Conditions} filesep connfile  '.mat'],'conn');
                                
                                sub_freqhigh{foih} = conn;
                                
                            end
                            
                            sub_freqlow{foil} = sub_freqhigh{foih};
                            param  = cellfun(@(c) c.param,sub_freqhigh,'UniformOutput',false);
                            sub_freqlow{foil}.param = cat(2,param{:});
                            
                            if ~strcmpi(conntype{conoi},'ispc')
                                freq_freqhigh = cellfun(@(c) c.freqhigh,sub_freqhigh,'UniformOutput',false);
                            else
                                freq_freqhigh = cellfun(@(c) c.freqlow,sub_freqhigh,'UniformOutput',false);
                            end
                            
                            sub_freqlow{foil}.freqhigh = cat(2,freq_freqhigh{:});
                            clear sub_freqhigh
                        end
                        
                        sub_chanh{choih} = sub_freqlow{foil};
                        param  = cellfun(@(c) c.param,sub_freqlow,'UniformOutput',false);
                        sub_chanh{choih}.param = cat(1,param{:});
                        freq_freqlow = cellfun(@(c) c.freqlow,sub_freqlow,'UniformOutput',false);
                        sub_chanh{choih}.freqlow = cat(1,freq_freqlow{:});
                        clear sub_freqlow
                        
                    end
                    
                    sub_chanl{choil} = sub_chanh{choih};
                    param  = cellfun(@(c) c.param,sub_chanh,'UniformOutput',false);
                    sub_chanl{choil}.param = cat(4,param{:});
                    chan_chanhigh = cellfun(@(c) c.chanhigh,sub_chanh,'UniformOutput',false);
                    sub_chanl{choil}.chanhigh = cat(1,chan_chanhigh{:})';
                    clear sub_chanh
                end
                
                sub_cond{subject_list,Conditions} = sub_chanl{choil};
                param  = cellfun(@(c) c.param,sub_chanl,'UniformOutput',false);
                sub_cond{subject_list,Conditions}.param = cat(5,param{:});
                chan_chanlow = cellfun(@(c) c.chanlow,sub_chanl,'UniformOutput',false);
                sub_cond{subject_list,Conditions}.chanlow = cat(1,chan_chanlow{:})';
                clear sub_chanl
                
                sub_cond{subject_list,Conditions}.param = permute(sub_cond{subject_list,Conditions}.param,[5 4 1 2 3]);
                sub_cond{subject_list,Conditions}.dimord = 'chanlow_chanhigh_freqlow_freqhigh_time';
                
                
                if isfield(sub_cond{subject_list,Conditions},'trllen')
                    trllen(subject_list,Conditions) = sub_cond{subject_list,Conditions}.trllen;
                else
                    trllen(subject_list,Conditions) = nan;
                end
                
                if isfield(sub_cond{subject_list,Conditions},'lockers')
                    LockNum = Deci.Analysis.Locks(ismember(Deci.Analysis.LocksTitle,Deci.Plot.Lock));
                    lockers(subject_list,Conditions,:) = sub_cond{subject_list,Conditions}.lockers - sub_cond{subject_list,Conditions}.lockers(LockNum);
                    
                else
                    lockers(subject_list,Conditions,:) = [];
                end
                
            end
        end
        
        %% Math
        if ~isempty(Deci.Plot.Math)
            
            display(' ')
            display(['Doing ' num2str(length(Deci.Plot.Math)) ' Maths'] )
            
            for conds = 1:length(Deci.Plot.Math)
                for subj = 1:size(sub_cond,1)
                    scfg.parameter = 'param';
                    scfg.operation = Deci.Plot.Math{conds};
                    
                    arginstr = sprintf('x%i,', 1:length([sub_cond{subj,:}]));
                    arginstr = arginstr(1:end-1); % remove the trailing ','
                    eval(sprintf('operation = @(x) %s;',regexprep( scfg.operation,'x(\d*)','x{$1}')));
                    
                    MathData{subj} =  sub_cond{subj,1};
                    MathData{subj}.param = feval(operation, cellfun(@(c) c.param,sub_cond(subj,:),'un',0));
                end
                
                sub_cond(:,size(sub_cond,2)+1) = MathData;
            end
        end
        
        %% Averaging
        
        if size(sub_cond,1) == 1
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
        
        
        for conds = 1:size(sub_cond,2)
            
            if Deci.Plot.GrandAverage
                
                
                Subjs = cellfun(@(c) c.param,sub_cond(:,conds),'UniformOutput',false);
                Subjs = cat(6,Subjs{:});
                %
                
                Subjs = permute(Subjs,[6 1 2 3 4 5]);
                StatsData{1,conds} = sub_cond{1,conds};
                StatsData{1,conds}.param = Subjs;
                StatsData{1,conds}.dimord = 'subj_chanlow_chanhigh_freqlow_freqhigh_time';
                
                Subjs = permute(nanmean(Subjs,1),[2 3 4 5 6 1]);
                ConnData{1,conds} = sub_cond{1,conds};
                ConnData{1,conds}.param = Subjs;
                clear Subjs;
            else
                Deci.Plot.Stat.do = false;
                ConnData(:,conds) = sub_cond(:,conds);
            end
            
            
            %% Data Management
            for subjs = 1:size(ConnData,1)
                
                if Params.FL_FH.do
                    FL_FH{subjs,conds} = ConnData{subjs,conds};
                    FL_FH{subjs,conds}.param = permute(mean(FL_FH{subjs,conds}.param,[1 2 5],'omitnan'),[3 4 1 2 5]);
                    FL_FH{subjs,conds}.dimord = 'freqlow_freqhigh';
                    
                    if Deci.Plot.Stat.do
                        FL_FH_Stat{subjs,conds} = StatsData{subjs,conds};
                        FL_FH_Stat{subjs,conds}.param = permute(mean(FL_FH_Stat{subjs,conds}.param,[2 3 6],'omitnan'),[1 4 5 2 3 6]);
                        FL_FH_Stat{subjs,conds}.dimord = 'subj_freqlow_freqhigh';
                    end
                end
                
                if Params.FL_time.do
                    FL_time{subjs,conds} = ConnData{subjs,conds};
                    FL_time{subjs,conds}.param = permute(mean(FL_time{subjs,conds}.param,[1 2 4],'omitnan'),[3 5 1 2 4]);
                    FL_time{subjs,conds}.dimord = 'freqlow_freqhigh';
                    
                    if Deci.Plot.Stat.do
                        FL_time_Stat{subjs,conds} = StatsData{subjs,conds};
                        FL_time_Stat{subjs,conds}.param = permute(mean(FL_time_Stat{subjs,conds}.param,[2 3 5],'omitnan'),[1 4 6 2 3 5]);
                        FL_time_Stat{subjs,conds}.dimord = 'subj_freqlow_freqhigh';
                    end
                end
                
                if Params.FH_time.do
                    FH_time{subjs,conds} = ConnData{subjs,conds};
                    FH_time{subjs,conds}.param = permute(mean(FH_time{subjs,conds}.param,[1 2 3],'omitnan'),[4 5 1 2 3]);
                    FH_time{subjs,conds}.dimord = 'freqlow_freqhigh';
                    
                    if Deci.Plot.Stat.do
                        FH_time_Stat{subjs,conds} = StatsData{subjs,conds};
                        FH_time_Stat{subjs,conds}.param = permute(mean(FH_time_Stat{subjs,conds}.param,[2 3 4],'omitnan'),[1 5 6 2 3 4]);
                        FH_time_Stat{subjs,conds}.dimord = 'subj_freqlow_freqhigh';
                    end
                end
                
            end
            
        end
        
        %% Stat
        
        %         if Deci.Plot.Stat.do
        %
        %             display(' ')
        %             display('Calculating Statistics')
        %
        %             neighbours       = load('easycap_neighbours','neighbours');
        %             Deci.Plot.Stat.neighbours = neighbours.neighbours;
        %             Deci.Plot.Stat.ivar = 1;
        %
        %             Deci.Plot.Stat.parameter = 'param';
        %
        %             for conds = 1:length(Deci.Plot.Draw)
        %                 design = [];
        %
        %                 switch Deci.Plot.Stat.Comp
        %                     case 'Conditions'
        %                         if length(Deci.Plot.Draw{conds}) ~= 1
        %                             Deci.Plot.Stat.uvar = 2;
        %
        %                             for subcond = 1:length(Deci.Plot.Draw{conds})
        %                                 for subj = 1:size(Subjects,1)
        %                                     design(1,subj+size(Subjects,1)*[subcond-1]) =  subcond;
        %                                     design(2,subj+size(Subjects,1)*[subcond-1]) = subj;
        %                                 end
        %                             end
        %
        %                             design = design';
        %
        %                             Deci.Plot.Stat.design = design;
        %
        %                             if length(Deci.Plot.Draw{conds}) > 2
        %                                 Deci.Plot.Stat.tail = 1;
        %                                 Deci.Plot.Stat.statistic = 'depsamplesFmultivariate';
        %                                 Deci.Plot.Stat.clustertail      = 1;
        %                             else
        %                                 Deci.Plot.Stat.statistic = 'depsamplesT';
        %                                 Deci.Plot.Stat.tail = 0;
        %                                 Deci.Plot.Stat.clustertail      = 0;
        %                             end
        %                             if  Deci.Plot.Wire.do
        %
        %                                 [wirestat{conds}] = ft_freqstatistics(Deci.Plot.Stat, WiretData{:,Deci.Plot.Draw{conds}});
        %                             end
        %
        %                             if Deci.Plot.Bar.do
        %                                 [barstat{conds}] = ft_freqstatistics(Deci.Plot.Stat, bartdata{:,Deci.Plot.Draw{conds}});
        %                             end
        %                         else
        %
        %                             if  Deci.Plot.Wire.do
        %                                 [wirestat{conds}.mask] = ones(size(ConnData{:,Deci.Plot.Draw{conds}}.param(1,:)));
        %                             end
        %
        %                             if Params.do
        %                                 [barstat{conds}.mask] = ones(size(bartdata{Deci.Plot.Draw{conds}}.param(1,:)));
        %                             end
        %
        %                         end
        %                     case 'Bsl'
        %                         Deci.Plot.Stat.tail = 0;
        %                         Deci.Plot.Stat.statistic = 'indepsamplesT';
        %
        %                         Deci.Plot.Stat.design = ones(size(Subjects,1));
        %
        %                         if  Deci.Plot.Wire.do
        %
        %                             allwiretdata = ft_timelockgrandaverage(struct('parameter','param'),ConnData{:,Deci.Plot.Draw{conds}});
        %                             allwiretdata.avg = squeeze(allwiretdata.avg);
        %                             allwiretdata.dimord = 'subj_time';
        %                             [wirestat{conds}] = ft_timelockstatistics(Deci.Plot.Stat, allwiretdata);
        %                         end
        %
        %                         if Deci.Plot.Bar.do
        %                             allbartdata = ft_freqgrandaverage(struct('parameter','param'),bartdata{:,Deci.Plot.Draw{conds}});
        %                             allbartdata.dimord = 'rpt_freq_freq_time';
        %                             [barstat{conds}] = ft_freqstatistics(Deci.Plot.Stat, allbartdata);
        %                         end
        %
        %                 end
        %
        %             end
        %         end
        
        %% Plot
        
        SubjectList = Deci.SubjectList;
        
        if Deci.Plot.GrandAverage
            SubjectList = {'Group Average'};
        end
        
        
        for cond = 1:length(Deci.Plot.Draw)
            clear flfh_fig fltime_fig fhtime_fig
            
            
            for subj = 1:size(ConnData,1)
                
                
                if Params.FL_FH.do
                    flfh_fig(subj) = figure;
                    
                    if Deci.Plot.Stat.do
                        %dc_pmask(flfh_fig)
                    end
                    suptitle([SubjectList{subj} ' ' Deci.Plot.Title{cond} ' ' conntype{conoi}  ' at time range '  regexprep(num2str(minmax(FL_FH{1,1}.time)),' +',' - ') 's']);
                    
                end
                
                if Params.FL_time.do
                    fltime_fig(subj) = figure;
                    if Deci.Plot.Stat.do
                        %dc_pmask(fltime_fig)
                    end
                    suptitle([SubjectList{subj} ' ' Deci.Plot.Title{cond} ' ' conntype{conoi} ' at Mean Freq High ' num2str(mean(FL_time{1,1}.freqhigh)) 'Hz']);
                end
                
                if Params.FH_time.do
                    fhtime_fig(subj) = figure;
                    if Deci.Plot.Stat.do
                        %dc_pmask(fhtime_fig)
                    end
                    suptitle([SubjectList{subj} ' ' Deci.Plot.Title{cond} ' ' conntype{conoi} ' at Mean Freq Low ' num2str(mean(FH_time{1,1}.freqlow)) 'Hz']);
                end
                
                for subcond = 1:length(Deci.Plot.Draw{cond})
                    
                    if Params.FL_FH.do
                        set(0, 'CurrentFigure', flfh_fig(subj) )
                        flfh_fig(subj).Visible = 'on';
                        subby_flfh(subj,subcond) = subplot(length(Deci.Plot.Draw{cond}),1,subcond );
                        
                        freq = FL_FH{subj,Deci.Plot.Draw{cond}(subcond)};
                        [~,ufreqlow] = unique(freq.freqlow);
                        [~,ufreqhigh] = unique(freq.freqhigh);
                        
                        contourf(subby_flfh(subj,subcond),unique(freq.freqlow),unique(freq.freqhigh),freq.param(ufreqlow,ufreqhigh)')
                        axis tight
                        colorbar;
                        
                        if Deci.Plot.Draw{cond}(subcond) <= size(lockers,2)
                            title([Deci.Plot.Subtitle{cond}{subcond} ' (' num2str(trllen(subj,Deci.Plot.Draw{cond}(subcond))) ')'],'Interpreter','none');
                        else
                            title([Deci.Plot.Subtitle{cond}{subcond}],'Interpreter','none');
                        end
                        
                        ylabel('Freq High')
                        if subcond == length(Deci.Plot.Draw{cond})
                            xlabel('Freq Low')
                        end
                    end
                    
                    if Params.FL_time.do
                        set(0, 'CurrentFigure', fltime_fig(subj) )
                        fltime_fig(subj).Visible = 'on';
                        subby_fltime(subj,subcond) = subplot(length(Deci.Plot.Draw{cond}),1,subcond );
                        
                        
                        freq = FL_time{subj,Deci.Plot.Draw{cond}(subcond)};
                        [~,ufreqlow] = unique(freq.freqlow);
                        
                        contourf(subby_fltime(subj,subcond),freq.time,unique(freq.freqlow),freq.param(ufreqlow,:))
                        colorbar;
                        if Deci.Plot.Draw{cond}(subcond) <= size(lockers,2)
                            xlims = xlim;
                            ylims = ylim;
                            
                            for locks = 1:length([lockers(subj,Deci.Plot.Draw{cond}(subcond),:)])
                                hold on
                                
                                locktime = [lockers(subj,Deci.Plot.Draw{cond}(subcond),locks)/1000];
                                
                                
                                if Deci.Plot.GrandAverage
                                    if locktime > xlims(1) && locktime < xlims(2)
                                        lockstd = [lockersstd(subj,Deci.Plot.Draw{cond}(subcond),locks)/1000];
                                        plotlock = line([locktime locktime], ylims,'LineWidth',2,'Color','k','LineStyle','--','HandleVisibility','off');
                                        
                                        
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
                            title([Deci.Plot.Subtitle{cond}{subcond} ' (' num2str(trllen(subj,Deci.Plot.Draw{cond}(subcond))) ')']);
                        else
                            title([Deci.Plot.Subtitle{cond}{subcond}]);
                        end
                        
                        ylabel('Freq Low')
                        if subcond == length(Deci.Plot.Draw{cond})
                            xlabel('Time')
                        end
                        
                    end
                    
                    if Params.FH_time.do
                        set(0, 'CurrentFigure', fhtime_fig(subj))
                        fhtime_fig(subj).Visible = 'on';
                        subby_fhtime(subj,subcond) = subplot(length(Deci.Plot.Draw{cond}),1,subcond );
                        
                        freq = FH_time{subj,Deci.Plot.Draw{cond}(subcond)};
                        [~,ufreqhigh] = unique(freq.freqhigh);
                        
                        contourf(subby_fhtime(subj,subcond),freq.time,unique(freq.freqhigh),freq.param(ufreqhigh,:))
                        colorbar;
                        if Deci.Plot.Draw{cond}(subcond) <= size(lockers,2)
                            xlims = xlim;
                            ylims = ylim;
                            
                            for locks = 1:length([lockers(subj,Deci.Plot.Draw{cond}(subcond),:)])
                                hold on
                                
                                locktime = [lockers(subj,Deci.Plot.Draw{cond}(subcond),locks)/1000];
                                
                                
                                if Deci.Plot.GrandAverage
                                    if locktime > xlims(1) && locktime < xlims(2)
                                        lockstd = [lockersstd(subj,Deci.Plot.Draw{cond}(subcond),locks)/1000];
                                        plotlock = line([locktime locktime], ylims,'LineWidth',2,'Color','k','LineStyle','--','HandleVisibility','off');
                                        
                                        
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
                            title([Deci.Plot.Subtitle{cond}{subcond} ' (' num2str(trllen(subj,Deci.Plot.Draw{cond}(subcond))) ')']);
                        else
                            title([Deci.Plot.Subtitle{cond}{subcond}]);
                        end
                        
                        ylabel('Freq High')
                        if subcond == length(Deci.Plot.Draw{cond})
                            xlabel('Time')
                        end
                        
                        
                    end
                    
                end
                
            end
            
            %% Normalize Plots
            
            Deci.Plot.Roi = 'maxmin';
            
            for subj = 1:size(ConnData,1)
                
                if Params.FL_FH.do
                    set(0, 'CurrentFigure', flfh_fig(subj) )
                    for r = 1:length(subby_flfh(:))
                        if length(Deci.Plot.Roi) == 2 && isnumeric(Deci.Plot.Roi)
                            subby_flfh(r).CLim = Deci.Plot.Roi;
                        elseif strcmp(Deci.Plot.Roi,'maxmin')
                            if ~isempty(subby_flfh(r).Children.UserData)
                                subby_flfh(r).CLim = [min(arrayfun(@(c) min(c.Children.UserData(:)),subby_flfh(:))) max(arrayfun(@(c) max(c.Children.UserData(:)),subby_flfh(:)))];
                            else
                                subby_flfh(r).CLim = [min(arrayfun(@(c) min(c.Children.ZData(:)),subby_flfh(:))) max(arrayfun(@(c) max(c.Children.ZData(:)),subby_flfh(:)))];
                            end

                        elseif strcmp(Deci.Plot.Roi,'maxabs')
                            if ~isempty(subby_flfh(r).Children.UserData)
                                subby_flfh(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.UserData(:))),subby_flfh(:))) max(arrayfun(@(c) max(abs(c.Children.UserData(:))),subby_flfh(:)))];
                            else
                                subby_flfh(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.ZData(:))),subby_flfh(:))) max(arrayfun(@(c) max(abs(c.Children.ZData(:))),subby_flfh(:)))];
                                
                            end
                            
                        end
                    end
                    
                end
                
                if Params.FL_time.do
                    set(0, 'CurrentFigure', fltime_fig(subj) )
                    for r = 1:length(subby_fltime(:))
                        if length(Deci.Plot.Roi) == 2 && isnumeric(Deci.Plot.Roi)
                            subby_fltime(r).CLim = Deci.Plot.Roi;
                        elseif strcmp(Deci.Plot.Roi,'maxmin')
                            if ~isempty(subby_fltime(r).Children.UserData)
                                subby_fltime(r).CLim = [min(arrayfun(@(c) min(c.Children.UserData(:)),subby_fltime(:))) max(arrayfun(@(c) max(c.Children.UserData(:)),subby_fltime(:)))];
                            else
                                subby_fltime(r).CLim = [min(arrayfun(@(c) min(c.Children.ZData(:)),subby_fltime(:))) max(arrayfun(@(c) max(c.Children.ZData(:)),subby_fltime(:)))];
                            end
                            
                            
                        elseif strcmp(Deci.Plot.Roi,'maxabs')
                            if ~isempty(subby_fltime(r).Children.UserData)
                                subby_fltime(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.UserData(:))),subby_fltime(:))) max(arrayfun(@(c) max(abs(c.Children.UserData(:))),subby_fltime(:)))];
                            else
                                subby_fltime(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.ZData(:))),subby_fltime(:))) max(arrayfun(@(c) max(abs(c.Children.ZData(:))),subby_fltime(:)))];
                                
                            end
                            
                        end
                    end
                    
                end
                
                if Params.FH_time.do
                    set(0, 'CurrentFigure', fhtime_fig(subj) )
                    for r = 1:length(subby_fhtime(:))
                        if length(Deci.Plot.Roi) == 2 && isnumeric(Deci.Plot.Roi)
                            subby_fhtime(r).CLim = Deci.Plot.Roi;
                        elseif strcmp(Deci.Plot.Roi,'maxmin')
                            if ~isempty(subby_fhtime(r).Children.UserData)
                                subby_fhtime(r).CLim = [min(arrayfun(@(c) min(c.Children.UserData(:)),subby_fhtime(:))) max(arrayfun(@(c) max(c.Children.UserData(:)),subby_fhtime(:)))];
                            else
                                subby_fhtime(r).CLim = [min(arrayfun(@(c) min(c.Children.ZData(:)),subby_fhtime(:))) max(arrayfun(@(c) max(c.Children.ZData(:)),subby_fhtime(:)))];
                            end
                            
                        elseif strcmp(Deci.Plot.Roi,'maxabs')
                            if ~isempty(subby_fhtime(r).Children.UserData)
                                subby_fhtime(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.UserData(:))),subby_fhtime(:))) max(arrayfun(@(c) max(abs(c.Children.UserData(:))),subby_fhtime(:)))];
                            else
                                subby_fhtime(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.ZData(:))),subby_fhtime(:))) max(arrayfun(@(c) max(abs(c.Children.ZData(:))),subby_fhtime(:)))];
                                
                            end
                            
                        end
                    end
                    
                end
                
            end
            
            clear subby_fhtime
            clear subby_flfh
            clear subby_fltime
            
        end
    end
    
end