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
        
        
        if isequal(chanl,chanh)
            chancmb = [chanl; chanh]';
        else
            chancmb = [chanl chanh];
            chancmb = chancmb(combvec(1:length(chanl),[1:length(chanh)]+length(chanl)))';
        end
        
        if isequal(freqlow,freqhigh)
            freqcmb = [freqlow; freqhigh]';
        else
            freqcmb = [freqlow freqhigh];
            freqcmb = freqcmb(combvec(1:length(freqlow),[1:length(freqhigh)]+length(freqlow)))';
        end
        
        if isequal(size(freqcmb), [2 1])
            freqcmb = freqcmb';
        end
        
        
        %% load 1 full conoi
        clear sub_cond
        
        for  subject_list = 1:length(Deci.SubjectList)
            
            display(['Loading Plottor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
            for Conditions = 1:length(Deci.Plot.CondTitle)
                
                for choicmb = 1:size(chancmb,1)
                    
                    for foicmb = 1:size(freqcmb,1)
                        
                        connfile = strjoin([chancmb(choicmb,1) chancmb(choicmb,2) freqcmb(foicmb,1) freqcmb(foicmb,2) conntype(conoi)],'_');
                        
                        load([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.CondTitle{Conditions} filesep connfile  '.mat'],'conn');
                        
                        sub_freq{foicmb} = conn;
                        
                        
                    end
                    
                    if ismember(conntype(conoi),{'plv','wpli_debiased','wpli'})
                        param = [conntype{conoi} 'spctrm'];
                        
                        sub_chan{choicmb} = sub_freq{foicmb};
                    else
                        param = 'crsspctrm';
                        
                        ufl = unique(freqcmb(:,1));
                        
                        for fl = 1:length(ufl)
                            fldat = cellfun(@(c) c.crsspctrm,sub_freq(ismember(freqcmb(:,1),ufl{fl})),'UniformOutput',false);
                            tempdata{fl} = cat(2,fldat{:});
                        end
                        clear fldat
                        
                        sub_chan{choicmb} = sub_freq{1};
                        sub_chan{choicmb}.crsspctrm = cat(1,tempdata{:});
                        
                    end
                    clear sub_freq
                end
                
                
                if isfield(conn,'trllength')
                    trllen(subject_list,Conditions) = conn.trllength;
                else
                    trllen(subject_list,Conditions) = nan;
                end
                
                if isfield(conn,'lockers')
                    LockNum = Deci.Analysis.Locks(ismember(Deci.Analysis.LocksTitle,Deci.Plot.Lock));
                    lockers(subject_list,Conditions,:) = conn.lockers - conn.lockers(LockNum);
                else
                    lockers(subject_list,Conditions,:) = nan;
                end
                
                if ismember(conntype(conoi),{'plv','wpli_debiased','wpli'})
                    param = [conntype{conoi} 'spctrm'];
                    
                    uchoi = unique(chancmb(:,1));
                    
                    for choi = 1:length(uchoi)
                        chandat = cellfun(@(c) c.(param),sub_chan(ismember(chancmb(:,1),uchoi{choi})),'UniformOutput',false);
                        tempdata{choi} = cat(3,chandat{:});
                    end
                    clear chandat
                    
                    sub_cond{subject_list,Conditions} = sub_chan{1};
                    sub_cond{subject_list,Conditions}.(param) = permute(cat(4,tempdata{:}),[4 3 1 2]);
                    sub_cond{subject_list,Conditions}.labelcmb = chancmb;
                    sub_cond{subject_list,Conditions}.dimord = ['chanlow_chanhigh_' sub_cond{subject_list,Conditions}.dimord];
                else
                    param = 'crsspctrm';
                    
                    chandat = cellfun(@(c) c.crsspctrm, sub_chan,'un',0);
                    
                    sub_cond{subject_list,Conditions} = sub_chan{1};
                    sub_cond{subject_list,Conditions}.crsspctrm = permute(cat(length(strsplit(sub_chan{1}.dimord,'_'))+1,chandat{:}),[length(strsplit(sub_chan{1}.dimord,'_'))+1 1:length(strsplit(sub_chan{1}.dimord,'_'))]);
                    sub_cond{subject_list,Conditions}.dimord = ['chan_' sub_cond{subject_list,Conditions}.dimord];
                    sub_cond{subject_list,Conditions}.label = chancmb(:,1)';
                end
                clear sub_chan
            end
        end
        
        %% Bsl Correction
        
        
        
        %% Math
        if ~isempty(Deci.Plot.Math)
            
            display(' ')
            display(['Doing ' num2str(length(Deci.Plot.Math)) ' Maths'] )
            
            for conds = 1:length(Deci.Plot.Math)
                for subj = 1:size(sub_cond,1)
                    scfg.parameter = param;
                    scfg.operation = Deci.Plot.Math{conds};
                    
                    arginstr = sprintf('x%i,', 1:length([sub_cond(subj,:)]));
                    arginstr = arginstr(1:end-1); % remove the trailing ','
                    eval(sprintf('operation = @(x) %s;',regexprep( scfg.operation,'x(\d*)','x{$1}')));
                    
                    MathData{subj} =  sub_cond{subj,1};
                    MathData{subj}.(param) = feval(operation, cellfun(@(c) c.(param),sub_cond(subj,:),'un',0));
                end
                
                sub_cond(:,size(sub_cond,2)+1) = MathData;
            end
        end
        
        %% Averaging
        
        if size(sub_cond,1) == 1
            Deci.Plot.GrandAverage = false;
        end
        
        
        for conds = 1:size(sub_cond,2)
            
            if Deci.Plot.GrandAverage
                
                
                if any(~isnan(trllen))
                    trlstd = nanstd(trllen,[],1);
                    trllen = nanmean(trllen,1);
                end
                
                if any(~isnan(lockers))
                    lockersstd = nanstd(lockers,[],1);
                    lockers = nanmean(lockers,1);
                end
                
                
                Subjs = cellfun(@(c) c.(param),sub_cond(:,conds),'UniformOutput',false);
                Subjs = cat(length(strsplit(sub_cond{1}.dimord,'_'))+1,Subjs{:});
                %
                Subjs = permute(Subjs,[length(strsplit(sub_cond{1}.dimord,'_'))+1 1:length(strsplit(sub_cond{1}.dimord,'_'))]);
                StatsData{1,conds} = sub_cond{1,conds};
                StatsData{1,conds}.(param) = Subjs;
                StatsData{1,conds}.dimord = ['subj_' sub_cond{1,conds}.dimord];
                
                Subjs = permute(nanmean(Subjs,1),[2:length(strsplit(sub_cond{1}.dimord,'_'))+1 1]);
                ConnData{1,conds} = sub_cond{1,conds};
                ConnData{1,conds}.(param) = Subjs;
                clear Subjs;
            else
                Deci.Plot.Stat.do = false;
                ConnData(:,conds) = sub_cond(:,conds);
            end
            
   %%         
            
%    if Deci.Plot.Stat.do
%        info.isfreq = 0;
%        info.isconn = 1;
%        StatData = dc_plotstat(Deci,SegStatdata,info);
%    end
   
            %% Data Management
            for subjs = 1:size(ConnData,1)
                
                dim = strsplit(ConnData{subjs,conds}.dimord,'_');
                
                if Deci.Plot.Stat.do
                    dimstat = strsplit(StatsData{subjs,conds}.dimord,'_');
                end
                
                if Params.FL_FH.do && all(ismember({'freqlow' 'freqhigh'},dim))
                    FL_FH{subjs,conds} = ConnData{subjs,conds};
                    FL_FH{subjs,conds}.(param) = permute(mean(FL_FH{subjs,conds}.(param),[find(~ismember(dim,{'freqlow' 'freqhigh'}))],'omitnan'),[find(ismember(dim,{'freqlow' 'freqhigh'})) find(~ismember(dim,{'freqlow' 'freqhigh'}))]);
                    FL_FH{subjs,conds}.dimord = 'freqlow_freqhigh';
                    
                    if Deci.Plot.Stat.do
                        FL_FH_Stat{subjs,conds} = StatsData{subjs,conds};
                        FL_FH_Stat{subjs,conds}.(param) = permute(mean(FL_FH_Stat{subjs,conds}.(param),[find(~ismember(dimstat,{'freqlow' 'freqhigh'}))],'omitnan'),[find(ismember(dimstat,{'freqlow' 'freqhigh'})) find(~ismember(dimstat,{'freqlow' 'freqhigh'}))]);
                        FL_FH_Stat{subjs,conds}.dimord = 'subj_freqlow_freqhigh';
                    end
                    
                    FLFH = true;
                else
                    FLFH = false;
                end
                
                if Params.FL_time.do && [all(ismember({'freq' 'time'},dim)) || all(ismember({'freqlow' 'time'},dim))]
                    FL_time{subjs,conds} = ConnData{subjs,conds};
                    FL_time{subjs,conds}.(param) = permute(mean(FL_time{subjs,conds}.(param),[find(~ismember(dim,{'freq' 'time' 'freqlow'}))],'omitnan'),[find(ismember(dim,{'freq' 'time' 'freqlow'})) find(~ismember(dim,{'freq' 'time' 'freqlow'}))]);
                    FL_time{subjs,conds}.dimord = 'freq_time';
                    
                    if Deci.Plot.Stat.do
                        FL_time_Stat{subjs,conds} = StatsData{subjs,conds};
                        FL_time_Stat{subjs,conds}.(param) = permute(mean(FL_time_Stat{subjs,conds}.(param),[find(~ismember(dimstat,{'freq' 'time' 'freqlow'}))],'omitnan'),[find(ismember(dimstat,{'freq' 'time' 'freqlow'})) find(~ismember(dimstat,{'freq' 'time' 'freqlow'}))]);
                        FL_time_Stat{subjs,conds}.dimord = 'subj_freq_time';
                    end
                    
                    FLtime = true;
                else
                    FLtime = false;
                end
                
                if Params.FH_time.do && ~all(ismember({'freq' 'time'},dim)) && all(ismember({'freqhigh' 'time'},dim))
                    FH_time{subjs,conds} = ConnData{subjs,conds};
                    FH_time{subjs,conds}.param = permute(mean(FH_time{subjs,conds}.(param),[find(~ismember(dim,{'freq' 'time' 'freqhigh'}))],'omitnan'),[find(ismember(dim,{'freq' 'time' 'freqhigh'})) find(~ismember(dim,{'freq' 'time' 'freqhigh'}))]);
                    FH_time{subjs,conds}.dimord = 'freq_time';
                    
                    if Deci.Plot.Stat.do
                        FH_time_Stat{subjs,conds} = StatsData{subjs,conds};
                        FH_time_Stat{subjs,conds}.(param) = permute(mean(FH_time_Stat{subjs,conds}.(param),[find(~ismember(dimstat,{'freq' 'time' 'freqhigh'}))],'omitnan'),[find(ismember(dimstat,{'freq' 'time' 'freqhigh'})) find(~ismember(dimstat,{'freq' 'time' 'freqhigh'}))]);
                        FH_time_Stat{subjs,conds}.dimord = 'subj_freq_time';
                    end
                    
                    FHtime = true;
                else
                    FHtime = false;
                end
                
                if Params.CL_time.do && [all(ismember({'chan' 'time'},dim)) || all(ismember({'chanlow' 'time'},dim))]
                    CL_time{subjs,conds} = ConnData{subjs,conds};
                    CL_time{subjs,conds}.(param) = permute(mean(CL_time{subjs,conds}.(param),[find(~ismember(dim,{'chan' 'time' 'chanlow'}))],'omitnan'),[find(ismember(dim,{'chan' 'time' 'chanlow'})) find(~ismember(dim,{'chan' 'time' 'chanlow'}))]);
                    CL_time{subjs,conds}.dimord = 'chan_time';
                    
                    if Deci.Plot.Stat.do
                        CL_time_Stat{subjs,conds} = StatsData{subjs,conds};
                        CL_time_Stat{subjs,conds}.(param) = permute(mean(CL_time_Stat{subjs,conds}.(param),[find(~ismember(dimstat,{'chan' 'time' 'chanlow'}))],'omitnan'),[find(ismember(dimstat,{'chan' 'time' 'chanlow'})) find(~ismember(dimstat,{'chan' 'time' 'chanlow'}))]);
                        CL_time_Stat{subjs,conds}.dimord = 'subj_chan_time';
                    end
                    
                    CLtime = true;
                else
                    CLtime = false;
                end
                
                if Params.CH_time.do && ~all(ismember({'chan' 'time'},dim)) && all(ismember({'chanhigh' 'time'},dim))
                    CH_time{subjs,conds} = ConnData{subjs,conds};
                    CH_time{subjs,conds}.(param) = permute(mean(CH_time{subjs,conds}.(param),[find(~ismember(dim,{'chan' 'time' 'chanhigh'}))],'omitnan'),[find(ismember(dim,{'chan' 'time' 'chanhigh'})) find(~ismember(dim,{'chan' 'time' 'chanhigh'}))]);
                    CH_time{subjs,conds}.dimord = 'chan_time';
                    
                    if Deci.Plot.Stat.do
                        CH_time_Stat{subjs,conds} = StatsData{subjs,conds};
                        CH_time_Stat{subjs,conds}.(param) = permute(mean(CH_time_Stat{subjs,conds}.(param),[find(~ismember(dimstat,{'chan' 'time' 'chanhigh'}))],'omitnan'),[find(ismember(dimstat,{'chan' 'time' 'chanhigh'})) find(~ismember(dimstat,{'chan' 'time' 'chanhigh'}))]);
                        CH_time_Stat{subjs,conds}.dimord = 'subj_chan_time';
                    end
                    
                    CHtime = true;
                else
                    CHtime = false;
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
                
                
                if FLFH
                    flfh_fig(subj) = figure;
                    
                    if Deci.Plot.Stat.do
                        %dc_pmask(flfh_fig)
                    end
                    suptitle([SubjectList{subj} ' ' Deci.Plot.Title{cond} ' ' conntype{conoi}  ]); %' at time range '  regexprep(num2str(minmax(FL_FH{1,1}.time)),' +',' - ') 's']);
                    
                end
                
                if FLtime
                    fltime_fig(subj) = figure;
                    if Deci.Plot.Stat.do
                        %dc_pmask(fltime_fig)
                    end
                    if ismember(dim, {'freqhigh'})
                        suptitle([SubjectList{subj} ' ' Deci.Plot.Title{cond} ' ' conntype{conoi} ' at Mean Freq High ' num2str(mean(FL_time{1,1}.freqhigh)) 'Hz']);
                    else
                        suptitle([SubjectList{subj} ' ' Deci.Plot.Title{cond} ' ' conntype{conoi}]);
                        
                        
                    end
                end
                
                if FHtime
                    fhtime_fig(subj) = figure;
                    if Deci.Plot.Stat.do
                        %dc_pmask(fhtime_fig)
                    end
                    suptitle([SubjectList{subj} ' ' Deci.Plot.Title{cond} ' ' conntype{conoi} ' at Mean Freq Low ' num2str(mean(FH_time{1,1}.freqlow)) 'Hz']);
                end
                
                if CLtime
                    cltime_fig(subj) = figure;
                    if Deci.Plot.Stat.do
                        %dc_pmask(fltime_fig)
                    end
                    suptitle([SubjectList{subj} ' ' Deci.Plot.Title{cond} ' ' conntype{conoi} ' at Mean Chan highs ']);
                end
                
                if CHtime
                    chtime_fig(subj) = figure;
                    if Deci.Plot.Stat.do
                        %dc_pmask(fhtime_fig)
                    end
                    suptitle([SubjectList{subj} ' ' Deci.Plot.Title{cond} ' ' conntype{conoi} ' at Mean Chan lows ']);
                end
                
                for subcond = 1:length(Deci.Plot.Draw{cond})
                    
                    if FLFH
                        set(0, 'CurrentFigure', flfh_fig(subj) )
                        flfh_fig(subj).Visible = 'on';
                        subby_flfh(subj,subcond) = subplot(length(Deci.Plot.Draw{cond}),1,subcond );
                        colormap('jet')
                        
                        freq = FL_FH{subj,Deci.Plot.Draw{cond}(subcond)};
                        [~,ufreqlow] = unique(freq.freqlow);
                        [~,ufreqhigh] = unique(freq.freqhigh);
                        
                        imagesc(subby_flfh(subj,subcond),unique(freq.freqlow),unique(freq.freqhigh),freq.(param)(ufreqlow,ufreqhigh)')
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
                    
                    if FLtime
                        set(0, 'CurrentFigure', fltime_fig(subj) )
                        fltime_fig(subj).Visible = 'on';
                        subby_fltime(subj,subcond) = subplot(length(Deci.Plot.Draw{cond}),1,subcond );
                        
                        
                        freq = FL_time{subj,Deci.Plot.Draw{cond}(subcond)};
                        [~,ufreqlow] = unique(freq.(dim{ismember(dim,{'freq' 'freqlow'})}));
                        
                        contourf(subby_fltime(subj,subcond),freq.time,unique(freq.(dim{ismember(dim,{'freq' 'freqlow'})})),freq.(param)(ufreqlow,:));
                        colorbar;
                        colormap('jet')
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
                        
                        ylabel('Frequency (Hz)')  %changed 'Freq Low' -> 'Frequency (Hz)'
                        if subcond == length(Deci.Plot.Draw{cond})
                            xlabel('Time')
                        end
                        
                    end
                    
                    if FHtime
                        set(0, 'CurrentFigure', fhtime_fig(subj))
                        fhtime_fig(subj).Visible = 'on';
                        subby_fhtime(subj,subcond) = subplot(length(Deci.Plot.Draw{cond}),1,subcond );
                        
                        freq = FH_time{subj,Deci.Plot.Draw{cond}(subcond)};
                        [~,ufreqhigh] = unique(freq.freqhigh);
                        
                        contourf(subby_fhtime(subj,subcond),freq.time,unique(freq.freqhigh),freq.param(ufreqhigh,:))
                        colorbar;
                        colormap('jet')
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
                    
                    %chans
                    
                    
                    if CLtime
                        set(0, 'CurrentFigure', cltime_fig(subj) )
                        cltime_fig(subj).Visible = 'on';
                        subby_cltime(subj,subcond) = subplot(length(Deci.Plot.Draw{cond}),1,subcond );
                        
                        
                        freq = CL_time{subj,Deci.Plot.Draw{cond}(subcond)};
                        
                        if isfield(freq,'labelcmb')
                        imagesc(subby_cltime(subj,subcond),freq.time,1:length(unique(freq.labelcmb(:,1))),freq.(param));
                                                
                        subby_cltime(subj,subcond).YTickLabel = unique(freq.labelcmb(:,1));
                        subby_cltime(subj,subcond).YTick = 1:length(unique(freq.labelcmb(:,1)));
                        else
                        imagesc(subby_cltime(subj,subcond),freq.time,1:length(freq.label),freq.(param));    
                        subby_cltime(subj,subcond).YTickLabel =freq.label;
                        subby_cltime(subj,subcond).YTick = 1:length(freq.label);
                        end

                        colormap('jet')
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
                        
                        ylabel('Channels Low')
                        if subcond == length(Deci.Plot.Draw{cond})
                            xlabel('Time')
                        end
                        
                    end
                    
                    if CHtime
                        set(0, 'CurrentFigure', chtime_fig(subj) )
                        chtime_fig(subj).Visible = 'on';
                        subby_chtime(subj,subcond) = subplot(length(Deci.Plot.Draw{cond}),1,subcond );
                        
                        
                        freq = CH_time{subj,Deci.Plot.Draw{cond}(subcond)};
                        
                        imagesc(subby_chtime(subj,subcond),freq.time,1:length(unique(freq.labelcmb(:,2))),freq.(param));
                        colormap('jet')
                        subby_chtime(subj,subcond).YTickLabel = unique(freq.labelcmb(:,2));
                        subby_chtime(subj,subcond).YTick = 1:length(unique(freq.labelcmb(:,2)));
                        
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
                        
                        ylabel('Channels High')
                        if subcond == length(Deci.Plot.Draw{cond})
                            xlabel('Time')
                        end
                        
                        
                    end
                    
                    
                    
                end
                
            end
            
            %% Normalize Plots
            
            Deci.Plot.Roi = 'maxmin';
            
            for subj = 1:size(ConnData,1)
                
                if FLFH
                    set(0, 'CurrentFigure', flfh_fig(subj) )
                    for r = 1:length(subby_flfh(:))
                        if length(Deci.Plot.Roi) == 2 && isnumeric(Deci.Plot.Roi)
                            subby_flfh(r).CLim = Deci.Plot.Roi;
                        elseif strcmp(Deci.Plot.Roi,'maxmin')
                            if ~isempty(subby_flfh(r).Children.UserData)
                                subby_flfh(r).CLim = [min(arrayfun(@(c) min(c.Children.UserData(:)),subby_flfh(:))) max(arrayfun(@(c) max(c.Children.UserData(:)),subby_flfh(:)))];
                            else
                                subby_flfh(r).CLim = [min(arrayfun(@(c) min(c.Children.CData(:)),subby_flfh(:))) max(arrayfun(@(c) max(c.Children.CData(:)),subby_flfh(:)))];
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
                
                if FLtime
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
                
                if FHtime
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
                
                if CLtime
                    set(0, 'CurrentFigure', cltime_fig(subj) )
                    for r = 1:length(subby_cltime(:))
                        if length(Deci.Plot.Roi) == 2 && isnumeric(Deci.Plot.Roi)
                            subby_cltime(r).CLim = Deci.Plot.Roi;
                        elseif strcmp(Deci.Plot.Roi,'maxmin')
                            if ~isempty(subby_cltime(r).Children.UserData)
                                subby_cltime(r).CLim = [min(arrayfun(@(c) min(c.Children.UserData(:)),subby_cltime(:))) max(arrayfun(@(c) max(c.Children.UserData(:)),subby_cltime(:)))];
                            else
                                subby_cltime(r).CLim = [min(arrayfun(@(c) min(c.Children.CData(:)),subby_cltime(:))) max(arrayfun(@(c) max(c.Children.CData(:)),subby_cltime(:)))];
                            end
                            
                        elseif strcmp(Deci.Plot.Roi,'maxabs')
                            if ~isempty(subby_cltime(r).Children.UserData)
                                subby_cltime(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.UserData(:))),subby_cltime(:))) max(arrayfun(@(c) max(abs(c.Children.UserData(:))),subby_cltime(:)))];
                            else
                                subby_cltime(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.CData(:))),subby_cltime(:))) max(arrayfun(@(c) max(abs(c.Children.CData(:))),subby_cltime(:)))];
                                
                            end
                            
                        end
                    end
                    
                end
                
                if CHtime
                    set(0, 'CurrentFigure', chtime_fig(subj) )
                    for r = 1:length(subby_chtime(:))
                        if length(Deci.Plot.Roi) == 2 && isnumeric(Deci.Plot.Roi)
                            subby_chtime(r).CLim = Deci.Plot.Roi;
                        elseif strcmp(Deci.Plot.Roi,'maxmin')
                            if ~isempty(subby_chtime(r).Children.UserData)
                                subby_chtime(r).CLim = [min(arrayfun(@(c) min(c.Children.UserData(:)),subby_chtime(:))) max(arrayfun(@(c) max(c.Children.UserData(:)),subby_chtime(:)))];
                            else
                                subby_chtime(r).CLim = [min(arrayfun(@(c) min(c.Children.CData(:)),subby_chtime(:))) max(arrayfun(@(c) max(c.Children.CData(:)),subby_chtime(:)))];
                            end
                            
                        elseif strcmp(Deci.Plot.Roi,'maxabs')
                            if ~isempty(subby_chtime(r).Children.UserData)
                                subby_chtime(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.UserData(:))),subby_chtime(:))) max(arrayfun(@(c) max(abs(c.Children.UserData(:))),subby_chtime(:)))];
                            else
                                subby_chtime(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.CData(:))),subby_chtime(:))) max(arrayfun(@(c) max(abs(c.Children.CData(:))),subby_chtime(:)))];
                                
                            end
                            
                        end
                    end
                    
                end
                
            end
            
            clear subby_fhtime
            clear subby_flfh
            clear subby_fltime
            clear subby_cltime
            clear subby_chtime
        end
    end
    
end