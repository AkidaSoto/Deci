function Plottor_Conn(Deci,Params)

%% Load
cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';

disp('----------------------');
display(' ')

display(['Plotting Conn']);
warning off
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
        
        if isequal(chanl,{'Reinhart-All'})
            chanl = dc_getchans('noeyes');
        end
        
        if isequal(chanh,{'Reinhart-All'})
            chanh = dc_getchans('noeyes');
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
            
            if size(chancmb,2) ~= 1
                chancmb = chancmb(combvec(1:length(chanl),[1:length(chanh)]+length(chanl)))';
            end
        end
        
        if ismember(conntype(conoi),{'plv','wpli_debiased','wpli','amplcorr','powcorr'})
            chancmb = chancmb(cellfun(@(a,b) ~isequal(a,b),chancmb(:,1),chancmb(:,2)),:);
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
        clear sub_cond Foi tempdata
        Deci.Plot.Extra.Conn = Exist(Deci.Plot.Extra.Conn,'toi',[-inf inf]);
        
        
        for  subject_list = 1:length(Deci.SubjectList)
            
            display(['Loading Plottor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
            for Conditions = 1:length(Deci.Plot.CondTitle)
                
                for choicmb = 1:size(chancmb,1)
                    
                    for foicmb = 1:size(freqcmb,1)
                        
                        connfile = strjoin([chancmb(choicmb,1) chancmb(choicmb,2) freqcmb(foicmb,1) freqcmb(foicmb,2) conntype(conoi)],'_');
                        
                        
                        load([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.CondTitle{Conditions} filesep connfile  '.mat'],'conn');
                        
                        sub_freq{foicmb} = conn;
                        
                        
                        if isfield(conn,'time')
                            toi = round(conn.time,4) >= Deci.Plot.Extra.Conn.toi(1) & round(conn.time,4) <= Deci.Plot.Extra.Conn.toi(2);
                        end
                        
                        %                         try
                        %                             Foi(subject_list,:) = conn.freq;
                        %                         catch
                        %                             error(['mismatch in frequencies at subject #' subject_list ', confirm and try reanalyzing?']);
                        %                         end
                        
                    end
                    
                    if ismember(conntype(conoi),{'plv','wpli_debiased','wpli','coh','amplcorr','powcorr'})
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
                
                if ismember(conntype(conoi),{'plv','wpli_debiased','wpli','coh','amplcorr','powcorr'})
                    param = [conntype{conoi} 'spctrm'];
                    
                    uchoil = unique(chancmb(:,1),'stable');
                    uchoih = unique(chancmb(:,2),'stable');
                    
                    tempdata.(param) = nan([length(uchoil) length(uchoih) length(sub_chan{1}.freq) length(sub_chan{1}.time)]);
                    
                    for choi = 1:size(chancmb,1)
                        
                        chanlow = ismember(uchoil,chancmb(choi,1));
                        chanhigh = ismember(uchoih,chancmb(choi,2));
                        
                            tempdata.(param)(chanlow,chanhigh,:,:) = sub_chan{choi}.(param);
                    end
                   
                    sub_cond{subject_list,Conditions} = sub_chan{1};
                    sub_cond{subject_list,Conditions}.(param) = tempdata.(param);
                    sub_cond{subject_list,Conditions}.labelcmb = chancmb;
                    
                    sub_cond{subject_list,Conditions}.chanlow = uchoil;
                    sub_cond{subject_list,Conditions}.chanhigh =  uchoih;
                    
                    
                    
                    sub_cond{subject_list,Conditions}.dimord = ['chanlow_chanhigh_' sub_cond{subject_list,Conditions}.dimord];
                else
                    param = 'crsspctrm';
                    
                    chandat = cellfun(@(c) c.crsspctrm, sub_chan,'un',0);
                    
                    sub_cond{subject_list,Conditions} = sub_chan{1};
                    sub_cond{subject_list,Conditions}.crsspctrm = permute(cat(length(strsplit(sub_chan{1}.dimord,'_'))+1,chandat{:}),[length(strsplit(sub_chan{1}.dimord,'_'))+1 1:length(strsplit(sub_chan{1}.dimord,'_'))]);
                    sub_cond{subject_list,Conditions}.dimord = ['chan_' sub_cond{subject_list,Conditions}.dimord];
                    sub_cond{subject_list,Conditions}.label = chancmb(:,1)';
                end
                
                if isfield(sub_cond(subject_list,Conditions),'dof')
                    sub_cond(subject_list,Conditions) = rmfield(sub_cond(subject_list,Conditions),'dof');
                end
                
            end
            clear sub_chan
        end
        
        
        %% Bsl Correction
        
        
        display(' ');
        display(['Using Lock: ' Deci.Plot.Lock]);
        display(['Using Ref: ' Deci.Plot.BslRef ' at times ' strrep(regexprep(num2str(Deci.Plot.Bsl),' +',' '),' ','-')]);
        
        
        if ~strcmp(Deci.Plot.BslType,'none') && isfield(sub_cond{subject_list,Conditions},'time')
            
            for  subject_list = 1:length(Deci.SubjectList)
                
                display(['Loading BSL for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
                for Conditions = 1:length(Deci.Plot.CondTitle)
                    
                    
                    if ~strcmpi(Deci.Plot.BslRef,Deci.Plot.Lock) || ~isempty(Deci.Plot.LockCond)
                        
                        
                        for choicmb = 1:size(chancmb,1)
                            
                            for foicmb = 1:size(freqcmb,1)
                                
                                if ~isempty(Deci.Plot.LockCond)
                                    BslCond =    Deci.Plot.CondTitle{Deci.Plot.LockCond(Conditions)};
                                else
                                    BslCond =    Deci.Plot.CondTitle{Conditions};
                                end
                                
                                
                                connfile = strjoin([chancmb(choicmb,1) chancmb(choicmb,2) freqcmb(foicmb,1) freqcmb(foicmb,2) conntype(conoi)],'_');
                                
                                load([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep BslCond filesep connfile  '.mat'],'conn');
                                
                                bslfreq{foicmb} = conn;
                                
                                
                            end
                            
                            if ismember(conntype(conoi),{'plv','wpli_debiased','wpli','amplcorr','powcorr'})
                                param = [conntype{conoi} 'spctrm'];
                                
                                bsl_chan{choicmb} = bslfreq{foicmb};
                            else
                                param = 'crsspctrm';
                                
                                ufl = unique(freqcmb(:,1));
                                
                                for fl = 1:length(ufl)
                                    fldat = cellfun(@(c) c.crsspctrm,bslfreq(ismember(freqcmb(:,1),ufl{fl})),'UniformOutput',false);
                                    tempdata{fl} = cat(2,fldat{:});
                                end
                                clear fldat
                                
                                bsl_chan{choicmb} = bslfreq{1};
                                bsl_chan{choicmb}.crsspctrm = cat(1,tempdata{:});
                                
                            end
                            clear bsl_freq
                        end
                        
                        
                        if ismember(conntype(conoi),{'plv','wpli_debiased','wpli','amplcorr','powcorr'})
                            param = [conntype{conoi} 'spctrm'];
                            
                            uchoi = unique(chancmb(:,1));
                            
                            for choi = 1:length(uchoi)
                                chandat = cellfun(@(c) c.(param),bsl_chan(ismember(chancmb(:,1),uchoi{choi})),'UniformOutput',false);
                                tempdata{choi} = cat(3,chandat{:});
                            end
                            clear chandat
                            
                            bsl_cond{subject_list,Conditions} = bsl_chan{1};
                            bsl_cond{subject_list,Conditions}.(param) = permute(cat(4,tempdata{:}),[4 3 1 2]);
                            bsl_cond{subject_list,Conditions}.labelcmb = chancmb;
                            bsl_cond{subject_list,Conditions}.dimord = ['chanlow_chanhigh_' bsl_cond{subject_list,Conditions}.dimord];
                        else
                            param = 'crsspctrm';
                            
                            chandat = cellfun(@(c) c.crsspctrm, bsl_chan,'un',0);
                            
                            bsl_cond{subject_list,Conditions} = bsl_chan{1};
                            bsl_cond{subject_list,Conditions}.crsspctrm = permute(cat(length(strsplit(bsl_chan{1}.dimord,'_'))+1,chandat{:}),[length(strsplit(bsl_chan{1}.dimord,'_'))+1 1:length(strsplit(bsl_chan{1}.dimord,'_'))]);
                            bsl_cond{subject_list,Conditions}.dimord = ['chan_' bsl_cond{subject_list,Conditions}.dimord];
                            bsl_cond{subject_list,Conditions}.label = chancmb(:,1)';
                        end
                        clear sub_chan
                        
                    else
                        bsl_cond{subject_list,Conditions} = sub_cond{subject_list,Conditions};
                    end
                    
                    
                    
                    toi2 = bsl_cond{subject_list,Conditions}.time >= round(Deci.Plot.Bsl(1),4) & bsl_cond{subject_list,Conditions}.time <= round(Deci.Plot.Bsl(2),4);
                    
                    bsl_cond{subject_list,Conditions}.(param) = nanmean(bsl_cond{subject_list,Conditions}.(param)(:,:,:,toi2),4);
                    bsl_cond{subject_list,Conditions}.(param) = repmat(bsl_cond{subject_list,Conditions}.(param),[1 1 1 size(sub_cond{subject_list,Conditions}.(param),4)]);
                    
                    switch Deci.Plot.BslType
                        case 'none'
                        case 'absolute'
                            sub_cond{subject_list,Conditions}.(param) =  sub_cond{subject_list,Conditions}.(param) - bsl_cond{subject_list,Conditions}.(param);
                        case 'relative'
                            sub_cond{subject_list,Conditions}.(param)=  sub_cond{subject_list,Conditions}.(param) ./ bsl_cond{subject_list,Conditions}.(param);
                        case 'relchange'
                            sub_cond{subject_list,Conditions}.(param) = ( sub_cond{subject_list,Conditions}.(param) - bsl_cond{subject_list,Conditions}.(param)) ./ bsl_cond{subject_list,Conditions}.(param);
                        case 'db'
                            sub_cond{subject_list,Conditions}.(param) = 10*log10( sub_cond{subject_list,Conditions}.(param) ./ bsl_cond{subject_list,Conditions}.(param));
                    end
                    
                    if isfield(conn,'time')
                        sub_cond{subject_list,Conditions}.time = sub_cond{subject_list,Conditions}.time(toi);
                        sub_cond{subject_list,Conditions}.(param) = sub_cond{subject_list,Conditions}.(param)(:,:,:,toi);
                    end
                end
                
            end
            
        elseif strcmp(Deci.Plot.BslType,'none') && isfield(sub_cond{subject_list,Conditions},'time')
            
            
            for  subject_list = 1:length(Deci.SubjectList)
                
                for Conditions = 1:length(Deci.Plot.CondTitle)
                    if isfield(conn,'time')
                        sub_cond{subject_list,Conditions}.time = sub_cond{subject_list,Conditions}.time(toi);
                        sub_cond{subject_list,Conditions}.(param) = sub_cond{subject_list,Conditions}.(param)(:,:,:,toi);
                    end
                    
                end
            end
            
        end
    end
    
    
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
    
    
    %% Hemifield
    %         if Deci.Plot.Hemiflip.do
    %             display(' ')
    %             display(['Applying Hemifield Flipping'] )
    %
    %             Deci.Plot.Hemiflip = Exist(Deci.Plot.Hemiflip,'Type','Subtraction');
    %
    %             for conds = 1:size(sub_cond,2)
    %                 for subj = 1:size(sub_cond,1)
    %
    %                     hcfg.parameter = 'avg';
    %                     hcfg.operation = 'x2 - x1';
    %
    %                     ContraCfg.channel = dc_getchans('even');
    %
    %                     ismember(sub_cond{subj,conds}.chanhigh,ContraCfg.channel);
    %
    %                     ContraData{subj,conds} = ft_selectdata(ContraCfg,sub_cond{subj,conds});
    %                     ContraData{subj,conds} = hemifieldflip(ContraData{subj,conds},);
    %
    %                     IpsiCfg.channel = dc_getchans('odd');
    %                     IpsiData{subj,conds} = ft_selectdata(IpsiCfg,sub_cond{subj,conds});
    %
    %                     if strcmpi(Deci.Plot.Hemiflip.Type,'Subtraction')
    %                         sub_cond{subj,conds} = ft_math(hcfg,IpsiData{subj,conds},ContraData{subj,conds});
    %                     end
    %                 end
    %
    %             end
    %
    %             if strcmpi(Deci.Plot.Hemiflip.Type,'Both')
    %                 sub_cond = cat(2,IpsiData,ContraData);
    %
    %                 %Deci.SubjectList = cat(2,cellfun(@(c) [c ' Ipsilateral'],Deci.SubjectList,'un',0),cellfun(@(c) [c ' Contralateral'],Deci.SubjectList,'un',0));
    %
    %                 drawlength =  length(Deci.Plot.Draw) + length(Deci.Plot.Math);
    %
    %                 Deci.Plot.Draw =  cellfun(@(c) [c arrayfun(@(d) d+drawlength,c,'un',1)],Deci.Plot.Draw,'un',0);
    %                 Deci.Plot.Subtitle =  cellfun(@(c) [cellfun(@(d) [d ' Ipsilateral'],c,'un',0) cellfun(@(d) [d ' Contralateral'],c,'un',0)],Deci.Plot.Subtitle,'un',0);
    %
    %             else
    %                 Deci.Plot.Title = cellfun(@(c) [c ' Contra - Ipsilateral'],Deci.Plot.Title,'un',0);
    %             end
    %         end
    %
    %
    %         if Deci.Plot.LeftRight
    %             for conds = 1:size(sub_cond,2)
    %                 for subj = 1:size(sub_cond,1)
    %
    %                     Right = ismember(sub_cond{subj,conds}.chanhigh,dc_getchans('even'));
    %                     Left = ismember(sub_cond{subj,conds}.chanhigh,dc_getchans('odd'));
    %
    %                     Rights = sub_cond{subj,conds};
    %                     Rights.(param) = Rights.(param)(:,Right,)
    %
    %
    %                     sub_cond = cat(2,IpsiData,ContraData);
    %
    %                     %Deci.SubjectList = cat(2,cellfun(@(c) [c ' Ipsilateral'],Deci.SubjectList,'un',0),cellfun(@(c) [c ' Contralateral'],Deci.SubjectList,'un',0));
    %
    %                     drawlength =  length(Deci.Plot.Draw) + length(Deci.Plot.Math);
    %
    %                     Deci.Plot.Draw =  cellfun(@(c) [c arrayfun(@(d) d+drawlength,c,'un',1)],Deci.Plot.Draw,'un',0);
    %                     Deci.Plot.Subtitle =  cellfun(@(c) [cellfun(@(d) [d ' Ipsilateral'],c,'un',0) cellfun(@(d) [d ' Contralateral'],c,'un',0)],Deci.Plot.Subtitle,'un',0);
    %
    %                     Deci.Plot.Title = cellfun(@(c) [c ' Contra - Ipsilateral'],Deci.Plot.Title,'un',0);
    %
    %                 end
    %             end
    %         end
    
    %%
    
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
        
        
        %% Data Management
        info.parameter = param;
        
        
        AllPlots = {Deci.Plot.Extra.Conn.FL_FH.do Deci.Plot.Extra.Conn.CL_CH.do Deci.Plot.Extra.Conn.FL_time.do ...
            Deci.Plot.Extra.Conn.FH_time.do Deci.Plot.Extra.Conn.CL_time.do Deci.Plot.Extra.Conn.CH_time.do ...
            Deci.Plot.Extra.Conn.CL.do Deci.Plot.Extra.Conn.CH.do};
        
        
        AllPlots_Dims  ={{'freqlow' 'freqhigh'}, {'chanlow' 'chanhigh'}, {'freq' 'time' 'freqlow'}, {'freq' 'time' 'freqhigh'} ...
            {'chan' 'time' 'chanlow'},{'chan' 'time' 'chanhigh'}, {'chan','chanlow'}, {'chan','chanhigh'}} ;
        
        
        for subjs = 1:size(ConnData,1)
            
            dim = strsplit(ConnData{subjs,conds}.dimord,'_');
            
            if Deci.Plot.Stat.do
                dimstat = strsplit(StatsData{subjs,conds}.dimord,'_');
            end
            
            for DataType = 1:length(AllPlots)
                
                if AllPlots{DataType} && sum(ismember(AllPlots_Dims{DataType},dim)) >= [length(AllPlots_Dims{DataType})-1]
                    
                    AllData{DataType}{subjs,conds} = ConnData{subjs,conds};
                    
                    AllData{DataType}{subjs,conds}.(param) = permute(mean(AllData{DataType}{subjs,conds}.(param),[find(~ismember(dim,AllPlots_Dims{DataType}))],'omitnan'),[find(ismember(dim,AllPlots_Dims{DataType})) find(~ismember(dim,AllPlots_Dims{DataType}))]);
                    AllData{DataType}{subjs,conds}.dimord = strjoin(dim(ismember(dim,AllPlots_Dims{DataType})),'_');
                    
                    
                    if Deci.Plot.Stat.do
                        AllData_Stat{DataType}{subjs,conds} = rmfield(StatsData{subjs,conds},[{'labelcmb'} dim(~ismember(dim,AllPlots_Dims{DataType}))]);
                        AllData_Stat{DataType}{subjs,conds}.(param) = permute(mean(AllData_Stat{DataType}{subjs,conds}.(param),[find(~ismember(dimstat,[{'subj'} AllPlots_Dims{DataType}]))],'omitnan'),[find(ismember(dimstat,[{'subj'} AllPlots_Dims{DataType}])) find(~ismember(dimstat,[{'subj'} AllPlots_Dims{DataType}]))]);
                        AllData_Stat{DataType}{subjs,conds}.dimord = strjoin([{'subj'} dim(ismember(dim,AllPlots_Dims{DataType}))],'_');
                    end
                    
                    AllPlots_do{DataType} = true;
                else
                    AllPlots_do{DataType} = false;
                    
                    
                end
            end
            
        end
        
    end
    
    
    %% Stats
    Plot.Draw = Deci.Plot.Draw;
    Plot.Title = Deci.Plot.Title;
    Plot.Subtitle = Deci.Plot.Subtitle;
    
    if Deci.Plot.Stat.do
        info.isfreq = 0;
        info.isconn = 1;
        
        Deci.Plot.Stat.parameter = param;
        
        
        for DataType = 1:length(AllPlots)
            
            if AllPlots_do{DataType}
                AllData_StatData{DataType} = dc_plotstat(Deci,AllData_Stat{DataType},info);
                
                if Deci.Plot.Stat.twoway.do
                    for cond = 1:length(Deci.Plot.Draw)
                        
                        if length(size(AllData_StatData{DataType}{cond}.prob)) == 2
                            
                        end
                        
                        if length(Deci.Plot.Draw{cond}) == 4
                            
                            anovadata = [AllData{DataType}{1,Deci.Plot.Draw{cond}}];
                            AllData{DataType}{1,end+1} = AllData{DataType}{1};
                            AllData{DataType}{1,end}.(param) = mean(cat(5,anovadata([3 4]).(param)),5) - mean(cat(5,anovadata([1 2]).(param)),5);
                            
                            AllData_StatData{DataType}{end+1} = structfun(@(c) c(:,:,:,1),AllData_StatData{DataType}{cond},'UniformOutput',false);
                            
                            AllData{DataType}{1,end+1} = AllData{DataType}{1};
                            AllData{DataType}{1,end}.(param) = mean(cat(5,anovadata([1 3]).(param)),5) - mean(cat(5,anovadata([2 4]).(param)),5);
                            
                            
                            AllData_StatData{DataType}{end+1} = structfun(@(c) c(:,:,:,2),AllData_StatData{DataType}{cond},'UniformOutput',false);
                            
                            AllData{DataType}{1,end+1} = AllData{DataType}{1};
                            AllData{DataType}{1,end}.(param) =   [[anovadata(3).(param)] - [anovadata(4).(param)]]  -  [[anovadata(1).(param)] - [anovadata(2).(param)]];
                            
                            AllData_StatData{DataType}{end+1} = structfun(@(c) c(:,:,:,3),AllData_StatData{DataType}{cond},'UniformOutput',false);
                            
                        end
                    end
                    
                    
                end
                
            end
        end
        
        if Deci.Plot.Stat.twoway.do

            
            Plot.Draw{end+1} = length(AllData{find([AllPlots_do{:}],1,'first')}) -2;
            Plot.Draw{end+1} = length(AllData{find([AllPlots_do{:}],1,'first')}) -1;
            Plot.Draw{end+1} = length(AllData{find([AllPlots_do{:}],1,'first')});
            
            Plot.Title(end+1:end+3)        = Deci.Plot.Stat.twoway.Title;
            Plot.Subtitle(end+1:end+3)   = Deci.Plot.Stat.twoway.Subtitle;
        end
    end
    
    
    
    %% Plot
    
    SubjectList = Deci.SubjectList;
    
    if Deci.Plot.GrandAverage
        SubjectList = {'Group Average'};
    end
    
    
    for cond = 1:length(Plot.Draw)
        clear flfh_fig fltime_fig fhtime_fig cl_fig
        
        
        for subj = 1:size(SubjectList,1)
            
            
            for DataType = 1:length(AllPlots)
                
                if AllPlots_do{DataType}
                    
                    AllPlots_fig{DataType}{cond}(subj) = figure;
                    
                    if Deci.Plot.Stat.do
                        dc_pmask(AllPlots_fig{DataType}{cond}(subj))
                    end
                    
                   
                    
                end
                
            end
            
            for subcond = 1:length(Plot.Draw{cond})
                
                for DataType = 1:length(AllPlots)
                    if AllPlots_do{DataType}
                        set(0, 'CurrentFigure',  AllPlots_fig{DataType}{cond}(subj) )
                        AllPlots_fig{DataType}{cond}(subj).Visible = 'on';
                        
                        AllPlots_subby{DataType}{cond}(subj,subcond) = subplot(length(Plot.Draw{cond}),1,subcond );
                        %colormap(Deci.Plot.ColorMap)
                        
                        pcfg = cfg;
                        pcfg.parameter = param;
                        
                        connplot = [];
                        connplot.(param) =  permute(AllData{DataType}{subj,Plot.Draw{cond}(subcond)}.(param),[3 1 2 4]);
                        if Deci.Plot.Stat.do
                            pcfg.clim = 'maxmin';
                            pcfg.maskparameter ='mask';
                            connplot.mask =  permute(AllData_StatData{DataType}{cond}.mask,[3 1 2 4]); %repmat(,[length(Segdata{subj,Plot.Draw{cond}(subcond)}.label) 1 1]);
                            connplot.mask = double(connplot.mask);
                            connplot.mask(connplot.mask == 0) = .2;
                        end

                        nonsingleton = sum(size(AllData{DataType}{subj,Plot.Draw{cond}(subcond)}.(param)) ~= 1);
                        
                        if nonsingleton == 2
                            pcfg.imagetype = 'imagesc';
                            pcfg.colormap = Deci.Plot.ColorMap;
                            
                            connplot.dimord = 'chan_freq_time';
                            connplot.label = {'dummy'};
                            
                            curdim = dim(ismember(dim,AllPlots_Dims{DataType}));
                            
                            connplot.freq = 1:length(AllData{DataType}{subj,Plot.Draw{cond}(subcond)}.(curdim{1}));
                            
                            if ~ismember({'time'},curdim)
                                connplot.time = 1:length(AllData{DataType}{subj,Plot.Draw{cond}(subcond)}.(curdim{2}));
                                
                            else
                                connplot.time = AllData{DataType}{subj,Plot.Draw{cond}(subcond)}.time;
                            end
                            
                            evalc('ft_singleplotTFR(pcfg,connplot)');
                            
                            axis tight
                            colorbar;
                            
                            if Plot.Draw{cond}(subcond) <= size(lockers,2)
                                title([Plot.Subtitle{cond}{subcond} ' (' num2str(trllen(subj,Plot.Draw{cond}(subcond))) ')'],'Interpreter','none');
                            else
                                title([Plot.Subtitle{cond}{subcond}],'Interpreter','none');
                            end
                            
                            AllPlots_subby{DataType}{cond}(subj,subcond).YTick = 1:length(AllData{DataType}{subj,Plot.Draw{cond}(subcond)}.(curdim{1}));
                            AllPlots_subby{DataType}{cond}(subj,subcond).YTickLabel = AllData{DataType}{subj,Plot.Draw{cond}(subcond)}.(curdim{1});
                            
                            
                            if ~ismember({'time'},curdim)
                                AllPlots_subby{DataType}{cond}(subj,subcond).XTick = 1:length(AllData{DataType}{subj,Plot.Draw{cond}(subcond)}.(curdim{2}));
                                AllPlots_subby{DataType}{cond}(subj,subcond).XTickLabel = AllData{DataType}{subj,Plot.Draw{cond}(subcond)}.(curdim{2});
                                
                            else
                                
                                if Plot.Draw{cond}(subcond) <= size(lockers,2)
                                    xlims = xlim;
                                    ylims = ylim;
                                    
                                    for locks = 1:length([lockers(subj,Plot.Draw{cond}(subcond),:)])
                                        hold on
                                        
                                        locktime = [lockers(subj,Plot.Draw{cond}(subcond),locks)/1000];
                                        
                                        
                                        if Deci.Plot.GrandAverage
                                            if locktime > xlims(1) && locktime < xlims(2)
                                                lockstd = [lockersstd(subj,Plot.Draw{cond}(subcond),locks)/1000];
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
                                    title([Plot.Subtitle{cond}{subcond} ' (' num2str(trllen(subj,Plot.Draw{cond}(subcond))) ')']);
                                else
                                    title([Plot.Subtitle{cond}{subcond}]);
                                end
                                
                                ylabel('Frequency (Hz)')  %changed 'Freq Low' -> 'Frequency (Hz)'
                                if subcond == length(Plot.Draw{cond})
                                    xlabel('Time')
                                end
                                
                            end
                            
                            
                        else
                            
                            connplot.(param) = permute(connplot.(param),[2 1]);
                            
                            if Deci.Plot.Stat.do
                            connplot.mask = permute(connplot.mask,[2 1 3 4]);
                            end
                            connplot.dimord = 'chan';
                            
                            labelcheck = ismember(AllPlots_Dims{DataType},fields(sub_cond{1}));
                            if any(labelcheck)
                              connplot.label = sub_cond{1}.(AllPlots_Dims{DataType}{labelcheck});
                              
                              switch AllPlots_Dims{DataType}{labelcheck}
                                  case 'chanhigh'
                                      tcfg.channel = uchoih(~ismember(uchoih,uchoil));
                                  case 'chanlow'
                                      tcfg.channel = uchoil(~ismember(uchoil,uchoih));
                              end
                              
                              connplot = ft_selectdata(tcfg,connplot);
                              
                            else
                              connplot.label = sub_cond{1}.label;
                            end

                            pcfg.imagetype = Deci.Plot.ImageType;
                            pcfg.comment = 'no';
                            pcfg.style = 'fill';
                            pcfg.markersymbol = '.';
                            pcfg.colormap = Deci.Plot.ColorMap;
                            pcfg.colorbar = 'yes';
                            pcfg.contournum = 15;
                            
                            ft_topoplotER(pcfg, connplot);
                            
                            if Plot.Draw{cond}(subcond) <= size(lockers,2)
                                title([Plot.Subtitle{cond}{subcond} ' (' num2str(trllen(subj,Plot.Draw{cond}(subcond))) ')'],'Interpreter','none');
                            else
                                title([Plot.Subtitle{cond}{subcond}],'Interpreter','none');
                            end
                        end
                    end
                end
                
            end
        end
        
        %% Normalize Plots
        
        for subj = 1:size(SubjectList,1)
            
            for DataType = 1:length(AllPlots)
                
                if AllPlots_do{DataType}
                    
                    set(0, 'CurrentFigure', AllPlots_fig{DataType}{cond}(subj) )
                    childs = AllPlots_fig{DataType}{cond}(subj).Children.findobj('Type','Axes');
                    
                    for r = 1:length(childs)
                        if length(Deci.Plot.Roi) == 2 && isnumeric(Deci.Plot.Roi)
                            AllPlots_subby{DataType}{cond}(r).CLim = Deci.Plot.Roi;
                        elseif strcmp(Deci.Plot.Roi,'maxmin')
                            if isempty(AllPlots_subby{DataType}{cond}(r).Children.findobj('Type','Image'))
                                AllPlots_subby{DataType}{cond}(r).CLim = [min(arrayfun(@(c) min(c.Children.findobj('Type','Image').UserData(:)),AllPlots_subby{DataType}{cond}(:))) max(arrayfun(@(c) max(c.Children.findobj('Type','Image').UserData(:)),AllPlots_subby{DataType}{cond}(:)))];
                            else
                                AllPlots_subby{DataType}{cond}(r).CLim = [min(arrayfun(@(c) min(c.Children.findobj('Type','Image').CData(:)),AllPlots_subby{DataType}{cond}(:))) max(arrayfun(@(c) max(c.Children.findobj('Type','Image').CData(:)),AllPlots_subby{DataType}{cond}(:)))];
                            end
                            
                        elseif strcmp(Deci.Plot.Roi,'maxabs')
                            if ~isempty(AllPlots_subby{DataType}{cond}(1).Children.findobj('Type','Contour'))
                                
                                childs(r).CLim = [-1*max(arrayfun(@(c) nanmax(abs(c.CLim(:))),childs(:))) nanmax(arrayfun(@(c) max(abs(c.CLim(:))),childs(:)))];
                            else
                                childs(r).CLim = [-1*max(arrayfun(@(c) nanmax(abs(c.Children.CData(:))),childs(:))) max(arrayfun(@(c) max(abs(c.Children.CData(:))),AllPlots_subby{DataType}{cond}(:)))];
                                
                            end
                            
                        end
                    end
                     
                    suptitle([SubjectList{subj} ' ' Plot.Title{cond} ' ' conntype{conoi}  ]);
                end
            end
            
        end
        
    end
end
end
