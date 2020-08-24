
function Plottor_Behv(Deci)

for subject_list = 1:length(Deci.SubjectList)
    
    data = [];
    
    switch Deci.Plot.Behv.Source
        case 'PostArt'
            
            info = [];
            load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'info');
            
            if isempty(info)
                load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data');
                data = rmfield(data,'trial');
                
                info = data;
                
                save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'info');
            else
                data = info;
            end
            
            if isfield(data,'condinfo')  %replacer starting 12/22, lets keep for ~4 months
                data.postart.locks = data.condinfo{1};
                data.postart.events = data.condinfo{2};
                data.postart.trlnum = data.condinfo{3};
                
                data.locks = data.preart{1};
                data.events = data.preart{2};
                data.trlnum = data.preart{3};
                
                data = rmfield(data,'condinfo');
                data = rmfield(data,'preart');
            end
            
            
            data.trialinfo = data.postart.locks;
            data.event = data.postart.events;
            data.trialnum = data.postart.trlnum;
            
        case 'Definition'
            load([Deci.Folder.Version filesep 'Definition' filesep Deci.SubjectList{subject_list}]);
            data = cfg;
            data.trialinfo = data.trl(:,end-length(Deci.DT.Locks)+1:end);
    end
    
    
    for fig = find(Deci.Plot.Behv.RT.Figure)
        
        if exist('RT') == 0 || length(RT) < fig
            RT{fig} = [];
        end
        
        
        if ~isempty(Deci.Plot.Behv.RT)
            
            if length(Deci.DT.Locks) <2 || length(Deci.Plot.Behv.RT.Locks) ~= 2
                error('DT.Locks seems to not have enough to calculate Locks(2)-Locks(1)')
            end
            
            
            RTevents = data.event;
            
            if isa(Deci.Plot.Behv.RT.Block,'function_handle')
                RTevents(:,find(data.event(1,:) < 1)) = Deci.Plot.Behv.RT.Block(RTevents(:,find(data.event(1,:) < 1)));
                
                RTBlock= unique(RTevents(:,find(RTevents(1,:) < 1)),'stable');
            elseif isempty(Deci.Plot.Behv.RT.Block)
                RTBlock = {-1};
            end
            
            
            for blk = 1:length(RTBlock)
                for draw = 1:length(Deci.Plot.Behv.RT.Draw{fig})
                    
                    
                    draws = Deci.Analysis.Conditions(Deci.Plot.Behv.RT.Draw{fig}{draw});
                    
                    if isfield(Deci.Plot.Behv,'Static')
                        
                        if sum(ismember(Deci.Plot.Behv.Static,[draws{:}])) == 1
                            
                            eve = RTevents(logical(sum(ismember(RTevents,Deci.Plot.Behv.Static(ismember(Deci.Plot.Behv.Static,[draws{:}]))),2)),:);
                            rteve = data.trialinfo(logical(sum(ismember(RTevents,Deci.Plot.Behv.Static(ismember(Deci.Plot.Behv.Static,[draws{:}]))),2)),:);
                        else
                            eve = RTevents;
                            rteve = data.trialinfo;
                        end
                    else
                        eve = RTevents;
                        rteve = data.trialinfo;
                    end
                    rteve = rteve(find(any(ismember(eve,RTBlock(blk)),2)),:);
                    eve = eve(find(any(ismember(eve,RTBlock(blk)),2)),:);
                    
                    
                    eveTotal = nan([1 length(find(any(ismember(eve,RTBlock(blk)),2)))]);
                    
                    
                    maxt = length(find(cellfun(@(c) any(ismember([draws{:}],c)), Deci.DT.Markers)));
                    trl = [sum(ismember(eve,[draws{:}]),2) == maxt];
                    
                    eveTotal(trl) =  [-rteve(trl,Deci.Plot.Behv.RT.Locks(2)) - -rteve(trl,Deci.Plot.Behv.RT.Locks(1))];
                    
                    if size(eveTotal,2) ~= size(RT{fig},4) && size(RT{fig},1) ~= 0
                        
                        if size(eveTotal,2) > size(RT{fig},4)
                            RT{fig} = cat(4, RT{fig}, nan(size(RT{fig},1),size(RT{fig},2),size(RT{fig},3),size(eveTotal,2)-size(RT{fig},4)));
                            RT{fig}(subject_list,draw,blk,:) = eveTotal;
                        else
                            eveTotal = [eveTotal nan([1 size(RT{fig},4)-size(eveTotal,2)])];
                            RT{fig}(subject_list,draw,blk,:) = eveTotal;
                        end
                        
                    else
                        RT{fig}(subject_list,draw,blk,:) = eveTotal;
                    end
                    
                end
            end
            
            %disp([Deci.SubjectList{subject_list} ' RT:']);
            %disp(RT(subject_list,:,:));
        end
    end
    
    for fig = find(Deci.Plot.Behv.Acc.Figure)
        
        
        if exist('Acc') == 0 || length(Acc) < fig
            Acc{fig} = [];
        end
        
        if ~isempty(Deci.Plot.Behv.Acc)
            Exist(Deci.Plot.Behv.Acc,'Total');
            Exist(Deci.Plot.Behv.Acc,'Subtotal');
            
            
            Accevents = data.event;
            
            
            if isa(Deci.Plot.Behv.Acc.Block,'function_handle')
                Accevents(:,find(data.event(1,:) < 1)) = Deci.Plot.Behv.Acc.Block(data.event(:,find(data.event(1,:) < 1)));
                
                AccBlock= unique(Accevents(:,find(Accevents(1,:) < 1)),'stable');
            elseif isempty(Deci.Plot.Behv.Acc.Block)
                AccBlock = {-1};
            end
            
            for blk = 1:length(AccBlock)
                for draw = 1:length(Deci.Plot.Behv.Acc.Total{fig})
                    
                    draws = Deci.Analysis.Conditions(Deci.Plot.Behv.Acc.Total{fig}{draw});
                    
                    if isfield(Deci.Plot.Behv,'Static')
                        
                        if sum(ismember(Deci.Plot.Behv.Static,[draws{:}])) == 1
                            
                            eve = Accevents(logical(sum(ismember(Accevents,Deci.Plot.Behv.Static(ismember(Deci.Plot.Behv.Static,[draws{:}]))),2)),:);
                        else
                            eve = Accevents;
                        end
                    else
                        eve = Accevents;
                    end
                    eve = eve(find(any(ismember(eve,AccBlock(blk)),2)),:);
                    
                    eveTotal = nan([1 length(find(any(ismember(eve,AccBlock(blk)),2)))]);
                    
                    maxt = length(find(cellfun(@(c) any(ismember([draws{:}],c)), Deci.DT.Markers)));
                    trl = [[sum(ismember(eve,[draws{:}]),2)] == maxt];
                    
                    eveTotal(trl) = 0;
                    
                    
                    subdraws = Deci.Analysis.Conditions(Deci.Plot.Behv.Acc.Subtotal{fig}{draw});
                    subtrl = [[sum(ismember(eve,[subdraws{:}]),2) ] == maxt];
                    eveTotal([subtrl & trl]) = 1;
                    
                    if size(eveTotal,2) ~= size(Acc{fig},4) && size(Acc{fig},1) ~= 0
                        if size(eveTotal,2) > size(Acc{fig},4)
                            Acc{fig} = cat(4, Acc{fig}, nan(size(Acc{fig},1),size(Acc{fig},2),size(Acc{fig},3),size(eveTotal,2)-size(Acc{fig},4)));
                            Acc{fig}(subject_list,draw,blk,:) = eveTotal;
                        else
                            eveTotal = [eveTotal nan([1 size(Acc{fig},4)-size(eveTotal,2)])];
                            Acc{fig}(subject_list,draw,blk,:) = eveTotal;
                        end
                        
                    else
                        Acc{fig}(subject_list,draw,blk,:) = eveTotal;
                    end
                    
                    Deci.Plot.Behv.Acc.Collapse = Exist(Deci.Plot.Behv.Acc.Collapse,'Movmean',[]);
                    
                    if ~isempty(Deci.Plot.Behv.Acc.Collapse.Movmean)
                        
                        if Deci.Plot.Behv.Acc.Collapse.Movmean(fig)
                            
                            if isfield(Deci.Plot.Behv.Acc.Collapse,'MovWindow')
                                
                                MovWindow = Deci.Plot.Behv.Acc.Collapse.MovWindow;
                                
                                if MovWindow == -1
                                    MovWindow = length( Acc{fig}(subject_list,draw,blk,:));
                                end
                                
                            else
                                MovWindow = length( Acc{fig}(subject_list,draw,blk,:));
                            end
                            
                            Acc{fig}(subject_list,draw,blk,:) = movmean(Acc{fig}(subject_list,draw,blk,:),[MovWindow 0],'omitnan');
                        end
                    end
                end
            end
            
            %disp([Deci.SubjectList{subject_list} ' Acc:']);
            %disp(Acc(subject_list,:,:));
        end
        
    end
end

if ~isempty(Deci.Plot.Behv.Acc) && ~isempty(find(Deci.Plot.Behv.Acc.Figure))
    mkdir([Deci.Folder.Version filesep 'Plot'])
    save([Deci.Folder.Version filesep 'Plot' filesep 'Acc'],'Acc');
end

if ~isempty(Deci.Plot.Behv.RT) &&  ~isempty(find(Deci.Plot.Behv.RT.Figure))
    mkdir([Deci.Folder.Version filesep 'Plot'])
    save([Deci.Folder.Version filesep 'Plot' filesep 'RT'],'RT');
end

%% sort

Deci.Plot.Behv.Acc = Exist(Deci.Plot.Behv.Acc,'BaW',false);
for fig = find(Deci.Plot.Behv.Acc.Figure)
    
    clear fAcc
    
    if ~isempty(Deci.Plot.Behv.Acc) && ~isempty(find(Deci.Plot.Behv.Acc.Figure))
        fullAcc{fig} = Acc{fig};
        
        if Deci.Plot.Behv.Acc.Collapse.Trial
            Accsem{fig} = nanstd(Acc{fig},[],4)/sqrt(size(Acc{fig},4));
            Acc{fig} = nanmean(Acc{fig},4);
            
        end
        
        if Deci.Plot.Behv.Acc.Collapse.Block
            Accsem{fig} = nanstd(Acc{fig},[],3)/sqrt(size(Acc{fig},3));
            Acc{fig} = nanmean(Acc{fig},3);
            
        end
        
        Sub.Acc = Deci.SubjectList;
        if Deci.Plot.Behv.Acc.Collapse.Subject

            
            if ~Deci.Plot.Behv.Acc.BaW
            Accsem{fig} =  nanstd(Acc{fig},[],1)/sqrt(size(Acc{fig},1));
            Acc{fig} =  nanmean(Acc{fig},1);
            Sub.Acc = {'SubjAvg'};
            else

%                 if Deci.Plot.Behv.Acc.Collapse.Block && Deci.Plot.Behv.Acc.Collapse.Trial
%                     error('Either Block or Trial has to be not collapsed to do Box and Whiskers')
%                 end
                
            end
        end
        
        %save([Deci.Folder.Version filesep 'Plot' filesep 'SimAcc'],'Acc','fullAcc');
    end
end

Deci.Plot.Behv.RT = Exist(Deci.Plot.Behv.RT,'BaW',false);
for fig = find(Deci.Plot.Behv.RT.Figure)
    
    
    clear fRT
    if ~isempty(Deci.Plot.Behv.RT) &&  ~isempty(find(Deci.Plot.Behv.RT.Figure))
        
        fullRT{fig} = RT{fig};
        
        if Deci.Plot.Behv.RT.Collapse.Trial
            RTsem{fig} = nanstd(RT{fig},[],4)/sqrt(size(RT{fig},4));
            RT{fig} = nanmean(RT{fig},4);
        end
        
        if Deci.Plot.Behv.RT.Collapse.Block
            RTsem{fig} = nanstd(RT{fig},[],3)/sqrt(size(RT{fig},3));
            RT{fig} = nanmean(RT{fig},3);
        end
        
        Sub.RT = Deci.SubjectList;
        if Deci.Plot.Behv.RT.Collapse.Subject
            
            
            
            if ~Deci.Plot.Behv.RT.BaW
                
                RTsem{fig} = nanstd(RT{fig},[],1)/sqrt(size(RT{fig},1));
                RT{fig} =  nanmean(RT{fig},1);
                Sub.RT = {'SubjAvg'};
            else
                
%                 if Deci.Plot.Behv.RT.Collapse.Block && Deci.Plot.Behv.RT.Collapse.Trial
%                     error('Either Block or Trial has to be not collapsed to do Box and Whiskers')
%                 end
                
            end
            
        end
        
        %save([Deci.Folder.Version filesep 'Plot' filesep 'SimRT'],'RT','fullRT');
    end
    
end

   


%% Table outputs for SPSS


Deci.Plot.Behv = Exist(Deci.Plot.Behv,'WriteExcel',false);

if Deci.Plot.Behv.WriteExcel
    
    for fig = find(Deci.Plot.Behv.Acc.Figure)
        
        
        %ExportExcel
        colnames = [];
        for blk = 1:size(fullAcc{fig},3)
            for cond = 1:size(fullAcc{fig},2)
                
                colnames{end+1} = [Deci.Plot.Behv.Acc.Subtitle{fig}{cond} ' : Block ' num2str(abs(AccBlock(blk)))];
            end
        end
        
        exceldata = reshape(nanmean(fullAcc{fig},4),[size(fullAcc{fig},1) prod([size(fullAcc{fig},2) size(fullAcc{fig},3)])]);
        
        exceldata = horzcat(Deci.SubjectList',num2cell(exceldata));
        exceldata = vertcat([{'Accuracy'} colnames],exceldata);
        
        if exist([Deci.Folder.Plot filesep  Deci.Plot.Behv.Acc.Title{fig} ' Behavioral Outputs' ]) == 2
            writematrix([],[Deci.Folder.Plot filesep  Deci.Plot.Behv.Acc.Title{fig} ' Behavioral Outputs' ],'FileType','spreadsheet','Sheet','TempSheet');
            %xls_delete_sheets([Deci.Folder.Plot filesep  Deci.Plot.Behv.Acc.Title{fig} ' Behavioral Outputs' ],'Accuracy_Full');
        end
        
        writecell(exceldata,[Deci.Folder.Plot filesep  Deci.Plot.Behv.Acc.Title{fig} ' Behavioral Outputs' ],'FileType','spreadsheet','Sheet','Accuracy_Full');
        %xls_delete_sheets([Deci.Folder.Plot filesep  Deci.Plot.Behv.Acc.Title{fig} ' Behavioral Outputs' ],'TempSheet');
        
        subs= [];
        conds = [];
        blks = [];
        trls = [];
        
        
        for trl = 1:size(Acc{fig},4)
            for blk = 1:size(Acc{fig},3)
                for cond = 1:size(Acc{fig},2)
                    for sub = 1:size(Acc{fig},1)
                        subs(end+1) =  sub;
                        conds(end+1) = cond;
                        blks(end+1) = blk;
                        trls(end+1) = trl;
                    end
                end
            end
        end
        
        
        sAcc = size(Acc{fig});
        
        fAcc = reshape(Acc{fig},[prod([sAcc]) 1]);
        fAccsem = reshape(Accsem{fig},[prod([sAcc]) 1]);
        
        
        
        if exist([Deci.Folder.Plot filesep  Deci.Plot.Behv.Acc.Title{fig} ' Behavioral Outputs' ]) == 2
            writematrix([],[Deci.Folder.Plot filesep  Deci.Plot.Behv.Acc.Title{fig} ' Behavioral Outputs' ],'FileType','spreadsheet','Sheet','TempSheet');
            %xls_delete_sheets([Deci.Folder.Plot filesep  Deci.Plot.Behv.Acc.Title{fig} ' Behavioral Outputs' ],'Accuracy_Summary');
        end
        
        excelAccdata = table(Sub.Acc(subs)',Deci.Plot.Behv.Acc.Subtitle{fig}(conds)',abs(AccBlock(blks)),trls',fAcc,fAccsem,'VariableNames',{'Subj' 'Cond' 'Blk'  'Trl' 'Accuracy' 'SEM'});
        writetable(excelAccdata,[Deci.Folder.Plot filesep Deci.Plot.Behv.Acc.Title{fig} ' Behavioral Outputs' ],'FileType','spreadsheet','Sheet','Accuracy_Summary');
        %xls_delete_sheets([Deci.Folder.Plot filesep  Deci.Plot.Behv.Acc.Title{fig} ' Behavioral Outputs' ],'TempSheet');
    end
    
    for fig = find(Deci.Plot.Behv.RT.Figure)
        
        %ExportExcel
        colnames = [];
        
        for blk = 1:size(fullRT{fig},3)
            for cond = 1:size(fullRT{fig},2)
                
                colnames{end+1} = [Deci.Plot.Behv.RT.Subtitle{fig}{cond} ' : Block ' num2str(abs(RTBlock(blk)))];
            end
        end
        
        exceldata = reshape(nanmean(fullRT{fig},4),[size(fullRT{fig},1) prod([size(fullRT{fig},2) size(fullRT{fig},3)])]);
        
        exceldata = horzcat(Deci.SubjectList',num2cell(exceldata));
        exceldata = vertcat([{'Reaction Time'} colnames],exceldata);
        
        if exist([Deci.Folder.Plot filesep  Deci.Plot.Behv.RT.Title{fig} ' Behavioral Outputs' ]) == 2
            writematrix([],[Deci.Folder.Plot filesep  Deci.Plot.Behv.RT.Title{fig} ' Behavioral Outputs' ],'FileType','spreadsheet','Sheet','TempSheet');
            %xls_delete_sheets([Deci.Folder.Plot filesep  Deci.Plot.Behv.RT.Title{fig} ' Behavioral Outputs' ],'RT_Full');
        end
        
        writecell(exceldata,[Deci.Folder.Plot filesep  Deci.Plot.Behv.RT.Title{fig} ' Behavioral Outputs' ],'FileType','spreadsheet','Sheet','RT_Full');
        %xls_delete_sheets([Deci.Folder.Plot filesep  Deci.Plot.Behv.RT.Title{fig} ' Behavioral Outputs' ],'TempSheet');
        
        
        
        subs= [];
        conds = [];
        blks = [];
        trls = [];
        
        
        
        for trl = 1:size(RT{fig},4)
            for blk = 1:size(RT{fig},3)
                for cond = 1:size(RT{fig},2)
                    for sub = 1:size(RT{fig},1)
                        subs(end+1) =  sub;
                        conds(end+1) = cond;
                        blks(end+1) = blk;
                        trls(end+1) = trl;
                    end
                end
            end
        end
        
        sRT = size(RT{fig});
        
        fRT = reshape(RT{fig},[prod([sRT]) 1]);
        fRTsem = reshape(RTsem{fig},[prod([sRT]) 1]);
        
        if exist([Deci.Folder.Plot filesep  Deci.Plot.Behv.RT.Title{fig} ' Behavioral Outputs' ]) == 2
            writematrix([],[Deci.Folder.Plot filesep  Deci.Plot.Behv.RT.Title{fig} ' Behavioral Outputs' ],'FileType','spreadsheet','Sheet','TempSheet');
            %xls_delete_sheets([Deci.Folder.Plot filesep  Deci.Plot.Behv.RT.Title{fig} ' Behavioral Outputs' ],'RT_Summary');
        end
        
        excelAccdata = table(Sub.RT(subs)',Deci.Plot.Behv.RT.Subtitle{fig}(conds)',abs(RTBlock(blks)),trls',fRT,fRTsem,'VariableNames',{'Subj' 'Cond' 'Blk'  'Trl' 'ReactionTime','SEM'});
        writetable(excelAccdata,[Deci.Folder.Plot filesep Deci.Plot.Behv.RT.Title{fig} ' Behavioral Outputs' ],'FileType','spreadsheet','Sheet','RT_Summary');
        %xls_delete_sheets([Deci.Folder.Plot filesep  Deci.Plot.Behv.RT.Title{fig} ' Behavioral Outputs' ],'TempSheet');
        
        
    end
end



%% plot

if ~isempty(Deci.Plot.Behv.Acc)
    
    if Deci.Plot.Behv.Acc.BaW
        Sub.Acc = {'SubjAvg'};
        
    end
    
    
    for fig = find(Deci.Plot.Behv.Acc.Figure)
        
        fignum = 1;
        
        for subj = 1:length(Sub.Acc)
            
            a(fignum,subj) = figure;
            a(fignum,subj).Visible = 'on';
            
            
            for draw = 1:size(Acc{fig},2)
                
                
                if ~Deci.Plot.Behv.Acc.BaW
                    
                    if length(find(size(squeeze(Acc{fig}(subj,draw,:,:))) ~= 1)) ==1
                        
                        top = squeeze(Acc{fig}(subj,draw,:,:))*100 + squeeze(Accsem{fig}(subj,draw,:,:))*100;
                        bot = squeeze(Acc{fig}(subj,draw,:,:))*100 - squeeze(Accsem{fig}(subj,draw,:,:))*100;
                        
                        shapes = [1:length(Acc{fig}(subj,draw,:,:)) length(Acc{fig}(subj,draw,:,:)):-1:1];
                        shapes(isnan([top' fliplr(bot')])) = nan;
                        
                        pgon = polyshape(shapes,[top' fliplr(bot')],'Simplify', false);
                        b = plot(pgon,'HandleVisibility','off');
                        hold on
                        b.EdgeAlpha = 0;
                        b.FaceAlpha = .15;
                        
                        h = plot(squeeze(Acc{fig}(subj,draw,:))*100);
                        h.Parent.YLim = [0 100];
                        h.Color = b.FaceColor;
                        h.LineWidth = 1;
                        hold on;
                        title(h.Parent,[Sub.Acc{subj} ' ' Deci.Plot.Behv.Acc.Title{fig}],'Interpreter','none');
                        
                        if ~Deci.Plot.Behv.Acc.Collapse.Block && length(find(size(squeeze(Acc{fig}(subj,draw,:,:))) ~= 1)) ==1
                            h.Parent.XLabel.String = 'Block #';
                            xticks(1:length(AccBlock))
                        elseif ~Deci.Plot.Behv.Acc.Collapse.Trial && length(find(size(squeeze(Acc{fig}(subj,draw,:,:))) ~= 1)) ==1
                            h.Parent.XLabel.String = 'Trial #';
                        end
                        
                        h.Parent.YLabel.String = 'Percent (%)';
                        legend(h.Parent,[ Deci.Plot.Behv.Acc.Subtitle{fig}]);
                        
                    elseif length(find(size(squeeze(Acc{fig}(subj,draw,:,:))) ~= 1)) ==2
                        
                        subplot(size(Acc{fig},2),1,draw)
                        h = imagesc(squeeze(Acc{fig}(subj,draw,:,:))*100);
                        h.Parent.CLim = [0 100];
                        h.Parent.YLabel.String = 'Block #';
                        h.Parent.XLabel.String = 'Trial #';
                        title(h.Parent,[Sub.Acc{subj} ' ' Deci.Plot.Behv.Acc.Subtitle{fig}{draw}],'Interpreter','none');
                        
                    else
                        disp(['Acc Total for ' Sub.Acc{subj} ' ' Deci.Plot.Behv.Acc.Title{fig} ' ' num2str(squeeze(Acc{fig}(subj,draw,:,:))*100) '%' ' +- ' num2str(squeeze(Accsem{fig}(subj,draw,:,:))*100)]);
                        
                        if draw == 1
                            CleanBars(Acc{fig}(subj,:,:,:),Accsem{fig}(subj,:,:,:))
                            ylim([0 1]);
                            title(['Acc Total for ' Sub.Acc{subj} ': ' Deci.Plot.Behv.Acc.Title{fig} ' '],'Interpreter','none')
                            legend([Deci.Plot.Behv.Acc.Subtitle{fig}])
                            xticklabels(Sub.Acc{subj})
                            ylabel('Accuracy (Percent)')
                            %disp(num2str(squeeze(Acc(subj,draw,:,:))*100))
                        end
                        
                    end
                    
                else
                    
                    if length(find(size(squeeze(Acc{fig}(subj,draw,:,:))) ~= 1)) ==2
                        h = subplot(size(Acc{fig},2),1,draw);
                        boxplot(squeeze(Acc{fig}(:,draw,:,:))*100);
                        h.CLim = [0 100];
                        
                        if ~Deci.Plot.Behv.Acc.Collapse.Block && length(find(size(squeeze(Acc{fig}(subj,draw,:,:))) ~= 1)) ==1
                            h.XLabel.String = 'Block #';
                            xticks(1:length(AccBlock))
                        elseif ~Deci.Plot.Behv.Acc.Collapse.Trial && length(find(size(squeeze(Acc{fig}(subj,draw,:,:))) ~= 1)) ==1
                            h.XLabel.String = 'Trial #';
                        end
                        
                        h.YLabel.String = 'Trial #';
                        title(h,[Sub.Acc{subj} ' ' Deci.Plot.Behv.Acc.Subtitle{fig}{draw}],'Interpreter','none');
                        hold on
                    else
                        h = gca;
                        boxplot(squeeze(Acc{fig})*100);
                        h.CLim = [0 100];
                        
                        title(['Acc Total for ' Sub.Acc{subj} ': ' Deci.Plot.Behv.Acc.Title{fig} ' '],'Interpreter','none')
                        xticklabels(Deci.Plot.Behv.Acc.Subtitle{fig})
                        ylabel('Accuracy (Percent)')
                        hold on
                    end
                    
                end
                
                
                
                
            end
            
            if length(find(size(squeeze(Acc{fig}(subj,draw,:,:))) ~= 1)) == 2
                suptitle([Sub.Acc{subj} ' ' Deci.Plot.Behv.Acc.Title{fig}]);
            end
            
            fignum = fignum + 1;
        end
        
        for lim = 1:numel(a)
            a(lim).Children(end).YLim = minmax(cell2mat(arrayfun(@(c) [c.Children(end).YLim],a(:)','UniformOutput',false)));
        end
        
    end
    
    
    
end


if ~isempty(Deci.Plot.Behv.RT) &&  ~isempty(find(Deci.Plot.Behv.RT.Figure))
    
    if Deci.Plot.Behv.RT.BaW
    Sub.RT = {'SubjAvg'};
   
end 
    
    for fig = find(Deci.Plot.Behv.RT.Figure)
        
        fignum = 1;
        
        for subj = 1:length(Sub.RT)
            
            
            d(fignum,subj) = figure;
            d(fignum,subj).Visible = 'on';
            
            
            
            for draw = 1:size(RT{fig},2)
                
                
             if ~Deci.Plot.Behv.RT.BaW
                if length(find(size(squeeze(RT{fig}(subj,draw,:,:))) ~= 1)) ==1
                    
                    
                    top = squeeze(RT{fig}(subj,draw,:,:)) + squeeze(RTsem{fig}(subj,draw,:,:));
                    bot = squeeze(RT{fig}(subj,draw,:,:)) - squeeze(RTsem{fig}(subj,draw,:,:));
                    
                    top = top(~isnan(top));
                    bot = bot(~isnan(bot));
                    
                    pgon = polyshape([1:length(find(~isnan(RT{fig}(subj,draw,:,:)))) length(find(~isnan(RT{fig}(subj,draw,:,:)))):-1:1],[top' fliplr(bot')],'Simplify', false);
                    b = plot(pgon,'HandleVisibility','off');
                    hold on
                    b.EdgeAlpha = 0;
                    b.FaceAlpha = .15;
                    
                    h = plot(squeeze(RT{fig}(subj,draw,:)));
                    h.Color = b.FaceColor;
                    h.LineWidth = 1;
                    hold on;
                    
                    title(h.Parent,[Sub.RT{subj} ' ' Deci.Plot.Behv.RT.Title{fig}],'Interpreter','none');
                    
                    if ~Deci.Plot.Behv.RT.Collapse.Block && length(find(size(squeeze(RT{fig}(subj,draw,:,:))) ~= 1)) ==1
                        h.Parent.XLabel.String = 'Block #';
                    elseif ~Deci.Plot.Behv.RT.Collapse.Trial && length(find(size(squeeze(RT{fig}(subj,draw,:,:))) ~= 1)) ==1
                        h.Parent.XLabel.String = 'Trial #';
                    end
                    
                    h.Parent.YLabel.String = 'Reaction Time (ms)';
                    legend(h.Parent,[Deci.Plot.Behv.RT.Subtitle{fig}]);
                    
                elseif length(find(size(squeeze(RT{fig}(subj,draw,:,:))) ~= 1)) ==2
                    
                    subplot(size(RT{fig},2),1,draw)
                    h = imagesc(squeeze(RT{fig}(subj,draw,:,:)));
                    h.Parent.YLabel.String = 'Block #';
                    h.Parent.XLabel.String = 'Trial #';
                    title(h.Parent,[Sub.RT{subj} ' ' Deci.Plot.Behv.RT.Subtitle{fig}{draw}],'Interpreter','none');
                    
                else
                    disp(['RT Total for ' Sub.RT{subj} ' ' Deci.Plot.Behv.RT.Subtitle{fig}{draw} ' ' num2str(squeeze(RT{fig}(subj,draw,:,:))) ' +- ' num2str(squeeze(RTsem{fig}(subj,draw,:,:)))])
                    
                    if draw == 1
                        CleanBars(RT{fig}(subj,:,:,:),RTsem{fig}(subj,:,:,:))
                        title(['RT Total for ' Sub.RT{subj} ': ' Deci.Plot.Behv.RT.Title{fig} ' '],'Interpreter','none')
                        legend([Deci.Plot.Behv.RT.Subtitle{fig}])
                        xticklabels(Sub.RT{subj})
                        ylabel('Reaction Time (ms)')
                    end
                end
             
             else
                 
                 if length(find(size(squeeze(RT{fig}(subj,draw,:,:))) ~= 1)) ==2
                     h = subplot(size(RT{fig},2),1,draw);
                     boxplot(squeeze(RT{fig}(:,draw,:,:))*100);
                     h.CLim = [0 100];
                     
                     if ~Deci.Plot.Behv.RT.Collapse.Block && length(find(size(squeeze(RT{fig}(subj,draw,:,:))) ~= 1)) ==1
                         h.XLabel.String = 'Block #';
                     elseif ~Deci.Plot.Behv.RT.Collapse.Trial && length(find(size(squeeze(RT{fig}(subj,draw,:,:))) ~= 1)) ==1
                         h.XLabel.String = 'Trial #';
                     end
                     
                     h.YLabel.String = 'Reaction Time (ms)';
                     title(h,[Sub.RT{subj} ' ' Deci.Plot.Behv.RT.Subtitle{fig}{draw}],'Interpreter','none');
                     hold on
                     
                     
                 else
                     h = gca;
                     boxplot(squeeze(RT{fig}));
                     h.CLim = [0 100];
                     
                     title(['Reaction Time Total for ' Sub.RT{subj} ': ' Deci.Plot.Behv.RT.Title{fig} ' '],'Interpreter','none')
                     xticklabels(Deci.Plot.Behv.RT.Subtitle{fig})
                     ylabel('Reaction Time (ms)')
                     hold on
                 end
             end
             
             
             
            end
            
            if length(find(size(squeeze(RT{fig}(subj,draw,:,:))) ~= 1)) ==2
                suptitle([Sub.RT{subj} ' ' Deci.Plot.Behv.RT.Title{fig}]);
            end
            
        end
        
        for lim = 1:numel(d)
            d(lim).Children(end).YLim = minmax(cell2mat(arrayfun(@(c) [c.Children(end).YLim],d(:)','UniformOutput',false)));
        end
    end
    
    fignum = fignum + 1;
    
end

end
