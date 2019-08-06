function Plottor_Behv(Deci)

for subject_list = 1:length(Deci.SubjectList)
    
    data = [];
    
    switch Deci.Plot.Behv.Source
        case 'PostArt'
            load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} '_info']);
            data.trialinfo = data.condinfo{1};
            data.event = data.condinfo{2};
            data.trialnum = data.condinfo{3};
            
            data.full = data.preart(2:3);
            
        case 'Definition'
            load([Deci.Folder.Version filesep 'Definition' filesep Deci.SubjectList{subject_list}]);
            data = cfg;
            data.trialinfo = data.trl(:,end-length(Deci.DT.Locks)+1:end);
            
            data.full{1} = data.event;
            data.full{2} = data.trialnum;
    end
    
    
    for fig = 1:length(find(Deci.Plot.Behv.RT.Figure))
        if ~isempty(Deci.Plot.Behv.RT)
            
            if length(Deci.DT.Locks) <2 || length(Deci.Plot.Behv.RT.Locks) ~= 2
                error('DT.Locks seems to not have enough to calculate Locks(2)-Locks(1)')
            end
            
            if isempty(Deci.Plot.Behv.RT.Block)
                Deci.Plot.Behv.RT.Block = {-1};
            end
            
            
            for blk = 1:length(Deci.Plot.Behv.RT.Block)
                for draw = 1:length(Deci.Plot.Behv.RT.Draw{fig})
                    
                    
                    if any(any(data.event < 0))
                        iblk = -1;
                    else
                        iblk = 1;
                    end
                    
                    draws = Deci.Analysis.Conditions(Deci.Plot.Behv.RT.Draw{fig}{draw});
                    maxt = max(sum(ismember(data.event,[draws{:}]),2));
                    trl = [sum(ismember(data.event,[draws{:}]),2) == maxt] & any(ismember(data.event,iblk*Deci.Plot.Behv.RT.Block(blk)),2);
                    
                    if strcmpi(Deci.Plot.Behv.RT.Collapse.Uneven,'positional:nans')
                        
                        mint =  [[sum(ismember(data.full{1},[draws{:}]),2) + sum(isnan(data.event),2)] == maxt] & any(ismember(data.full{1},iblk*Deci.Plot.Behv.RT.Block(blk)),2);
                        mintrl = data.full{2}(mint);
                        
                        RT{fig}{subject_list,draw,blk} = nan([length(find(mintrl)) 1]);
                        RT{fig}{subject_list,draw,blk}(ismember(mintrl,find(trl))) =  [-data.trialinfo(trl,Deci.Plot.Behv.RT.Locks(2)) - -data.trialinfo(trl,Deci.Plot.Behv.RT.Locks(1))];
                    else
                        RT{fig}{subject_list,draw,blk} =  [-data.trialinfo(trl,Deci.Plot.Behv.RT.Locks(2)) - -data.trialinfo(trl,Deci.Plot.Behv.RT.Locks(1))];
                    end
                    
                end
            end
            
            %disp([Deci.SubjectList{subject_list} ' RT:']);
            %disp(RT(subject_list,:,:));
        end
    end
    
    for fig = 1:length(find(Deci.Plot.Behv.Acc.Figure))
        
        if ~isempty(Deci.Plot.Behv.Acc)
            Exist(Deci.Plot.Behv.Acc,'Total');
            Exist(Deci.Plot.Behv.Acc,'Subtotal');
            
            if isempty(Deci.Plot.Behv.Acc.Block)
                Deci.Plot.Behv.Acc.Block = {-1};
            end
            
            
            for blk = 1:length(Deci.Plot.Behv.Acc.Block)
                for draw = 1:length(Deci.Plot.Behv.Acc.Total{fig})
                    
                    if any(any(data.event < 0))
                        iblk = -1;
                    else
                        iblk = 1;
                    end
                    
                    draws = Deci.Analysis.Conditions(Deci.Plot.Behv.Acc.Total{fig}{draw});
                    maxt = max(sum(ismember(data.event,[draws{:}]),2));
                    trl = [sum(ismember(data.event,[draws{:}]),2) == maxt] & any(ismember(data.event,iblk*Deci.Plot.Behv.Acc.Block(blk)),2);
                    
                    Total = data.event(trl,:);
                    
                    subdraws = Deci.Analysis.Conditions(Deci.Plot.Behv.Acc.Subtotal{fig}{draw});
                    maxt2 = max(sum(ismember(data.event,[subdraws{:}]),2));
                    subtrl = [sum(ismember(Total,[subdraws{:}]),2) == maxt2] & any(ismember(Total,iblk*Deci.Plot.Behv.Acc.Block(blk)),2);
                    
                    if strcmpi(Deci.Plot.Behv.Acc.Collapse.Uneven,'positional:nans')

                        %+ sum(isnan(data.event),2)
                        
                        mint =  [[sum(ismember(data.full{1},[draws{:}]),2) ] == maxt] &  any(ismember(data.full{1},[iblk*Deci.Plot.Behv.Acc.Block(blk)]),2);
                        mintrl = data.full{2}(mint);
                        
                        Acc{fig}{subject_list,draw,blk} = nan([length(find(mintrl)) 1]);
                        Acc{fig}{subject_list,draw,blk}(ismember(mintrl,find(trl))) =  subtrl;
                        
                    else
                        Acc{fig}{subject_list,draw,blk} = subtrl;
                    end
                    
                    
                    if Deci.Plot.Behv.Acc.Collapse.Movmean
                        Acc{fig}{subject_list,draw,blk} = movmean(Acc{fig}{subject_list,draw,blk},[length( Acc{fig}{subject_list,draw,blk}) 0],'omitnan');
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
    clear fAcc
    clear fRT
for fig = 1:length(find(Deci.Plot.Behv.Acc.Figure))
    
    if ~isempty(Deci.Plot.Behv.Acc) && ~isempty(find(Deci.Plot.Behv.Acc.Figure))
        if ismember(Deci.Plot.Behv.Acc.Collapse.Uneven,{'maxlength:nans','positional:nans'})
            
            for subject_list = 1:length(Deci.SubjectList)
                for draw = 1:size(Acc{fig},2)
                    fAcc(subject_list,draw,:,:) =  squeeze(cell2mat(cellfun(@(c) [c;nan(max(cellfun(@length,Acc{fig}(:)))-length(c),1)],Acc{fig}(subject_list,draw,:),'un',0)))';
                end
            end
            Acc{fig} = fAcc;
            Accsem{fig} = [];
        end
        
        if Deci.Plot.Behv.Acc.Collapse.Trial
            Accsem{fig} = nanstd(Acc{fig},4)/sqrt(size(Acc{fig},4));
            Acc{fig} = nanmean(Acc{fig},4);
            
        end
        
        if Deci.Plot.Behv.Acc.Collapse.Block
            Accsem{fig} = nanstd(Acc{fig},3)/sqrt(size(Acc{fig},3));
            Acc{fig} = nanmean(Acc{fig},3);
            
        end
        
        Sub.Acc = Deci.SubjectList;
        if Deci.Plot.Behv.Acc.Collapse.Subject
            Accsem{fig} =  nanstd(Acc{fig},1)/sqrt(size(Acc{fig},1));
            Acc{fig} =  nanmean(Acc{fig},1);
            
            Sub.Acc = {'SubjAvg'};
        end
        save([Deci.Folder.Version filesep 'Plot' filesep 'AccSem'],'Accsem');
    end
end

for fig = 1:length(find(Deci.Plot.Behv.RT.Figure)) 
    if ~isempty(Deci.Plot.Behv.RT) &&  ~isempty(find(Deci.Plot.Behv.RT.Figure))
        if ismember(Deci.Plot.Behv.RT.Collapse.Uneven,{'maxlength:nans','positional:nans'})
            
            for subject_list = 1:length(Deci.SubjectList)
                for draw = 1:size(RT{fig},2)
                    fRT(subject_list,draw,:,:) =  squeeze(cell2mat(cellfun(@(c) [c;nan(max(cellfun(@length,RT{fig}(:)))-length(c),1)],RT{fig}(subject_list,draw,:),'un',0)))';
                end
            end
            RT{fig} = fRT;
            RTsem = [];
        end
        
        if Deci.Plot.Behv.RT.Collapse.Trial
            RTsem{fig} = nanstd(RT{fig},4)/sqrt(size(RT{fig},4));
            RT{fig} = nanmean(RT{fig},4);
        end
        
        if Deci.Plot.Behv.RT.Collapse.Block
            RTsem{fig} = nanstd(RT{fig},3)/sqrt(size(RT{fig},3));
            RT{fig} = nanmean(RT{fig},3);
        end
        
        Sub.RT = Deci.SubjectList;
        if Deci.Plot.Behv.RT.Collapse.Subject
            RTsem{fig} = nanstd(RT{fig},1)/sqrt(size(RT{fig},1));
            RT{fig} =  nanmean(RT{fig},1);
            Sub.RT = {'SubjAvg'};
        end
        save([Deci.Folder.Version filesep 'Plot' filesep 'RTSem'],'RTsem');
    end
    
end
%% plot


for fig = 1:length(find(Deci.Plot.Behv.Acc.Figure))
    
    
    if ~isempty(Deci.Plot.Behv.Acc)
        for subj = 1:size(Acc{fig},1)
            
            a = figure;
            a.Visible = 'on';
            
            
            for draw = 1:size(Acc{fig},2)
                
                
                if length(find(size(squeeze(Acc{fig}(subj,draw,:,:))) ~= 1)) ==1
                    
                    top = squeeze(Acc{fig}(subj,draw,:,:))*100 + squeeze(Accsem{fig}(subj,draw,:,:))*100;
                    bot = squeeze(Acc{fig}(subj,draw,:,:))*100 - squeeze(Accsem{fig}(subj,draw,:,:))*100;
                    
                    pgon = polyshape([1:length(Acc{fig}(subj,draw,:,:)) length(Acc{fig}(subj,draw,:,:)):-1:1],[top' fliplr(bot')],'Simplify', false);
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
                    title(h.Parent,[Sub.Acc{subj} ' ' Deci.Plot.Behv.Acc.Subtitle{draw}],'Interpreter','none');
                    
                else
                    disp(['Acc Total for ' Sub.Acc{subj} ' ' Deci.Plot.Behv.Acc.Title{fig} ' ' num2str(squeeze(Acc{fig}(subj,draw,:,:))*100) '%' ' +- ' num2str(squeeze(Accsem{fig}(subj,draw,:,:))*100)]);
                    
                    if draw == 1
                        CleanBars(Acc{fig}(subj,:,:,:),Accsem{fig}(subj,:,:,:))
                        title(['Acc Total for ' Sub.Acc{subj} ': ' Deci.Plot.Behv.Acc.Title{fig} ' '],'Interpreter','none')
                        legend([Deci.Plot.Behv.Acc.Subtitle{draw}])
                        xticklabels(Sub.Acc{subj})
                        ylabel('Accuracy (Percent)')
                        %disp(num2str(squeeze(Acc(subj,draw,:,:))*100))
                    end
                    
                end
                
                
            end
            
            if length(find(size(squeeze(Acc{fig}(subj,draw,:,:))) ~= 1)) == 2
                suptitle([Sub.Acc{subj} ' ' Deci.Plot.Behv.Acc.Title{fig}]);
            end
        end
        
        
    end
end


for fig = 1:length(find(Deci.Plot.Behv.RT.Figure))
    
    if ~isempty(Deci.Plot.Behv.RT) &&  ~isempty(find(Deci.Plot.Behv.RT.Figure))
        for subj = 1:size(RT{fig},1)
            
            
            b = figure;
            b.Visible = 'on';
            
            
            for draw = 1:size(RT{fig},2)
                
                
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
                    title(h.Parent,[Sub.RT{subj} ' ' Deci.Plot.Behv.RT.Subtitle{draw}],'Interpreter','none');
                    
                else
                    disp(['RT Total for ' Sub.RT{subj} ' ' Deci.Plot.Behv.RT.Subtitle{fig}{draw} ' ' num2str(squeeze(RT{fig}(subj,draw,:,:))) ' +- ' num2str(squeeze(RTsem{fig}(subj,draw,:,:)))])
                    
                    if draw == 1
                        CleanBars(RT{fig}(subj,:,:,:),RTsem{fig}(subj,:,:,:))
                        title(['Acc Total for ' Sub.Acc{subj} ': ' Deci.Plot.Behv.RT.Title{fig} ' '],'Interpreter','none')
                        legend([Deci.Plot.Behv.RT.Subtitle{draw}])
                        xticklabels(Sub.RT{subj})
                        ylabel('Reaction Time (ms)')
                    end
                end
                
            end
            
            if length(find(size(squeeze(RT{fig}(subj,draw,:,:))) ~= 1)) ==2
                suptitle([Sub.RT{subj} ' ' Deci.Plot.Behv.RT.Title{fig}]);
            end
            
        end
        
        
    end
end

end
