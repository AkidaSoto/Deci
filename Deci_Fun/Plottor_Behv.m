function Plottor_Behv(Deci)

for subject_list = 1:length(Deci.SubjectList)
    
    data = [];
    
    switch Deci.Plot.Behv.Source
        case 'PostArt'
            load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} '_info']);
            data.trialinfo = data.condinfo{1};
            data.event = data.condinfo{2};
        case 'Definition'
            load([Deci.Folder.Version filesep 'Definition' filesep Deci.SubjectList{subject_list}]);
            data = cfg;
            data.trialinfo = data.trl(:,end-length(Deci.DT.Locks)+1:end);
    end
    

    if ~isempty(Deci.Plot.Behv.RT)
        
        if length(Deci.DT.Locks) <2 || length(Deci.Plot.Behv.RT.Locks) ~= 2
            error('DT.Locks seems to not have enough to calculate Locks(2)-Locks(1)')
        end
        
        if isempty(Deci.Plot.Behv.RT.Block)
            Deci.Plot.Behv.RT.Block = {-1};
        end
        
        
        for blk = 1:length(Deci.Plot.Behv.RT.Block)
            for draw = 1:length(Deci.Plot.Behv.RT.Draw)
                
                
                if any(any(data.event < 0))
                    iblk = -1;
                else
                    iblk = 1;
                end
                
                draws = Deci.Analysis.Conditions(Deci.Plot.Behv.RT.Draw{draw});
                maxt = max(sum(ismember(data.event,[draws{:} iblk*Deci.Plot.Behv.RT.Block(blk)]),2));
                trl = sum(ismember(data.event,[draws{:}  iblk*Deci.Plot.Behv.RT.Block(blk)]),2) == maxt;
                
                if strcmpi(Deci.Plot.Behv.RT.Collapse.Uneven,'positional:nans')
                    RT{subject_list,draw,blk} = nan([size(data.trialinfo(trl,:),1) 1]);
                    RT{subject_list,draw,blk}(~isnan(mean(data.trialinfo(trl,:),2))) =  [-data.trialinfo(trl,Deci.Plot.Behv.RT.Locks(2)) - -data.trialinfo(trl,Deci.Plot.Behv.RT.Locks(1))];
                else
                    RT{subject_list,draw,blk} =  [-data.trialinfo(trl,Deci.Plot.Behv.RT.Locks(2)) - -data.trialinfo(trl,Deci.Plot.Behv.RT.Locks(1))];
                end
                
            end
        end
        
        %disp([Deci.SubjectList{subject_list} ' RT:']);
        %disp(RT(subject_list,:,:));
    end
    
    if ~isempty(Deci.Plot.Behv.Acc)
        Exist(Deci.Plot.Behv.Acc,'Total');
        Exist(Deci.Plot.Behv.Acc,'Subtotal');
        
        if isempty(Deci.Plot.Behv.Acc.Block)
            Deci.Plot.Behv.Acc.Block = {-1};
        end
        
        
        for blk = 1:length(Deci.Plot.Behv.Acc.Block)
            for draw = 1:length(Deci.Plot.Behv.Acc.Total)
                
                if any(any(data.event < 0))
                    iblk = -1;
                else
                    iblk = 1;
                end
                
                draws = Deci.Analysis.Conditions(Deci.Plot.Behv.Acc.Total{draw});
                maxt = max(sum(ismember(data.event,[draws{:} iblk*Deci.Plot.Behv.Acc.Block(blk)]),2));
                trl = sum(ismember(data.event,[draws{:}  iblk*Deci.Plot.Behv.Acc.Block(blk)]),2) == maxt;
                Total = data.event(trl,:);
                Trlpos = data.trialnum(trl,:);
                
                subdraws = Deci.Analysis.Conditions(Deci.Plot.Behv.Acc.Subtotal{draw});
                maxt = max(sum(ismember(data.event,[subdraws{:} Deci.Plot.Behv.Acc.Block(blk)]),2));
                subtrl = sum(ismember(Total,[subdraws{:} Deci.Plot.Behv.Acc.Block(blk)]),2) == maxt;
                
                
                if strcmpi(Deci.Plot.Behv.Acc.Collapse.Uneven,'positional:nans')
                    Acc{subject_list,draw,blk} = nan([length(Total) 1]);
                    Acc{subject_list,draw,blk}(~isnan(mean(data.trialinfo(trl,:),2))) =   subtrl;
                else
                    Acc{subject_list,draw,blk} = subtrl;
                end
                
                
                if Deci.Plot.Behv.Acc.Collapse.Movmean
                    Acc{subject_list,draw,blk} = movmean(Acc{subject_list,draw,blk},[length( Acc{subject_list,draw,blk}) 0]);
                end
                
            end
        end
        
        %disp([Deci.SubjectList{subject_list} ' Acc:']);
        %disp(Acc(subject_list,:,:));
    end
    
    
end

 if ~isempty(Deci.Plot.Behv.Acc)
     mkdir([Deci.Folder.Version filesep 'Plot'])
     save([Deci.Folder.Version filesep 'Plot' filesep 'Acc'],'Acc');
 end
 
 if ~isempty(Deci.Plot.Behv.RT)
     mkdir([Deci.Folder.Version filesep 'Plot'])
     save([Deci.Folder.Version filesep 'Plot' filesep 'RT'],'RT');
 end

%% sort
if ~isempty(Deci.Plot.Behv.Acc)
    if strcmpi(Deci.Plot.Behv.Acc.Collapse.Uneven,'maxlength:nans')
        
        for subject_list = 1:length(Deci.SubjectList)
            for draw = 1:size(Acc,2)
                fAcc(subject_list,draw,:,:) =  squeeze(cell2mat(cellfun(@(c) [c;nan(max(cellfun(@length,Acc(:)))-length(c),1)],Acc(subject_list,draw,:),'un',0)))';
            end
        end
        Acc = fAcc;
        Accstd = [];
    end
    
    if Deci.Plot.Behv.Acc.Collapse.Trial
        Accstd = nanstd(Acc,4);
        Acc = nanmean(Acc,4);
        
    end
    
    if Deci.Plot.Behv.Acc.Collapse.Block
        Accstd = nanstd(Acc,3);
        Acc = nanmean(Acc,3);
        
    end
    
    Sub.Acc = Deci.SubjectList;
    if Deci.Plot.Behv.Acc.Collapse.Subject
        Accstd =  nanstd(Acc,1);
        Acc =  nanmean(Acc,1);
        
        Sub.Acc = {'SubjAvg'};
    end
    save([Deci.Folder.Version filesep 'Plot' filesep 'RTstd'],'Accstd');    
end

if ~isempty(Deci.Plot.Behv.RT)
    if strcmpi(Deci.Plot.Behv.RT.Collapse.Uneven,'maxlength:nans')
        
        for subject_list = 1:length(Deci.SubjectList)
            for draw = 1:size(RT,2)
                fRT(subject_list,draw,:,:) =  squeeze(cell2mat(cellfun(@(c) [c;nan(max(cellfun(@length,RT(:)))-length(c),1)],RT(subject_list,draw,:),'un',0)))';
            end
        end
        RT = fRT;
        RTstd = [];
    end
    
    if Deci.Plot.Behv.RT.Collapse.Trial
        RTstd = nanstd(RT,4);
        RT = nanmean(RT,4);
    end
    
    if Deci.Plot.Behv.RT.Collapse.Block
        RTstd = nanstd(RT,3);
        RT = nanmean(RT,3);
    end
    
    Sub.RT = Deci.SubjectList;
    if Deci.Plot.Behv.RT.Collapse.Subject
        RTstd = nanstd(RT,1);
        RT =  nanmean(RT,1);
        Sub.RT = {'SubjAvg'};
    end
    save([Deci.Folder.Version filesep 'Plot' filesep 'RTstd'],'RTstd');
end


%% plot
if ~isempty(Deci.Plot.Behv.Acc)
    for subj = 1:size(Acc,1)
        
        if ~all([Deci.Plot.Behv.Acc.Collapse.Trial Deci.Plot.Behv.Acc.Collapse.Block])
        a = figure;
        a.Visible = 'on';
        end
        
        for draw = 1:size(Acc,2)
            
            
            if length(find(size(squeeze(Acc(subj,draw,:,:))) ~= 1)) ==1
                
                h = plot(squeeze(Acc(subj,draw,:))*100);
                h.Parent.YLim = [0 100];
                hold on;
                title(h.Parent,[Sub.Acc{subj} ' ' Deci.Plot.Behv.Acc.Title],'Interpreter','none');
                
                if ~Deci.Plot.Behv.Acc.Collapse.Block && length(find(size(squeeze(Acc(subj,draw,:,:))) ~= 1)) ==1
                    h.Parent.XLabel.String = 'Block #';
                elseif ~Deci.Plot.Behv.Acc.Collapse.Trial && length(find(size(squeeze(Acc(subj,draw,:,:))) ~= 1)) ==1
                    h.Parent.XLabel.String = 'Trial #';
                end
                
                h.Parent.YLabel.String = 'Percent (%)';
                legend(h.Parent,[ Deci.Plot.Behv.Acc.Subtitle]);
                
            elseif length(find(size(squeeze(Acc(subj,draw,:,:))) ~= 1)) ==2
                
                subplot(size(Acc,2),1,draw)
                h = imagesc(squeeze(Acc(subj,draw,:,:))*100);
                h.Parent.CLim = [0 100];
                h.Parent.YLabel.String = 'Block #';
                h.Parent.XLabel.String = 'Trial #';
                title(h.Parent,[Sub.Acc{subj} ' ' Deci.Plot.Behv.Acc.Subtitle{draw}],'Interpreter','none');
                
            else
                    disp(['Acc Total for ' Sub.Acc{subj} ' ' Deci.Plot.Behv.Acc.Title{draw} ' ' num2str(squeeze(Acc(subj,draw,:,:))*100) '%' ' +- ' num2str(squeeze(Accstd(subj,draw,:,:))*100)]);
                    %disp(num2str(squeeze(Acc(subj,draw,:,:))*100))

            end
            
            
        end
        
        if length(find(size(squeeze(Acc(subj,draw,:,:))) ~= 1)) == 2
            suptitle([Sub.Acc{subj} ' ' Deci.Plot.Behv.Acc.Title]);
        end
    end
    
    
end

if ~isempty(Deci.Plot.Behv.RT)
    for subj = 1:size(RT,1)
        
        if ~all([Deci.Plot.Behv.RT.Collapse.Trial Deci.Plot.Behv.RT.Collapse.Block])
        b = figure;
        b.Visible = 'on';
        end
        
        for draw = 1:size(RT,2)
            
            
            if length(find(size(squeeze(RT(subj,draw,:,:))) ~= 1)) ==1
                
                h = plot(squeeze(RT(subj,draw,:)));
                hold on;
                
                title(h.Parent,[Sub.RT{subj} ' ' Deci.Plot.Behv.RT.Title],'Interpreter','none');
                
                if ~Deci.Plot.Behv.RT.Collapse.Block && length(find(size(squeeze(RT(subj,draw,:,:))) ~= 1)) ==1
                    h.Parent.XLabel.String = 'Block #';
                elseif ~Deci.Plot.Behv.RT.Collapse.Trial && length(find(size(squeeze(RT(subj,draw,:,:))) ~= 1)) ==1
                    h.Parent.XLabel.String = 'Trial #';
                end
                
                h.Parent.YLabel.String = 'Reaction Time (ms)';
                legend(h.Parent,[Deci.Plot.Behv.RT.Subtitle]);
                
            elseif length(find(size(squeeze(RT(subj,draw,:,:))) ~= 1)) ==2
                
                subplot(size(RT,2),1,draw)
                h = imagesc(squeeze(RT(subj,draw,:,:)));
                h.Parent.YLabel.String = 'Block #';
                h.Parent.XLabel.String = 'Trial #';
                title(h.Parent,[Sub.RT{subj} ' ' Deci.Plot.Behv.RT.Subtitle{draw}],'Interpreter','none');
                
            else
                disp(['RT Total for ' Sub.RT{subj} ' ' Deci.Plot.Behv.RT.Subtitle{draw} ' ' num2str(squeeze(RT(subj,draw,:,:))) ' +- ' num2str(squeeze(RTstd(subj,draw,:,:)))])
            end
            
        end
        
        if length(find(size(squeeze(RT(subj,draw,:,:))) ~= 1)) ==2
            suptitle([Sub.RT{subj} ' ' Deci.Plot.Behv.RT.Title]);
        end
        
    end
    
    
end


end
