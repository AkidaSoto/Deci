function Plottor_scatter(Deci,params)


params = Exist(params,'std','std');
params = Exist(params,'type',@mean);

for subject_list = 1:length(Deci.SubjectList)
    
    data = [];
    
    switch Deci.Plot.Behv.Source
        case 'PostArt'
            info = [];
            load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'info');
            data = info;
            
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
                        else
                            eve = RTevents;
                        end
                    else
                        eve = RTevents;
                    end
                    eve = eve(find(any(ismember(eve,RTBlock(blk)),2)),:);
                    
                    
                    eveTotal = nan([1 length(find(any(ismember(eve,RTBlock(blk)),2)))]);
                    
                    
                    maxt = length(find(cellfun(@(c) any(ismember([draws{:}],c)), Deci.DT.Markers)));
                    trl = [sum(ismember(eve,[draws{:}]),2) == maxt];
                    
                    eveTotal(trl) =  [-data.trialinfo(trl,Deci.Plot.Behv.RT.Locks(2)) - -data.trialinfo(trl,Deci.Plot.Behv.RT.Locks(1))];
                    
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
            
            
            if isa(Deci.Plot.Behv.Acc.Block,'function_handle')
                data.event(:,find(data.event(1,:) < 1)) = Deci.Plot.Behv.Acc.Block(data.event(:,find(data.event(1,:) < 1)));
                
                AccBlock= unique(data.event(:,find(data.event(1,:) < 1)),'stable');
            elseif isempty(Deci.Plot.Behv.Acc.Block)
                AccBlock = {-1};
            end
            
            for blk = 1:length(AccBlock)
                for draw = 1:length(Deci.Plot.Behv.Acc.Total{fig})
                    
                    draws = Deci.Analysis.Conditions(Deci.Plot.Behv.Acc.Total{fig}{draw});
                    
                    if isfield(Deci.Plot.Behv,'Static')
                        
                        if sum(ismember(Deci.Plot.Behv.Static,[draws{:}])) == 1
                            
                            eve = data.event(logical(sum(ismember(data.event,Deci.Plot.Behv.Static(ismember(Deci.Plot.Behv.Static,[draws{:}]))),2)),:);
                        else
                            eve = data.event;
                        end
                    else
                        eve = data.event;
                    end
                    eve = eve(find(any(ismember(eve,AccBlock(blk)),2)),:);
                    
                    eveTotal = nan([1 length(find(any(ismember(eve,AccBlock(blk)),2)))]);
                    
                    maxt = length(find(cellfun(@(c) any(ismember([draws{:}],c)), Deci.DT.Markers)));
                    trl = [[sum(ismember(eve,[draws{:}]),2)] == maxt];
                    
                    eveTotal(trl) = 0;
                    
                    
                    subdraws = Deci.Analysis.Conditions(Deci.Plot.Behv.Acc.Subtotal{fig}{draw});
                    subtrl = [[sum(ismember(eve,[subdraws{:}]),2) ] == maxt];
                    eveTotal([subtrl & trl]) = 1;
                    
                    eveTotal(subtrl) = 1;
                    
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
                end
            end
            
        end
        
    end
end

%% sort

for fig = find(Deci.Plot.Behv.Acc.Figure)
    
    clear fAcc
    
    if ~isempty(Deci.Plot.Behv.Acc) && ~isempty(find(Deci.Plot.Behv.Acc.Figure))
        fullAcc{fig} = Acc{fig};
        
        Acc{fig} = nanmean(Acc{fig},4);
        if Deci.Plot.Behv.Acc.Collapse.Block
            Acc{fig} = nanmean(Acc{fig},3);
        end
        
        AccMean{fig} = params.type(Acc{fig},1);
        Accstd{fig} =  nanstd(Acc{fig},[],1);
        
        Sub.Acc = {'SubjAvg'}
        
    end
end

for fig = find(Deci.Plot.Behv.RT.Figure)
    clear fRT
    if ~isempty(Deci.Plot.Behv.RT) &&  ~isempty(find(Deci.Plot.Behv.RT.Figure))
        fullRT{fig} = RT{fig};

        RT{fig} = nanmean(RT{fig},4);
        if Deci.Plot.Behv.RT.Collapse.Block
            RT{fig} = nanmean(RT{fig},3);
        end
        
        
        RTstd{fig} = nanstd(RT{fig},[],1);
        RTMean{fig} =  params.type(RT{fig},1);
        Sub.RT = {'SubjAvg'};
        
        %save([Deci.Folder.Version filesep 'Plot' filesep 'SimRT'],'RT','fullRT');
    end
    
end
%% plot
% 
% if Deci.Whatever
%     
%     
%     for subject_list = 1:length(Deci.SubjectList)
%         
%         mkdir([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{subject_list}]);
%         save([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{subject_list}  filesep 'Acc.mat'],'Acc');
%         Deci.SubjectList{subject_list} = Deci.SubjectList{subject_list}([1:6 8:end]);
%         mkdir([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{subject_list}]);
%         save([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{subject_list}  filesep 'Acc.mat'],'Acc');
%         
%     end
% end
    
for fig = find(Deci.Plot.Behv.Acc.Figure)
    
    liner = lines(size(Acc{fig},2));
    if ~isempty(Deci.Plot.Behv.Acc)
        
        for draw = 1:size(Acc{fig},2)
            
            a = figure;
            a.Visible = 'on';
            
            scatty = reshape(Acc{fig}(:,draw,:),[size(Acc{fig},1)*size(Acc{fig},3) 1]);
            
            scatter(sort(repmat([1:size(Acc{fig},3)],[1 size(Acc{fig},1)])),reshape(Acc{fig}(:,draw,:),[size(Acc{fig},1)*size(Acc{fig},3) 1]),'MarkerEdgeColor',liner(draw,:));
            hold on
            
            dx = 0.1; dy = 0; % displacement so the text does not overlay the data points
            
            
                for subj = 1:size(Acc{fig},1)
                    for blk = 1:size(Acc{fig},3)
                        
                        %if Acc{fig}(subj,draw,blk) >= AccMean{fig}(:,draw,blk)+params.std*Accstd{fig}(:,draw,blk) || Acc{fig}(subj,draw,blk) <= AccMean{fig}(:,draw,blk)-params.std*Accstd{fig}(:,draw,blk)
                            text(blk+dx, Acc{fig}(subj,draw,blk)+dy, Deci.SubjectList(subj),'Interpreter','none','FontSize',9);
                       % end
                    end
                end

            
            
            for blk = 1:size(Acc{fig},3)
                plot([.75+blk-1 1.25+blk-1],[AccMean{fig}(:,draw,blk)+params.std*Accstd{fig}(:,draw,blk) AccMean{fig}(:,draw,blk)+params.std*Accstd{fig}(:,draw,blk)],'LineWidth',2,'Color',liner(draw,:),'HandleVisibility','off');
                plot([.75+blk-1 1.25+blk-1],[AccMean{fig}(:,draw,blk)-params.std*Accstd{fig}(:,draw,blk) AccMean{fig}(:,draw,blk)-params.std*Accstd{fig}(:,draw,blk)],'LineWidth',2,'Color',liner(draw,:),'HandleVisibility','off');
            end
            
            %Changed to params.stdSTD
            
            xticks(1:size(Acc{fig},3))
            title(['Std Scatter Total for ' Sub.Acc{1} ': ' Deci.Plot.Behv.Acc.Title{fig} ' '],'Interpreter','none')
            legend([Deci.Plot.Behv.Acc.Subtitle{fig}(draw)])
            xticklabels(abs(AccBlock))
            ylabel('Accuracy (Percent)')
            xlabel('Block #');
            %disp(num2str(squeeze(Acc(subj,draw,:,:))*100))
            
        end
        
    end
end

for fig = find(Deci.Plot.Behv.RT.Figure)
    liner = lines(size(RT{fig},2));
    
    if ~isempty(Deci.Plot.Behv.RT)
        
        for draw = 1:size(RT{fig},2)
            
            a = figure;
            a.Visible = 'on';
            
            scatty = reshape(RT{fig}(:,draw,:),[size(RT{fig},1)*size(RT{fig},3) 1]);
            
            scatter(sort(repmat([1:size(RT{fig},3)],[1 size(RT{fig},1)])),reshape(RT{fig}(:,draw,:),[size(RT{fig},1)*size(RT{fig},3) 1]),'MarkerEdgeColor',liner(draw,:));
            hold on
            
            dx = 0.1; dy = 0; % displacement so the text does not overlay the data points
            
            for subj = 1:size(RT{fig},1)
                for blk = 1:size(RT{fig},3)
                    
                    %if RT{fig}(subj,draw,blk) >= RTMean{fig}(:,draw,blk)+params.std*RTstd{fig}(:,draw,blk) || RT{fig}(subj,draw,blk) <= RTMean{fig}(:,draw,blk)-params.std*RTstd{fig}(:,draw,blk)
                         text(blk+dx, RT{fig}(subj,draw,blk)+dy, Deci.SubjectList(subj),'Interpreter','none','FontSize',9);
                     %end
                end
            end
            
            for blk = 1:size(RT{fig},3)
                plot([.75+blk-1 1.25+blk-1],[RTMean{fig}(:,draw,blk)+params.std*RTstd{fig}(:,draw,blk) RTMean{fig}(:,draw,blk)+params.std*RTstd{fig}(:,draw,blk)],'LineWidth',2,'Color',liner(draw,:),'HandleVisibility','off');
                plot([.75+blk-1 1.25+blk-1],[RTMean{fig}(:,draw,blk)-params.std*RTstd{fig}(:,draw,blk) RTMean{fig}(:,draw,blk)-params.std*RTstd{fig}(:,draw,blk)],'LineWidth',2,'Color',liner(draw,:),'HandleVisibility','off');
            end
            
            %Changed to params.stdSTD
            
            xticks(1:size(RT{fig},3))
            title(['Std Scatter Total for ' Sub.RT{1} ': ' Deci.Plot.Behv.RT.Title{fig} ' '],'Interpreter','none')
            legend([Deci.Plot.Behv.RT.Subtitle{fig}(draw)])
            xticklabels(abs(RTBlock))
            ylabel('RT (ms)')
            xlabel('Block #');
            %disp(num2str(squeeze(RT(subj,draw,:,:))*100))
            
        end
        
    end
end



end