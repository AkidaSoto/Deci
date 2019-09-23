function Deci_Displacement(Deci,Params)


for subject_list = 1:length(Deci.SubjectList)
    load([Deci.Folder.Version filesep 'Definition' filesep Deci.SubjectList{subject_list}]);
    data = cfg;
    data.trialinfo = data.trl(:,end-length(Deci.DT.Locks)+1:end);

    
    if isempty(Params.Block)
        Deci.Plot.Behv.RT.Block = {-1};
    end
   
    
    combn = [];
    combntitles =[];
    
    for dura = 1:length([Params.RepDim])
        
        mrks = Params.Markers(Params.RepDim{dura});
        mrktitles = Params.MarkersTitle(Params.RepDim{dura});
        
        mrks = CombVec(mrks{:});
        mrktitles = CombVec(mrktitles{:})';
        
        combn{dura} = mrks+1000*[dura-1];
        combntitles{dura} = arrayfun(@(c) strjoin(mrktitles(c,:),'-'),1:size(mrktitles,1),'UniformOutput',false);
    end
    
    combn = CombVec(combn{:})';
    combntitles = CombVec(combntitles{:})';
    combntitles = arrayfun(@(c) strjoin(combntitles(c,:),', '),1:size(combntitles,1),'UniformOutput',false)';
    
    
    for blk = 1:length(Params.Block)
        
        for stat = 1:length(Params.Static)
            
            
            statictrialinfo = cfg.event([logical(sum(ismember(cfg.event,Params.Static(stat)),2)) & any(ismember(cfg.event,-1*Params.Block(blk)),2)],:);
            statictrl = cfg.trl([logical(sum(ismember(cfg.event,Params.Static(stat)),2)) & any(ismember(cfg.event,-1*Params.Block(blk)),2)],:);
            
            for draw = 1:size(combn,1)
                
                mrks = [Params.Static(stat) combn(draw,:)];
                
                maxt = max(sum(ismember(data.event,[mrks(:)]),2));
                
                trl = [sum(ismember(statictrialinfo,[mrks(:)]),2) == maxt];
                
                
                for dura =  1:length([Params.RepDim])
                    
                    temptrl = find(trl) + dura - 1;
                    
                    checktrl = find(trl) + length([Params.RepDim])-1;
                    
                    temptrl = temptrl(checktrl > 0 & checktrl <=  size(trl,1));
                    
                    temptrl = temptrl(temptrl > 0 & temptrl <=  size(trl,1));
                    
                    RT(subject_list,blk,stat,draw,dura) =  nanmean([-statictrl(temptrl,Params.Locks(2)+3) - -statictrl(temptrl,Params.Locks(1)+3)]);
                end
                    RTcount(subject_list,blk,stat,draw) = numel(temptrl);
                
            end
            
        end
    end
    
end

% RT Dimensions are Subj_Blk_Condition_SubCondition_Durations

if Params.Collapse.Block
    RT = nanmean(RT,2);
    RTcount = nanmean(RTcount,2);
    blktitle = {'all'};
else
    blktitle = strsplit(num2str(1:size(RT,2)),' ');
end

RTsz = length(size(RT));

MeanRT = permute(nanmean(RT,1),[2:RTsz 1]);
SemRT = permute(nanstd(RT,1),[2:RTsz 1])/sqrt(size(RT,1));

MeanRTcount = permute(nanmean(RTcount,1),[2:RTsz 1]);
StdRTcount = permute(nanstd(RTcount,1),[2:RTsz 1])/sqrt(size(RTcount,1));


% One Color for each Subcondition.
cmap = jet(size(MeanRT,3))/1.5;


% MeanRT Dimensions are Blk_Conditions_SubConditions_Duration
% MeanRTCount Dimensions are Blk_Conditions_Subconditions

for blk = 1:size(MeanRT,1)
    
    for stat = 1:size(MeanRT,2)
        h(blk,stat) = figure;
        hold on
        h(blk,stat).Visible = 'on';
        
        for draw = 1:size(MeanRT,3)
            errorbar(squeeze(MeanRT(blk,stat,draw,:))',squeeze(SemRT(blk,stat,draw,:))','Color', cmap(draw, :));
        end
        
        title([Params.StaticTitle{stat} ' blk ' blktitle{blk}]);
        legend(arrayfun(@(a,b,c)  [a{1} ' [trl ' num2str(b) '+-' num2str(c) ']'],combntitles,squeeze(MeanRTcount(blk,stat,:)),squeeze(StdRTcount(blk,stat,:)),'un',0));
        
        
        xLim = xlim;
        xticks([1:size(MeanRT,4)]);
        xticklabels(arrayfun(@(c) ['n+' num2str(c-1)],1:size(MeanRT,4),'UniformOutput',false));
        
        xlim([xLim(1)*.9 xLim(2)*1.1]);
        
        ylabel('reaction time (ms)');
    end
    
    for draw = 1:size(MeanRT,3)
        
        g(blk,draw) = figure;
        g(blk,draw).Visible = 'on';
        errorbar(squeeze(MeanRT(blk,:,draw,:))',squeeze(SemRT(blk,:,draw,:))');
        
        title([combntitles{draw}  '- blk ' blktitle{blk}]);
        
        legend(arrayfun(@(a,b,c)  [a{1} ' [trl ' num2str(b) '+-' num2str(c) ']'],[Params.StaticTitle],squeeze(MeanRTcount(blk,:,draw)),squeeze(StdRTcount(blk,:,draw)),'un',0));
        xLim = xlim;
        
        xticks([1:size(MeanRT,4)]);
        xticklabels(arrayfun(@(c) ['n+' num2str(c-1)],1:size(MeanRT,4),'UniformOutput',false));
        
        xlim([xLim(1)*.9 xLim(2)*1.1]);
        
        ylabel('time (ms)');
    end
    
    sMeanRT = squeeze(nanmean(diff(RT(:,blk,:,:,:),[],5),1));
    sSemRT = squeeze(nanstd(diff(RT(:,blk,:,:,:),[],5),1))/sqrt(size(RT(:,blk,:,:,:),1));
    
    for stat = 1:size(sMeanRT,1)
        
        
        f(blk,stat) = figure;
        f(blk,stat).Visible = 'on';
        CleanBars(squeeze(sMeanRT(stat,:,:)),squeeze(sSemRT(stat,:,:)));
        title([Params.StaticTitle{stat} ' Difference' ' blk ' blktitle{blk}]);
        
        xticks([1:size(sMeanRT,3)])
        legend(arrayfun(@(c) ['[n+' num2str(c) ' - n+'  num2str(c-1) ']'],1:size(MeanRT,4)-1,'UniformOutput',false));
        
        xticklabels(combntitles);
        ylabel('delta time (ms)');
        
    end
    
end

for lim = 1:numel(h)
   h(lim).Children(end).YLim = minmax(cell2mat(arrayfun(@(c) [c.Children(end).YLim],h,'UniformOutput',false))); 
end

for lim = 1:numel(f)
   f(lim).Children(end).YLim = minmax(cell2mat(arrayfun(@(c) [c.Children(end).YLim],f,'UniformOutput',false))); 
end

for lim = 1:numel(h)
   h(lim).Children(end).YLim = minmax(cell2mat(arrayfun(@(c) [c.Children(end).YLim],h,'UniformOutput',false))); 
end

end
