
function RLModelPlot(Deci,params)

disp('Starting Model Stat')
QL = [];
mdl = repmat({[]},[1 2]);
PseudoR = [];
PR = repmat({[]},[1 2]);

PrE = repmat({[]},[1 2]);
EV = repmat({[]},[1 2]);

for subj = 1:length(Deci.SubjectList)
    for Lock = 1:length(Deci.Analysis.LocksTitle)
        for Cond = 1:length(Deci.Analysis.CondTitle)
            
            load([Deci.Folder.Analysis filesep 'Extra' filesep 'QL' filesep Deci.SubjectList{subj}  filesep Deci.Analysis.CondTitle{Cond}],'QL');
            for QLs = 1:length(QL)
                mdl{QLs}(subj,:,Cond) = QL{QLs};
            end
            
            load([Deci.Folder.Analysis filesep 'Extra' filesep 'PseudoR' filesep Deci.SubjectList{subj}  filesep Deci.Analysis.CondTitle{Cond}],'PseudoR');
            for PRs = 1:length(PseudoR)
                PR{PRs}(subj,:,Cond) = PseudoR{PRs};
            end
            
            load([Deci.Folder.Analysis filesep 'Extra' filesep 'Pe_all' filesep Deci.SubjectList{subj}  filesep Deci.Analysis.CondTitle{Cond}],'Pe_all');
            PrE{subj,Cond} = vertcat(Pe_all{:});

            
            load([Deci.Folder.Analysis filesep 'Extra' filesep 'Q_all' filesep Deci.SubjectList{subj}  filesep Deci.Analysis.CondTitle{Cond}],'Q_all');
            EV{subj,Cond} = vertcat(Q_all{:});
            
        end
        
    end
end

subs = [];
conds = [];

%%


exceldata = [];

for model = 1:length(mdl)
    exceldata = [exceldata reshape(permute(mdl{model},[1 3 2]),[size(permute(mdl{model},[1 3 2]),1)*size(permute(mdl{model},[1 3 2]),2) size(permute(mdl{model},[1 3 2]),3)])];
end

for cond = 1:size(mdl{model},3)
    for sub = 1:size(mdl{model},1)
        subs(end+1) =  sub;
        conds(end+1) = cond;
    end
end

excelmdldata = mat2cell([exceldata],[size(mdl{model},1)*size(mdl{model},3)], [ones([1 size(exceldata,2)])]);

mkdir([Deci.Folder.Plot]);

excelmdldata = table(excelmdldata{:},Deci.SubjectList(subs)', Deci.Analysis.CondTitle(conds)','VariableNames',{'m1_Q0Opt' 'm1_Q02Wor' 'm1_aR' 'm1_aP1' 'm1_b' 'm2_aR' 'm2_aP' 'm2_b' 'm3_a' 'm3_b' 'Sub' 'Cond'});
writetable(excelmdldata,[Deci.Folder.Plot filesep 'Modeloutputs'],'FileType','spreadsheet','Sheet','Model Parameters')


for modl = 1:size(mdl,2)
        data =   squeeze(nanmean(mdl{modl}(:,:,:),1));
        datastd =  squeeze(nanstd(mdl{modl}(:,:,:),[],1));
        a = figure; a.Visible = 'on';
        CleanBars(data,datastd);
        
        title(['Model' num2str(modl)]);
        legend(Deci.Analysis.CondTitle)
        xticklabels(params.ParamTitles{modl}(:))
       
end


%%

exceldata = [];

for PRs = 1:length(PR)
    tt = squeeze([PR{PRs}]);
    
    sub = repmat(1:size(tt,1),[1 size(tt,2)])';
    cond = repmat(1:size(tt,2),[size(tt,1) 1]);
    
    %         [PRstat(PRs,:),tbl,stats,terms] = anovan(tt(:),{sub(:),cond(:)},'varnames',{'Subject','Condition'},'display','off');
    exceldata = cat(2,exceldata,tt(:));
end

subs = [];
conds = [];

exceldata = [];

for model = 1:length(PR)
    exceldata = [exceldata reshape(permute(PR{model},[1 3 2]),[size(permute(PR{model},[1 3 2]),1)*size(permute(PR{model},[1 3 2]),2) size(permute(PR{model},[1 3 2]),3)])];
end

for cond = 1:size(PR{model},3)
    for sub = 1:size(PR{model},1)
        subs(end+1) =  sub;
        conds(end+1) = cond;
    end
end

excelPRdata =mat2cell([exceldata],[size(PR{model},1)*size(PR{model},3)], [ones([1 size(exceldata,2)])]);

excelPRdata = table(excelPRdata{:},Deci.SubjectList(subs)', Deci.Analysis.CondTitle(conds)','VariableNames',{'m1_PseudoR' 'm2_PseudoR' 'm3_PseudoR' 'Sub' 'Cond'});

writetable(excelPRdata,[Deci.Folder.Plot filesep 'Modeloutputs'],'FileType','spreadsheet','Sheet','Model Fitness')

for mdl = 1:size(PR,2)
    data =   squeeze(nanmean(PR{mdl}(:,:,:),1) - 1);
    datastd =  squeeze(nanstd(PR{mdl}(:,:,:)-1,[],1));
   a = figure; a.Visible = 'on';
   CleanBars(data',datastd');
    
    title(['Model' num2str(mdl) ' Pseudo-R Squared']);
    legend(Deci.Analysis.CondTitle)
end

%%

PredictionError= [];
ExpectedValue =[];



for subj = 1:size(PrE,1)
    for cond = 1:size(PrE,2)
        
        for mdl = 1:size(PrE{subj,cond},1)
            for blk = 1:size(PrE{subj,cond},2)
                
                PredictionError(subj,cond,blk,:,mdl) = nan([1 params.MaxTrials]);
                PredictionError(subj,cond,blk,1:length(PrE{subj,cond}{mdl,blk}),mdl) = PrE{subj,cond}{mdl,blk};
                
                ExpectedValue(subj,cond,blk,:,mdl) = nan([1 params.MaxTrials]);
                ExpectedValue(subj,cond,blk,1:length(PrE{subj,cond}{mdl,blk}),mdl) = EV{subj,cond}{mdl,blk};
                
            end
        end
    end
end

MPredictionError = permute(nanmean(PredictionError,4),[1 2 3 5 4]);
MExpectedValue = permute(nanmean(ExpectedValue,4),[1 2 3 5 4]);


% if params.CollapseSubjects
%    PredictionError = nanmean(PredictionError,1);
%    ExpectedValue = nanmean(ExpectedValue,1);
%    
%    Deci.SubjectList = {'GrandAverage'};
% end


for mdl = 1:size(PredictionError,5)
    fig = figure;
    fig.Visible = 'on';
    liner = lines(size(PredictionError,2)*2);
    for cond = 1:size(PredictionError,2)
        
        data =   permute(nanmean(nanmean(PredictionError(:,cond,:,:,mdl),3),1),[2 4 5 1 3]);
        datastd =  permute(nanstd(nanmean(PredictionError(:,cond,:,:,mdl),3),[],1),[2 4 5 1 3]);
        
        top = data + datastd;
        bot = data - datastd;
        
        pgon = polyshape([1:params.MaxTrials params.MaxTrials:-1:1],[top fliplr(bot)],'Simplify', false);
        b = plot(pgon,'HandleVisibility','off','FaceColor',liner(cond,:));
        hold on
        b.EdgeAlpha = 0;
        b.FaceAlpha = .15;
        
        a = plot(1:params.MaxTrials,data,'Color',liner(cond,:));
        a.LineWidth = 1;
        
    end
    title(['Model' num2str(mdl) ': Prediction Error']);
    axis tight
    hold on
    plot([1 params.MaxTrials], [0 0], 'k--','HandleVisibility','off') % hor. line
    plot([0 0], ylim, 'k--') % vert. l
    
    legend(Deci.Analysis.CondTitle)
    
    fig = figure;
    fig.Visible = 'on';
    
    
    for cond = 1:size(MExpectedValue,2)
        
        data =   permute(nanmean(nanmean(ExpectedValue(:,cond,:,:,mdl),3),1),[2 4 5 1 3]);
        datastd =  permute(nanstd(nanmean(ExpectedValue(:,cond,:,:,mdl),3),[],1),[2 4 5 1 3]);
        
        top = data + datastd;
        bot = data - datastd;
        
        pgon = polyshape([1:params.MaxTrials params.MaxTrials:-1:1],[top fliplr(bot)],'Simplify', false);
        b = plot(pgon,'HandleVisibility','off','FaceColor',liner(cond+2,:));
        hold on
        b.EdgeAlpha = 0;
        b.FaceAlpha = .15;
        
        a = plot(1:params.MaxTrials,data,'Color',liner(cond+2,:));
        a.LineWidth = 1;
        a.Color = b.FaceColor;
        
    end
    title(['Model' num2str(mdl) ': Expected Value']);
    axis tight
    hold on
    plot([1 params.MaxTrials], [0 0], 'k--') % hor. line
    plot([0 0], ylim, 'k--','HandleVisibility','off') % vert. l
    
    legend(Deci.Analysis.CondTitle)
    
end

MPredictionError = reshape(MPredictionError,[size(MPredictionError,1)*size(MPredictionError,2) size(MPredictionError,3)*size(MPredictionError,4)]);

excelPRdata =mat2cell([exceldata],[size(PR{model},1)*size(PR{model},3)], [ones([1 size(exceldata,2)])]);

excelPRdata = table(excelPRdata{:},Deci.SubjectList(subs)', Deci.Analysis.CondTitle(conds)','VariableNames',{'m1_PseudoR' 'm2_PseudoR' 'm3_PseudoR' 'Sub' 'Cond'});

writetable(excelPRdata,[Deci.Folder.Plot filesep 'Modeloutputs'],'FileType','spreadsheet','Sheet','Model Fitness')




disp('Finished Model Stat')

end

