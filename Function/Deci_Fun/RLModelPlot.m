
function RLModelPlot(Deci,params)

disp('Starting Model Stat')
QL = [];
mdl = repmat({[]},[1 2]);
PseudoR = [];
PR = repmat({[]},[1 2]);
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
            
        end
        
    end
end

subs = [];
conds = [];

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

excelmdldata = mat2cell([exceldata],[size(mdl{model},1)*2], [ones([1 size(exceldata,2)])]);

excelmdldata = table(excelmdldata{:},Deci.SubjectList(subs)', Deci.Analysis.CondTitle(conds),'VariableNames',{'m1_Q0Opt' 'm1_Q02Wor' 'm1_aR' 'm1_aP1' 'm1_b' 'm2_aR' 'm2_aP' 'm2_b' 'm3_a' 'm3_b' 'Sub' 'Cond'});
writetable(excelmdldata,[Deci.Folder.Plot filesep 'Modeloutputs'],'FileType','spreadsheet','Sheet','Model Parameters')

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

excelPRdata =mat2cell([exceldata],[size(PR{model},1)*2], [ones([1 size(exceldata,2)])]);

excelPRdata = table(excelPRdata{:},Deci.SubjectList(subs)', Deci.Analysis.CondTitle(conds),'VariableNames',{'m1_PseudoR' 'm2_PseudoR' 'm3_PseudoR' 'Sub' 'Cond'});

writetable(excelPRdata,[Deci.Folder.Plot filesep 'Modeloutputs'],'FileType','spreadsheet','Sheet','Model Fitness')

disp('Finished Model Stat')

end

