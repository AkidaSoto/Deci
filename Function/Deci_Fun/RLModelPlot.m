
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

stat = [];

exceldata = [];

for mdls = 1:length(mdl)
    
    excelparam = [];
    
    for param = 1:size(mdl{mdls},2)
        
        
        
        tt = squeeze([mdl{mdls}(:,param,:)]);
        
        sub = repmat(1:size(tt,1),[1 size(tt,2)])';
        cond = repmat(1:size(tt,2),[size(tt,1) 1]);
        
        excelparam(:,param) = tt(:);
        
    end
    
    exceldata = cat(2,exceldata,excelparam);
    
end

excelmdldata = mat2cell([exceldata sub(:) cond(:)],[length(sub(:))],[ones([1 size(exceldata,2)+2])]);

excelmdldata = table(excelmdldata{:},'VariableNames',{'m1_Q0Opt' 'm1_Q02Wor' 'm1_a' 'm1_b' 'm2_Q0Opt' 'm2_Q0Wor' 'm2_aR' 'm2_aP' 'm2_b' 'm3_a' 'm3_b' 'Sub' 'Cond'});

writetable(excelmdldata,[Deci.Folder.Plot filesep 'Modeloutputs'],'FileType','spreadsheet','Sheet','Model Parameters')

exceldata = [];

for PRs = 1:length(PR)
    tt = squeeze([PR{PRs}]);
    
    sub = repmat(1:size(tt,1),[1 size(tt,2)])';
    cond = repmat(1:size(tt,2),[size(tt,1) 1]);
    
    %         [PRstat(PRs,:),tbl,stats,terms] = anovan(tt(:),{sub(:),cond(:)},'varnames',{'Subject','Condition'},'display','off');
    exceldata = cat(2,exceldata,tt(:));
end


excelPRdata =mat2cell([exceldata sub(:) cond(:)],[length(sub(:))],[ones([1 size(exceldata,2)+2])]);

excelPRdata = table(excelPRdata{:},'VariableNames',{'m1_PseudoR' 'm2_PseudoR' 'm3_PseudoR' 'Sub' 'Cond'});

writetable(excelPRdata,[Deci.Folder.Plot filesep 'Modeloutputs'],'FileType','spreadsheet','Sheet','Model Fitness')

disp('Finished Model Stat')

end

