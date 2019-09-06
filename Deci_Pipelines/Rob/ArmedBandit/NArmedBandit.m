%% 1.) General Info,
%You run this once before you run any other step

Deci = [];
%Deci.Folder.Raw         = ['/Users/RobReinhart/Documents/Data/Behavioral'];              % Raw Data Folder Directory
Deci.Folder.Raw         = ['C:\Users\User\Desktop\John\AB2\Behav'];

Deci.SubjectList        = 'gui';                                                  % Cell Array of strings, 'all' or 'gui'( for subject selection)                                          % Which Step to implement 1-TD, 2-PP, 3-Art, 4-Analysis, 5-Plot
%Deci.Folder.Version     = ['/Users/RobReinhart/Documents/Data/Processed'];        % Output Folder Directory
Deci.Folder.Version     = ['C:\Users\User\Desktop\John\AB2\ProcessedData'];

Deci.Folder.Definition   = [Deci.Folder.Version filesep 'Definition'];
Deci.Folder.Plot         = [Deci.Folder.Version filesep 'Plot'];
Deci.Folder.Analysis     = [Deci.Folder.Version filesep 'Analysis'];



% Finds all Files in the Behavioral Folder
if strcmp(Deci.SubjectList,'all')
    Deci.SubjectList = cellfun(@(c) strsplit(c,'.'),CleanDir(Deci.Folder.Raw),'un',0);
    Deci.SubjectList = unique(cellfun(@(c) c{1},Deci.SubjectList,'un',0));
elseif strcmp(Deci.SubjectList,'gui')
    Deci.SubjectList = cellfun(@(c) strsplit(c,'.'),CleanDir(Deci.Folder.Raw),'un',0);
    Deci.SubjectList = unique(cellfun(@(c) c{1},Deci.SubjectList,'un',0));
    
    fakeUI = figure;
    fakeUI.UserData = Deci.SubjectList;
    fakeUI.Visible =  'off';
    select_labels(fakeUI,[],Deci.SubjectList);
    waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
    Deci.SubjectList = fakeUI.UserData;
    close(fakeUI);
end


% Check to make sure all subject has a mat fle
for subject_list = 1:length(Deci.SubjectList)
files_ending = {'.mat'};
for file_ending = 1:length(files_ending)
    if ~exist([Deci.Folder.Raw filesep Deci.SubjectList{subject_list} files_ending{file_ending}],'file')
        error([Deci.SubjectList{subject_list} ' does not have files for ' files_ending{file_ending} ' in raw data folder']);
    end
end
end


% Trial definition infomaton
Deci.DT.Type = 'Manual';                                                                        % Change if need alternative file, otherwise 'Manual' (expfunor2)
Deci.DT.Starts     = {9};                                                                       % Cell Array of Markers for Start codes.
Deci.DT.Ends       = {10};                                                                      % Cell Array of Markers for End codes.
Deci.DT.Markers    = {[21 23] [27 28] [31 32] [51 52]};                                   % Cell Array of Markers for Trial Defs
Deci.DT.Locks      = [14 30 50];                                                                % Num Array for each timelock (usually Stim, Rsp and Fdb Onset )
Deci.DT.Toi        = [-2 3];                                                                    % Time of Interest, be sure to include larger window for freq
Deci.DT.Block.Start   = {11};                                                                     % Cell Array of Markers for Block Starts
Deci.DT.Block.End     = {12};                                                                       % Cell Array of Markers for Block Starts
Deci.DT.Block.Markers = [];      

% Just for saving plot info
Deci.Analysis.Conditions    = {[21] [23]};
Deci.Analysis.LocksTitle = {'None'};
info.Lock = 1;
Deci.Analysis.CondTitle = {'Correct' 'Incorrect'}';


%% 2.) Run Trial Definition

% Run through matfiles and find the trials

for subject_list = 1:length(Deci.SubjectList)
    
    cfg = [];
    cfg.dataset = [Deci.Folder.Raw filesep Deci.SubjectList{subject_list} files_ending{1}];
    cfg.DT = Deci.DT;
    cfg.file_ending = files_ending{1};
    

    cfg.trialfun = 'expfunor2';
    cfg.Raw = Deci.Folder.Version;
    cfg.Subject = Deci.SubjectList{subject_list};
    evalc('cfg = ft_definetrial(cfg);'); % will return cfg.trl, the segmented data
    
    cfg.trialnum = cfg.trl(:,4);
    cfg.trl = cfg.trl(:,[1:3 5:end]);
    
    trllength = num2str(length(find(~isnan(mean(cfg.trl,2)))));
    disp(['Found ' num2str(trllength) ' trials out of ' num2str(size(cfg.trl,1)) ' for ' Deci.SubjectList{subject_list}]);
    
    mkdir([Deci.Folder.Definition]);
    disp('Saving Preprocessing Data');
    save([Deci.Folder.Definition filesep Deci.SubjectList{subject_list}],'cfg')
    
end
    
    

%% 3.) Run Models
% Only need to run this once


        for subject_list = 1:length(Deci.SubjectList)
            
            load([Deci.Folder.Definition filesep Deci.SubjectList{subject_list}],'cfg')
    
            info.condinfo = {cfg.trl cfg.event cfg.trialnum};
            
            
            
            Deci.Analysis.Extra.QL.States = [21 23];
            Deci.Analysis.Extra.QL.Actions = [31 32];
            Deci.Analysis.Extra.QL.Reward  = [51 52];
            Deci.Analysis.Extra.QL.Value  = {[10 0] [0 -10]};
            Deci.Analysis.Extra.QL.Start = [0 0];
            
            for Cond = 1:length(Deci.Analysis.Conditions)
                maxt = max(sum(ismember(info.condinfo{2},Deci.Analysis.Conditions{Cond}),2));
                info.alltrials = sum(ismember(info.condinfo{2},Deci.Analysis.Conditions{Cond}),2) == maxt;
                info.allnonnans = ~isnan(mean(info.condinfo{1},2));
                ccfg.trials = info.alltrials & info.allnonnans;
                
                info.subject_list = subject_list;
                info.Cond = Cond;
                
                QL3(Deci,info,info,Deci.Analysis.Extra.QL);
                
            end
        end
%% 4.) Model Analysis 
Deci.Analysis.Conditions    = {[21] [23]};
Deci.Analysis.LocksTitle = {'None'};
info.Lock = 1;
Deci.Analysis.CondTitle = {'Correct' 'Incorrect'}';        


QL = [];
mdl = repmat({[]},[1 2]);
PseudoR = [];
PR = repmat({[]},[1 2]);
for subj = 1:length(Deci.SubjectList)
    for Lock = 1:length(Deci.Analysis.LocksTitle)
        for Cond = 1:length(Deci.Analysis.CondTitle)
            
            load([Deci.Folder.Analysis filesep 'Extra' filesep 'QL' filesep Deci.SubjectList{subj}  filesep filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond}],'QL');
            for QLs = 1:length(QL)
                mdl{QLs}(subj,:,Cond) = QL{QLs};   
            end

            load([Deci.Folder.Analysis filesep 'Extra' filesep 'PseudoR' filesep Deci.SubjectList{subj}  filesep filesep Deci.Analysis.LocksTitle{Lock} filesep Deci.Analysis.CondTitle{Cond}],'PseudoR');
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
        
%         [stat(mdls,param,:),tbl,stats,terms] = anovan(tt(:),{sub(:),cond(:)},'varnames',{'Subject','Condition'},'display','off');
        
        
        excelparam(:,param) = tt(:);
        
    end
    
    exceldata = cat(2,exceldata,excelparam);
    
end

excelmdldata = mat2cell([exceldata sub(:) cond(:)],[length(sub(:))],[ones([1 size(exceldata,2)+2])]);

excelmdldata = table(excelmdldata{:},'VariableNames',{'m1_Q0Opt' 'm1_Q02Wor' 'm1_a' 'm1_b' 'm2_Q0Opt' 'm2_Q0Wor' 'm2_aR' 'm2_aP' 'm2_b' 'm3_a' 'm3_b' 'Sub' 'Cond'});

writetable(excelmdldata,[Deci.Folder.Plot filesep 'Modeloutputs'],'FileType','spreadsheet','Sheet','Model Parameters')

% disp('---------')
% disp('QL model Stats')
% disp('Model (row) x Parameters (col)(4,5,2) x Stat Dimension [Subjects (N), Conditions(2)]');
% stat

% 


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

% 
% disp('---------')
% disp('PseudoR Stat')
% disp('Model (row) x Dimension (col) [Subjects, Conditions]')
% PRstat
% 
% disp('---------')
% disp(['Reward Psuedo-R is ' num2str(mean(PR{PRs}(:,:,1),1)) ' (+-' num2str(std(PR{PRs}(:,:,1),[],1)) ') , Punishment Psuedo-R is ' num2str(mean(PR{PRs}(:,:,2),1)) ' (+-' num2str(std(PR{PRs}(:,:,2),[],1)) ')'])
% disp(['Mean difference is ' num2str(diff(mean(PR{PRs},1)))])

% 

 %% 5.) Behavioral Plots
 
Deci.Analysis.Conditions    = {[21 31 51] [23 31 51] [21 31 52] [23 31 52] ...
                               [21 32 51] [23 32 51] [21 32 52] [23 32 52]};
Deci.Analysis.CondTitle     = {'Opt AB Corr' 'Opt CD Corr' ...
                               'Opt AB Inc'  'Opt CD Inc' ...
                               'Wst AB Corr' 'Wst CD Corr' ...
                               'Wst AB Inc'  'Wst CD Inc' };
 
Deci.Plot.Behv = [];
Deci.Plot.Behv.Source = 'Definition';

Deci.Plot.Behv.Acc.Figure = [true true];
Deci.Plot.Behv.Acc.Total = {{[1 3 5 7] [2 4 6 8]} {[1 3 5 7] [2 4 6 8]}};
Deci.Plot.Behv.Acc.Subtotal = {{[1 3] [2 4]} {[1 5] [2 6]}};
% Change [2 4] to [6 8] to flip punishment
% have to be [2 4] when you run the behavioral stat section (6)


Deci.Plot.Behv.Acc.Subtitle = {{'Rew Percent' 'Pun Percent'} {'Rew Correct' 'Pun Correct'}};

Deci.Plot.Behv.Acc.Block = [1:3];
Deci.Plot.Behv.Acc.Collapse.Trial =false;
% Turn this to true false to turn from wire and bar by collapsing trials

Deci.Plot.Behv.Acc.Collapse.Block = true;
Deci.Plot.Behv.Acc.Collapse.Subject = true;             % turn to false if you want to plot individual subjects
Deci.Plot.Behv.Acc.Collapse.Uneven = 'positional:nans';
Deci.Plot.Behv.Acc.Collapse.Movmean =  [false true];
% Change this to true to activate cummulative scores for Feedback correct
% Change this to false to turn off cummulative for Optimal Choice

Deci.Plot.Behv.Acc.Collapse.MovWindow = 10;
Deci.Plot.Behv.Acc.Title = {'Percent of Optimal Choice' ['Percent of Feedback Correct with sliding window ' num2str(Deci.Plot.Behv.Acc.Collapse.MovWindow)]};
% Change the sliding scale window range. If -1, get complete cummulative

Deci.Plot.Behv.RT.Figure = [false true];
Deci.Plot.Behv.RT.Draw = {{[1 3] [5 7] [2 4] [6 8]} {[1 3 5 7] [2 4 6 8]}};
Deci.Plot.Behv.RT.Title = {'All Trials RT by CondxChoice' 'All Trials RT by Cond'};
Deci.Plot.Behv.RT.Subtitle = {{'Rew-Opt' 'Rew-Wor' 'Pun-Opt' 'Pun-Worst'} {'Rew Percent' 'Pun Percent'}};
Deci.Plot.Behv.RT.Locks = [1 2];

Deci.Plot.Behv.RT.Block = [3];
Deci.Plot.Behv.RT.Collapse.Trial =false;
Deci.Plot.Behv.RT.Collapse.Block =  true;
% Change here for wires to bar, but across Blocks
Deci.Plot.Behv.RT.Collapse.Subject = true;
Deci.Plot.Behv.RT.Collapse.Uneven = 'maxlength:nans';

Plottor_Behv(Deci); 

%% (6) Behv Stats
% have to be [2 4] when you run the behavioral stat 

load([Deci.Folder.Version filesep 'Plot' filesep 'SimAcc'],'fullAcc');

Acc = fullAcc{1};

                subs= [];
                conds = [];
                blks = []; 
                trls = [];
                
             sumsubs = [];
             sumconds= [];
             sumblks= [];
             sumacc = [];
                
for sub = 1:size(Acc,1)
    for cond = 1:size(Acc,2)
        for blk = 1:size(Acc,3)
            for trl = 1:size(Acc,4)
                subs(end+1) =  sub;
                conds(end+1) = cond;
                blks(end+1) = blk;
                trls(end+1) = trl;
            end
            
             sumsubs(end+1) =  sub;
             sumconds(end+1)= cond;
             sumblks(end+1)= blk;
             sumacc(sub,cond,blk) = squeeze(nanmean(Acc(sub,cond,blk,:),4));
        end
    end
end

excelAccdata = table(subs',conds',blks',trls',Acc(:),'VariableNames',{'Subj' 'Cond' 'Blk' 'Trl','Optimal_Choice'});
writetable(excelAccdata,[Deci.Folder.Plot filesep 'Behavioraloutputs'],'FileType','spreadsheet','Sheet','Accuracy');

excelAccdata = table(sumsubs',sumconds',sumblks',sumacc(:),'VariableNames',{'Subj' 'Cond' 'Blk' 'Optimal_Choice'});
writetable(excelAccdata,[Deci.Folder.Plot filesep 'Behavioraloutputs'],'FileType','spreadsheet','Sheet','Accuracy_Summary');

% cond = [];
% for Conds = 1:size(Acc,2)
%     cond{Conds} = cat(1,Acc{:,Conds});
%     condtitle{Conds} = ones(size(cond{Conds}))*Conds;
% end
% [Accstat,tbl,stats,terms] = anovan(cat(1,cond{:}),{cat(1,condtitle{:})},'varnames',{'Condition'},'display','off');
% 
% for Subjs = 1:size(Acc,1)
%    for Conds = 1:size(Acc,2)
%         meanAcc(Subjs,Conds) = mean(Acc{Subjs,Conds});
%         %numelAcc(Subjs,Conds) = length(find(Acc{Subjs,Conds}));
%    end
% end
% disp('-----------')
% disp('Behavioral Data')
% disp(['Choice Accuracy: Reward ' num2str(mean(meanAcc(:,1),1)*100) '(+-' num2str(std(meanAcc(:,1),[],1)*100) ')' ', Punishment ' num2str(mean(meanAcc(:,2),1)*100) '(+-' num2str(std(meanAcc(:,2),[],1)*100) ')']);
% disp(['p value : ' num2str(Accstat)])

load([Deci.Folder.Version filesep 'Plot' filesep 'SimAcc'],'fullAcc');
Acc = fullAcc{2};

                subs= [];
                conds = [];
                blks = []; 
                trls = [];
                
                sumsubs = [];
                sumconds= [];
                sumblks= [];
                sumacc =[];
                
for sub = 1:size(Acc,1)
    for cond = 1:size(Acc,2)
        for blk = 1:size(Acc,3)
            for trl = 1:size(Acc,4)
                subs(end+1) =  sub;
                conds(end+1) = cond;
                blks(end+1) = blk;
                trls(end+1) = trl;
            end
            
            sumsubs(end+1) =  sub;
            sumconds(end+1)= cond;
            sumblks(end+1)= blk;
            sumacc(sub,cond,blk) = squeeze(nanmean(Acc(sub,cond,blk,:),4));
        end
    end
end

excelAccdata = table(subs',conds',blks',trls',Acc(:),'VariableNames',{'Subj' 'Cond' 'Blk' 'Trl','Correct_Feedback'});

writetable(excelAccdata,[Deci.Folder.Plot filesep 'Behavioraloutputs'],'FileType','spreadsheet','Sheet','Correct_Cummulative')

excelAccdata = table(sumsubs',sumconds',sumblks',sumacc(:),'VariableNames',{'Subj' 'Cond' 'Blk' 'Optimal_Choice'});
writetable(excelAccdata,[Deci.Folder.Plot filesep 'Behavioraloutputs'],'FileType','spreadsheet','Sheet','Correct_Cummulative_Summary');


% for Subjs = 1:size(Acc,1)
%    for Conds = 1:size(Acc,2)
%         meanAcc(Subjs,Conds) = mean(Acc{Subjs,Conds});
%         %numelAcc(Subjs,Conds) = length(find(Acc{Subjs,Conds}));
%    end
% end
% disp(['Number Correct: Reward ' num2str(mean(meanAcc(:,1),1)*100) '(+-' num2str(std(meanAcc(:,1),[],1)*100) ')' ', Punishment ' num2str(mean(meanAcc(:,2),1)*100) '(+-' num2str(std(meanAcc(:,2),[],1)*100) ')']);

% load([Deci.Folder.Version filesep 'Plot' filesep 'RT'],'RT');
% 
% for Conds = 1:size(RT,2)
%     cond{Conds} = cat(1,RT{:,Conds});
%     condtitle{Conds} = ones(size(cond{Conds}))*Conds;
% end
% [RTstat,tbl,stats,terms] = anovan(cat(1,cond{:}),{cat(1,condtitle{:})},'varnames',{'Condition'},'display','off');
% 
% 

load([Deci.Folder.Version filesep 'Plot' filesep 'SimRT'],'fullRT');
RT = fullRT{2};

                subs= [];
                conds = [];
                blks = []; 
                trls = [];
                
                sumsubs = [];
                sumconds= [];
                sumblks= [];
                sumRT =[];
                
for sub = 1:size(RT,1)
    for cond = 1:size(RT,2)
        for blk = 1:size(RT,3)
            for trl = 1:size(RT,4)
                subs(end+1) =  sub;
                conds(end+1) = cond;
                blks(end+1) = blk;
                trls(end+1) = trl;
            end
            
            sumsubs(end+1) =  sub;
            sumconds(end+1)= cond;
            sumblks(end+1)= blk;
            sumRT(sub,cond,blk) = squeeze(nanmean(RT(sub,cond,blk,:),4));
        end
    end
end

excelAccdata = table(subs',conds',blks',trls',RT(:),'VariableNames',{'Subj' 'Cond' 'Blk' 'Trl','RT'});

writetable(excelAccdata,[Deci.Folder.Plot filesep 'Behavioraloutputs'],'FileType','spreadsheet','Sheet','RT')

excelAccdata = table(sumsubs',sumconds',sumblks',sumRT(:),'VariableNames',{'Subj' 'Cond' 'Blk' 'Optimal_Choice'});
writetable(excelAccdata,[Deci.Folder.Plot filesep 'Behavioraloutputs'],'FileType','spreadsheet','Sheet','RT_Summary');

   