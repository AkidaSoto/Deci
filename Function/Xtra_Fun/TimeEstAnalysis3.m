clearvars;

includenoRT = false;

% False - Improvements and Adjustment Plots (~20 subj)
% True - Only Adjustment Plots (~80 Subj)

Files = CleanDir('/Users/REINHARTLAB/Desktop/TimeEst/TimeEst_All');

Deci.SubjectList = cellfun(@(c) strsplit(c,'.'),Files,'un',0);
Deci.SubjectList = unique(cellfun(@(c) c{1},Deci.SubjectList,'un',0));

fakeUI = figure;
fakeUI.UserData = Deci.SubjectList;
fakeUI.Visible =  'off';
dc_select_labels(fakeUI,[],Deci.SubjectList);
waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
Deci.SubjectList = fakeUI.UserData;
close(fakeUI);

Block_start = 5;  % new block signal
Block_stop = 6;  % end of block signal
Block_New = 7;
Exp_end = 2; %experiment end


Block_Learned = 33; %denotes whether the preceding block was learned (i.e. when success counter breaks first)
Block_UnLearned = 34; %if preceding blocked was not learned (i.e. when error counter breaks first)

Trial_start = 10; %start of each trial  (600 trials total). essentially the "i" variable
Trial_end = 11; %end of each trial


Resp_early = 19; %early response
Resp_late = 20;    %late response
Resp_correct = 21; %correct response

Stim_Onset = 12; %stim onset - cue box
Resp_all  = 17; % General Response Marker (when a participant responds)
Resp_none = 18; % no response

NewSubs = [];

if ~includenoRT

for Sub = 1:length(Files)
    
    load(['/Users/REINHARTLAB/Desktop/TimeEst/TimeEst_All' filesep Files{Sub}]);
    
    if ~includenoRT
    if ~exist('RT_lower')
        continue;
    end
    end
    
    events =  StandardizeEventMarkers(events);
    
    Block_beg = find(ismember({events.value},num2str(Block_start)));
    Block_end = Block_beg(2:end);
    Block_end = [Block_end find(ismember({events.value},num2str(Exp_end)))];
    
    if length(RT_lower) ~= length(Block_beg)
        continue;
    end
    
    NewSubs{end+1} = Files{Sub};
    
end
else
   NewSubs = Files;
end
%% Run Here

   SubDeviationAllL = [];
    SubDeviationAllU = [];

for Sub = 1:length(NewSubs)
    
    Learnness = [];
    BlockLength = [];
    
    load(['/Users/REINHARTLAB/Desktop/TimeEst/TimeEst_All' filesep NewSubs{Sub}]);
    
    events =  StandardizeEventMarkers(events);
    
    Block_beg = find(ismember({events.value},num2str(Block_start)));
    Block_end = Block_beg(2:end);
    Block_end = [Block_end find(ismember({events.value},num2str(Exp_end)))];
    
    AllRT = [];
    AllTT = [];
    AllCorrectness = [];
    MeanofAbsDeviation = [];
    AbsofMeanDeviation = [];
    
    AdjustmentLearned = repmat({[]},[1 2]);
    AdjustmentCLearned = repmat({[]},[1 2]);
    TrialCount = repmat({[]},[1 2]);
    DeviationAll = [];
 
    for Type = 1:length(Block_beg)-1
        
        Block = events(Block_beg(Type):Block_end(Type));
        
        Values = cellfun(@str2num,{Block.value});
        Times = [Block.sample];
        Learn = any(find(ismember(Values,Block_Learned)));
        UnLearn = any(find(ismember(Values,Block_UnLearned)));
        
        if Learn  == UnLearn
            error('does this even happen?')
        end
        Learnness(Type) = Learn;
        
        
        
        Trial = [];
        Trial(:,1) =  find(ismember(Values,Trial_start));
        Trial(:,2) =  find(ismember(Values,Trial_end));
        
        TrialCount{UnLearn+1}(end+1) = size(Trial,1);
        
        BlockLength{Sub}(Type) = size(Trial,1);
        
        Correctness = [];
        Deviation = [];
        wDeviation = [];
        RT = [];
        Adjustment = [];
        Improvement =[];
        AdjustmentC = [];
        AdjustmentT = [];
        
        for Trl = 1:size(Trial,1)
            
            TrialValues =  Values(Trial(Trl,1):Trial(Trl,2));
            TrialTimes = Times(Trial(Trl,1):Trial(Trl,2));
            RspAll = ismember(Resp_all,TrialValues);
            
            if RspAll
                
                % Find what type of RspAll trial it was
                RspCor = ismember(Resp_correct, Values(Trial(Trl,1):Trial(Trl,2)));
                RspErr = any(ismember([Resp_early Resp_late], Values(Trial(Trl,1):Trial(Trl,2))));
       
                if RspCor == RspErr
                    error('does this even happen?')
                end
                
                Correctness(end+1) = double(RspCor);
                AllCorrectness(end+1) = double(RspCor);
                
                %Find Stim Onset Time
                StimTime = TrialTimes(ismember(TrialValues,Stim_Onset));
                
                %Find Rsp Onset Time
                RspTime = TrialTimes(ismember(TrialValues,Resp_all));
               
                
                RT(end+1) = [RspTime - StimTime];
                AllRT(Sub,end+1) =  [RspTime - StimTime];

                
                
                if ~includenoRT
                    Deviation(end+1) = RT(end) - Target_Time(Type)*1000;
                    AllTT(Sub,end+1) = Target_Time(Type)*1000;
                    if RT(end) < RT_lower(Type)*1000 % early
                        wDeviation(end+1) = RT(end) - RT_lower(Type)*1000;
                        
                        if RspCor && diff([RT(end) RT_lower(Type)*1000]) > 1.01
                            error('does this even happen?')
                        end
                        
                    elseif RT(end) > RT_upper(Type)*1000 % late
                        wDeviation(end+1) = RT(end) - RT_upper(Type)*1000;
                        
                        if RspCor && diff([RT(end) RT_lower(Type)*1000]) > 1.01
                            error('does this even happen?')
                        end
                        
                    else
                        wDeviation(end+1) = 0;
                        
                        if RspErr && diff([RT(end) RT_lower(Type)*1000]) > 1.01
                            error('does this even happen?')
                        end
                    end
                end
                
                if length(RT) > 1
                    Adjust = abs(RT(end) - RT(end-1));
                    
                    if diff(Correctness(end-1:end)) == 1 % Error to Correct
                       Adjustment{1}(end+1) = Adjust;
                       AdjustmentC{1}(end+1) = 1;
                       
                       if length(RT) > 2 
                           if diff(Correctness(end-2:end-1)) == -1
                                AdjustmentT{1}(end+1) =  abs(RT(end) - RT(end-2));   
                           end
                       end
                       
                    elseif diff(Correctness(end-1:end)) == -1 %Correct to Error
                       Adjustment{2}(end+1) = Adjust;
                       AdjustmentC{2}(end+1) = 1;

                    elseif diff(Correctness(end-1:end)) == 0 && Correctness(end) %Correct to Correct
                       Adjustment{3}(end+1) = Adjust;
                       AdjustmentC{3}(end+1) = 1;
                       
                       if length(RT) > 2
                           if diff(Correctness(end-2:end-1)) == 0 && Correctness(end-1)
                               AdjustmentT{2}(end+1)  =   abs(RT(end) - RT(end-2));   
                           end
                       end
                    elseif diff(Correctness(end-1:end)) == 0 && ~Correctness(end) %Error to Error
                       Adjustment{4}(end+1) = Adjust;
                       AdjustmentC{4}(end+1) = 1;
                    end
                    
                    if ~includenoRT
                        Improve = abs(wDeviation(end)) - abs(wDeviation(end-1));
                        
                        if diff(Correctness(end-1:end)) == 1 % Error to Correct
                            Improvement{1}(end+1) = Improve;
                        elseif diff(Correctness(end-1:end)) == -1 %Correct to Error
                            Improvement{2}(end+1) = Improve;
                        elseif diff(Correctness(end-1:end)) == 0 && Correctness(end) %Correct to Correct
                            Improvement{3}(end+1) = Improve;
                        elseif diff(Correctness(end-1:end)) == 0 && ~Correctness(end) %Error to Error
                            Improvement{4}(end+1) =Improve;
                        end
                    end
                    
                else
                    AdjustmentT{1} = [];
                    AdjustmentT{2} = [];
                    
                    Adjustment{1} = [];
                    Adjustment{2} = [];
                    Adjustment{3} = [];
                    Adjustment{4} = [];
                    
                    AdjustmentC{1} = [];
                    AdjustmentC{2} = [];
                    AdjustmentC{3} = [];
                    AdjustmentC{4} = [];
                    
                    if ~includenoRT
                    Improvement{1} = [];
                    Improvement{2} = [];
                    Improvement{3} = [];
                    Improvement{4} = [];
                    end
                end
                
            else
                continue;
                
            end
            
            
        end
        
        RLT = rem(length(Deviation),4);
        
        AllDeviation(Type) = mean(abs(wDeviation),2);
        
        
        if size(DeviationAll,2) >= length(Deviation)
            DeviationAll(Type,:) = nan([1 size(DeviationAll,2)]);
            DeviationAll(Type,1:length(Deviation))  = Deviation;
        else
            DeviationAll(:,end+1:end+[length(Deviation)-size(DeviationAll,2)]) = nan;
            DeviationAll(Type,:) = Deviation;
        end
        
        if RLT~=0
        Deviation = reshape([Deviation nan(1,4-RLT)],4,ceil(length(Deviation)/4));
        else
        Deviation = reshape([Deviation],4,ceil(length(Deviation)/4));
        end

        MeanofAbsDeviation(Type,:) = nanmean(abs(Deviation),2);   
        AbsofMeanDeviation(Type,:) = nanmean(Deviation,2);
        AdjustmentBlock(Type,:) = cellfun(@mean,Adjustment);
        
        if ~isempty(AdjustmentT)
            AdjustmentBlockT(Type,:) = cellfun(@mean,AdjustmentT);
        else
           AdjustmentBlockT(Type,:) = [nan nan];
        end
        
        AdjustmentLearned{UnLearn+1}(end+1,:) = AdjustmentBlock(Type,:);
        AdjustmentCLearned{UnLearn+1}(end+1,:) =cellfun(@length,Adjustment);
        
        if ~includenoRT
        RLT = rem(length(wDeviation),4);
        if RLT~=0
        wDeviation = reshape([wDeviation nan(1,4-RLT)],4,ceil(length(wDeviation)/4));
        else
        wDeviation = reshape([wDeviation],4,ceil(length(wDeviation)/4));
        end

        MeanofAbswDeviation(Type,:) = nanmean(abs(wDeviation),2);
        AbsofMeanwDeviation(Type,:) = nanmean(wDeviation,2);

        
        if any(any(isnan( MeanofAbsDeviation)))
           k = 0; 
        end
        ImprovementBlock(Type,:) = cellfun(@mean,Improvement);
        end

    end
    
    DeviationAll = abs(DeviationAll);
    if size(SubDeviationAllL,2) >= length(nanmean(DeviationAll(logical(Learnness),:),1))
        SubDeviationAllL(Sub,:) = nan([1 size(SubDeviationAllL,2)]);
        SubDeviationAllL(Sub,1:length(nanmean(DeviationAll(logical(Learnness),:),1)))  = nanmean(DeviationAll(logical(Learnness),:),1);
    else
        SubDeviationAllL(:,end+1:end+[length(nanmean(DeviationAll(logical(Learnness),:),1))-size(SubDeviationAllL,2)]) = nan;
        SubDeviationAllL(Sub,:) = nanmean(DeviationAll(logical(Learnness),:),1);
    end
    
    if size(SubDeviationAllU,2) >= length(nanmean(DeviationAll(~logical(Learnness),:),1))
        SubDeviationAllU(Sub,:) = nan([1 size(SubDeviationAllU,2)]);
        SubDeviationAllU(Sub,1:length(nanmean(DeviationAll(~logical(Learnness),:),1)))  = nanmean(DeviationAll(~logical(Learnness),:),1);
    else
        SubDeviationAllU(:,end+1:end+[length(nanmean(DeviationAll(~logical(Learnness),:),1))-size(SubDeviationAllU,2)]) = nan;
        SubDeviationAllU(Sub,:) = nanmean(DeviationAll(~logical(Learnness),:),1);
    end
    
    
    SubTrialCount(Sub,:) = [mean(TrialCount{1}) mean(TrialCount{2})];
    
    SubAdjustment(Sub,:) = nanmean(AdjustmentBlock,1);
    SubAdjustmentL(Sub,:) = nanmean(AdjustmentLearned{1},1);
    SubAdjustmentU(Sub,:) = nanmean(AdjustmentLearned{2},1);
    
    SubAdjustmentT(Sub,:) = nanmean(AdjustmentBlockT,1);
    
    SubAdjustmentCL(Sub,:) = nanmean(AdjustmentCLearned{1},1);
    SubAdjustmentCU(Sub,:) = nanmean(AdjustmentCLearned{2},1);
    
    LearnInfo(Sub,:,1) =  nanmean(MeanofAbsDeviation(logical(Learnness),:),1);
    LearnInfo(Sub,:,2) =  abs(nanmean(AbsofMeanDeviation(logical(Learnness),:),1));
     
        
    AllLearnInfo(Sub,1) =  nanmean(AllDeviation(logical(Learnness)));
    AllUnLearnInfo(Sub,1) =  nanmean(AllDeviation(logical(~Learnness)));

    UnLearnInfo(Sub,:,1) =  nanmean(MeanofAbsDeviation(logical(~Learnness),:),1);
    UnLearnInfo(Sub,:,2) =  abs(nanmean(AbsofMeanDeviation(logical(~Learnness),:),1));
    
    if ~includenoRT
    SubImprovement(Sub,:) =nanmean(ImprovementBlock,1);
        
    wLearnInfo(Sub,:,1) =  nanmean(MeanofAbswDeviation(logical(Learnness),:),1);
    wLearnInfo(Sub,:,2) =  abs(nanmean(AbsofMeanwDeviation(logical(Learnness),:),1));
     
    wUnLearnInfo(Sub,:,1) =  nanmean(MeanofAbswDeviation(logical(~Learnness),:),1);
    wUnLearnInfo(Sub,:,2) =  abs(nanmean(AbsofMeanwDeviation(logical(~Learnness),:),1));
    end
%     figure;
%     plot(1:length(AllRT),AllRT);
%     hold on
%     plot(1:length(AllTT),AllTT);
    disp(['--------------------'])
    disp(Files{Sub})
    disp(['Percent of Blocks that are Learned: '  num2str(mean(Learnness)*100) '%']);
    disp(['Percent of Trials Correct: '  num2str(mean(AllCorrectness)*100) '%']);
end

figure;
CleanBars(mean(SubTrialCount),sem(SubTrialCount))
legend({'Learned','Not Learned'})
title('# of Trials to learn block type')

figure;
plot(nanmean(SubDeviationAllL,1))
hold on
plot(nanmean(SubDeviationAllU,1))
legend({'Learned','Not Learned'})
title('All Trials Learned')


figure;
plot(mean(LearnInfo(:,:,1),1))
hold on
plot(mean(UnLearnInfo(:,:,1),1))
legend({'Learned','Not Learned'})
title('Mean of Abs Deviation')

figure;
plot(mean(LearnInfo(:,:,2),1))
hold on
plot(mean(UnLearnInfo(:,:,2),1))
legend({'Learned','Not Learned'})
title('Abs of Mean Deviation')

figure;
CleanBars([mean(AllLearnInfo(:,:,1),1) mean(AllUnLearnInfo(:,:,1),1)])
legend({'Learned','Not Learned'})
title('Mean of Abs Deviation')



if ~includenoRT
figure;
plot(mean(wLearnInfo(:,:,1),1))
hold on
plot(mean(wUnLearnInfo(:,:,1),1))
legend({'Learned','Not Learned'})
title('Mean of Abs Deviation from Window')

figure;
plot(mean(wLearnInfo(:,:,2),1))
hold on
plot(mean(wUnLearnInfo(:,:,2),1))
legend({'Learned','Not Learned'})
title('Abs of Mean Deviation from Window')
end

figure
CleanBars(nanmean(SubAdjustment,1))
title('Doublet (Adjustment in ms)')
legend({'E-C','C-E','C-C','E-E'})

figure
CleanBars([nanmean(SubAdjustmentL,1); nanmean(SubAdjustmentU,1)]')
title('Doublet (Adjustment in ms)')
xticklabels({'E-C','C-E','C-C','E-E'})
legend({'Learned' 'Unlearned'})

figure
CleanBars([nanmean(SubAdjustmentCL,1); nanmean(SubAdjustmentCU,1)]')
title('Doublet (Count in ms)')
xticklabels({'E-C','C-E','C-C','E-E'})
legend({'Learned' 'Unlearned'})

if ~includenoRT
figure
CleanBars(nanmean(SubImprovement,1))
title('Doublet (Improvement in ms)')
legend({'E-C','C-E','C-C','E-E'})
end

figure
CleanBars(nanmean(SubAdjustmentT,1),sem(SubAdjustmentT))
title('Triplet (Adjustment in ms)')
legend({'C-E-C','C-C-C'})
