clearvars;
Files = CleanDir('C:\Users\User\Desktop\Radhika\TimeEst_Raw');

Stim = cellfun(@(c) ~rem(str2num(c(end-6)),2),Files);


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

Resp_early = 19; %early response
Resp_late = 20;    %late response
Resp_correct = 21; %correct response

UTrialDev = nan([length(Files) 100 120]);
LTrialDev = nan([length(Files) 100 120]);
UTrialDevtt = nan([length(Files) 100 120]);
LTrialDevtt = nan([length(Files) 100 120]);
LBlockDev = nan([length(Files) 100]);
ULBlockDev = nan([length(Files) 100]);
LBlockwDev = nan([length(Files) 100]);
ULBlockwDev = nan([length(Files) 100]);
LPair = [repmat({[]},[length(Files),9])];
ULPair = [repmat({[]},[length(Files),9])];
TotalPair = [repmat({[]},[length(Files),9])];

LBinsFHalf =[];
LBinsLHalf =[];
UBinsFHalf =[];
UBinsLHalf =[];
 
LPairFHalf =[];
LPairLHalf =[];
UPairFHalf =[];
UPairLHalf =[];
%%

for Sub = 1:length(Files)
    
    load(['C:\Users\User\Desktop\Radhika\TimeEst_Raw' filesep Files{Sub}]);
    
    if ~exist('RT_lower')
        continue;
    end
    
%     
    LBinsFHalfCOUNTER = 1;
    LBinsLHalfCOUNTER = 1;
    UBinsFHalfCOUNTER = 1;
    UBinsLHalfCOUNTER = 1;
    
    LPairFHalfCOUNTER =1;
    LPairLHalfCOUNTER =1;
    UPairFHalfCOUNTER =1;
    UPairLHalfCOUNTER =1;
%     
    LBNum = 1;
    ULBNum = 1;
    DevPlot = [];
    TTPlot = [];
    
    events =  StandardizeEventMarkers(events);
    
    Block_beg = find(ismember({events.value},num2str(Block_start)));
    Block_end = Block_beg(2:end);
    Block_end = [Block_end find(ismember({events.value},num2str(Exp_end)))];
    
    for Type = 1:length(Block_beg)
        
        Exam = {events(Block_beg(Type):Block_end(Type)).value};
        Time = {events(Block_beg(Type):Block_end(Type)).sample};
        Learn = any(find(ismember({events(Block_beg(Type):Block_end(Type)).value},num2str(Block_Learned))));
        UnLearn = any(find(ismember({events(Block_beg(Type):Block_end(Type)).value},num2str(Block_UnLearned))));
        
        totalblocks(Sub) = size(Block_end,2);
        avgblocks = mean(totalblocks);
        sdblocks = std(totalblocks);
        
        if ~Learn&& ~UnLearn && ~ismember(num2str(Exp_end),{events(Block_beg(Type):Block_end(Type)).value})
            error('No Block Markers');
        end
        
        clear Trial
        
        Trial(:,1) =  find(ismember(Exam,num2str(Trial_start)));
        Trial(:,2) =  find(ismember(Exam,num2str(Trial_end)));
        
        Correct = [];
        wDeviation = [];
        Deviation = [];
        CTL = [];
        
        for Trl = 1:size(Trial,1)
            
            RspAll = ismember(num2str(Resp_all), Exam(Trial(Trl,1):Trial(Trl,2)));
            
            if RspAll
                
                % Find what type of RspAll trial it was
                RspCor = ismember(num2str(Resp_correct), Exam(Trial(Trl,1):Trial(Trl,2)));
                RspLat = ismember(num2str(Resp_early), Exam(Trial(Trl,1):Trial(Trl,2)));
                RspEar = ismember(num2str(Resp_late), Exam(Trial(Trl,1):Trial(Trl,2)));
                
                
                TrialTimes = Time(Trial(Trl,1):Trial(Trl,2));
                
                %Find Stim Onset Time
                StimTime = TrialTimes{ismember(Exam(Trial(Trl,1):Trial(Trl,2)),num2str(Stim_Onset))};
                
                %Find Rsp Onset Time
                RspTime = TrialTimes{ismember(Exam(Trial(Trl,1):Trial(Trl,2)),num2str(Resp_all))};
                
                RT = [RspTime - StimTime];
                
                if RspEar %% && [RspTime - StimTime]/1000 <= RT_lower(Type)
                    Correct(Trl) = 2;
                    wDeviation(Trl) = RT - RT_lower(Type)*1000;
                    Deviation(Trl) = RT - Target_Time(Type)*1000;
                    
                elseif RspLat %% && [RspTime - StimTime]/1000 >= RT_upper(Type)
                    Correct(Trl) = 3;
                    wDeviation(Trl) = RT - RT_upper(Type)*1000;
                    Deviation(Trl) = RT - Target_Time(Type)*1000;
                    
                elseif RspCor %% && [RspTime - StimTime]/1000 >= RT_lower(Type) && [RspTime - StimTime]/1000 <= RT_upper(Type)
                    Correct(Trl) = 1;
                    wDeviation(Trl) = 0;
                    Deviation(Trl) = RT - Target_Time(Type)*1000;
                    
                end
                
                if Trl ~= 1
                    
                    if ~isnan(Correct(Trl-1))
                        
                        LRspTime = Time(Trial(Trl-1,1):Trial(Trl-1,2));
                        LRspTime = LRspTime{ismember(Exam(Trial(Trl-1,1):Trial(Trl-1,2)),num2str(Resp_all))};
                        
                        PMath = [Correct(Trl-1)-1]*3 + Correct(Trl);
                        
                        if Learn
                            LPair{Sub,PMath}(end+1) = Deviation(Trl) - Deviation(Trl-1);
                            
                            if Type < ceil(length(Block_beg)/2)
                                LPairFHalf{Sub,PMath}(LPairFHalfCOUNTER) =  LPair{Sub,PMath}(end);
                                LPairFHalfCOUNTER =LPairFHalfCOUNTER +1;
                            else
                                LPairLHalf{Sub,PMath}(LPairLHalfCOUNTER) = LPair{Sub,PMath}(end);
                                LPairLHalfCOUNTER =LPairLHalfCOUNTER +1;
                            end
                            
                        elseif  UnLearn
                            ULPair{Sub,PMath}(end+1) = Deviation(Trl) - Deviation(Trl-1);
                            
                            if Type < ceil(length(Block_beg)/2)
                                UPairFHalf{Sub,PMath}(UPairFHalfCOUNTER) = ULPair{Sub,PMath}(end);
                                UPairFHalfCOUNTER =UPairFHalfCOUNTER +1;
                            else
                                UPairLHalf{Sub,PMath}(UPairLHalfCOUNTER) = ULPair{Sub,PMath}(end);
                                UPairLHalfCOUNTER =UPairLHalfCOUNTER +1;
                            end
                            
                        end
                        
                        TotalPair{Sub,PMath}(end+1) = Deviation(Trl) - Deviation(Trl-1);
                        
                    end
                    
                end
                
            else
                Correct(Trl) = nan;
                wDeviation(Trl) = nan;
                Deviation(Trl) = nan;
            end
            
            DevPlot(end+1) = [RspTime - StimTime]/1000;
            TTPlot(end+1) =  Target_Time(Type);
            
            
        end
        
        %         BlockDev(Sub,Type) = nanmean(abs(Deviation));
        %         BlockwDev(Sub,Type) = nanmean(abs(wDeviation));
        
        if Learn || UnLearn
            
            RLT = rem(length(Deviation),floor(length(Deviation)/4));
            
            Dev = Deviation(1:end-RLT);
            wDev = wDeviation(1:end-RLT);
            
            
            if Learn
                
                LBlockDev(Sub,LBNum) = nanmean(abs(Deviation));
                LBlockwDev(Sub,LBNum) = nanmean(abs(wDeviation));
                
                LTrialDevtt(Sub,end+1,1:length(wDeviation)) = Deviation;
                LTrialDev(Sub,end+1,1:length(wDeviation)) = wDeviation;
                
                clear LearnedBin
                clear wLearnedBin
                
                v = 1;
                for k = 1:4
                    y = k > [ceil([length(Deviation)+.01]/4)*4 - length(Deviation)];
                    
                    LearnedBin(k,1:ceil(length(Deviation)/4)) = nan;
                    LearnedBin(k,1:[floor(length(Deviation)/4)]+y) = Deviation(v:v+y+[floor(length(Deviation)/4)]-1);
                    
                    wLearnedBin(k,1:ceil(length(wDeviation)/4)) = nan;
                    wLearnedBin(k,1:[floor(length(wDeviation)/4)]+y) = wDeviation(v:v+y+[floor(length(wDeviation)/4)]-1);
                    
                    v = v + y + [floor(length(Deviation)/4)];
                    
                end
                
                LBins{Sub}(LBNum,:,:) = [abs(nanmean(LearnedBin,2)) abs(nanmean(wLearnedBin,2)) nanmean(abs(LearnedBin),2)  nanmean(abs(wLearnedBin),2)];
                
                LBNum = LBNum+1;
                
                
                else
                
                ULBlockDev(Sub,ULBNum) = nanmean(abs(Deviation));
                ULBlockwDev(Sub,ULBNum) = nanmean(abs(wDeviation));
                
                UTrialDevtt(Sub,end+1,1:length(wDeviation)) = Deviation;
                UTrialDev(Sub,end+1,1:length(wDeviation)) = wDeviation;
                
                clear ULearnedBin
                clear wULearnedBin
                
                v = 1;
                for k = 1:4
                    y = k > [ceil([length(Deviation)+.01]/4)*4 - length(Deviation)];
                    
                    ULearnedBin(k,1:ceil(length(Deviation)/4)) = nan;
                    ULearnedBin(k,1:[floor(length(Deviation)/4)]+y) = Deviation(v:v+y+[floor(length(Deviation)/4)]-1);
                    
                    wULearnedBin(k,1:ceil(length(wDeviation)/4)) = nan;
                    wULearnedBin(k,1:[floor(length(wDeviation)/4)]+y) = wDeviation(v:v+y+[floor(length(wDeviation)/4)]-1);
                    
                    v = v + y + [floor(length(Deviation)/4)];
                    
                end
                
                
                UBins{Sub}(ULBNum,:,:) = [abs(nanmean(ULearnedBin,2)) abs(nanmean(wULearnedBin,2)) nanmean(abs(ULearnedBin),2)  nanmean(abs(wULearnedBin),2)];
                ULBNum = ULBNum +1;
            end
            
        end
        
        if Learn
            if Type < ceil(length(Block_beg)/2)
                LBinsFHalf{Sub}(LBinsFHalfCOUNTER,:,:) = squeeze(mean(LBins{Sub}(1:ceil(size(LBins{Sub},1)/2),:,:),1));
                LBinsFHalfCOUNTER = LBinsFHalfCOUNTER + 1;
                
                
            else
                LBinsLHalf{Sub}(LBinsLHalfCOUNTER,:,:) = squeeze(mean(LBins{Sub}(ceil(size(LBins{Sub},1)/2)+1:size(LBins{Sub},1),:,:),1));
                LBinsLHalfCOUNTER = LBinsLHalfCOUNTER + 1;
                
            end
            
            
            
        elseif UnLearn
            if Type < ceil(length(Block_beg)/2)
                 UBinsFHalf{Sub}(UBinsFHalfCOUNTER,:,:) = squeeze(mean(UBins{Sub}(1:ceil(size(UBins{Sub},1)/2),:,:),1));
                 UBinsFHalfCOUNTER = UBinsFHalfCOUNTER + 1;
            else
                 UBinsLHalf{Sub}(UBinsLHalfCOUNTER,:,:) = squeeze(mean(UBins{Sub}(ceil(size(UBins{Sub},1)/2)+1:size(UBins{Sub},1),:,:),1));
                 UBinsLHalfCOUNTER = UBinsLHalfCOUNTER + 1;
            end
        end
        
    end
 
     tlblocks(Sub) = LBNum -1;
     avglblocks = mean(tlblocks);
     sdlblocks = std(tlblocks);
        
    TotalLearned = length(find(ismember({events.value},num2str(Block_Learned))))/length(find(ismember({events.value},num2str(Block_New))));
    
    nanmean(TTPlot - DevPlot);
    disp(['--------------------'])
    disp(Files{Sub})
    disp(['Percent of Blocks that are Learned: '  num2str(TotalLearned*100) '%']);
    disp(['Average abs Deviation from Target Time '  num2str(nanmean(abs(TTPlot - DevPlot))) 'ms']);
    disp(['Average % Deviation from Target Time '  num2str(nanmean([DevPlot./TTPlot]-1)*100) '%']);
    disp(['--------------------'])
    
end
%%

totalstblocks = totalblocks(Stim);
avgstb = mean(totalstblocks);
sdstblocks = std(totalstblocks);

totalshblocks = totalblocks(~Stim);
avgshb = mean(totalshblocks);
sdshblocks = std(totalshblocks);

tLstblocks = tlblocks(Stim);       
avglstb = mean(tLstblocks);
sdlstblocks = std(tLstblocks);

tLshblocks = tlblocks(~Stim);
avglshb = mean(tLshblocks);
sdlshblocks = std(tLshblocks);

%Find size of LearnedPair and UnlearnedPair for S.E.M.
% learnedsize1 = sum(cellfun(@numel,LPair(Stim,:)),1); 
% learnedsize2 = sum(cellfun(@numel,LPair(~Stim,:)),1);
% learnedsize3 = sum(cellfun(@numel,ULPair(Stim,:)),1);
% learnedsize4 = sum(cellfun(@numel,ULPair(~Stim,:)),1);

%Find mean, SD, SEM of LearnedPair 
LearnedPair(:,1) = nanmean(cellfun(@nanmean,LPair(Stim,:)),1);
LearnedPair(:,2)= nanmean(cellfun(@nanmean,LPair(~Stim,:)),1);
sdLearnedPair(:,1) = nanstd(cellfun(@nanmean,LPair(Stim,:)),1);
sdLearnedPair(:,2) = nanstd(cellfun(@nanmean,LPair(~Stim,:)),1);
semLearnedPair(:,1) = sdLearnedPair(:,1)/sqrt(13);
semLearnedPair(:,2) = sdLearnedPair(:,2)/sqrt(13);
%Find mean, SD, SEM of UnlearnedPair 
UnlearnedPair(:,1) = nanmean(cellfun(@nanmean,ULPair(Stim,:)),1);
UnlearnedPair(:,2) = nanmean(cellfun(@nanmean,ULPair(~Stim,:)),1);
sdUnlearnedPair(:,1) = nanstd(cellfun(@nanmean,ULPair(Stim,:)),1);
sdUnlearnedPair(:,2) = nanstd(cellfun(@nanmean,ULPair(~Stim,:)),1);
semUnlearnedPair(:,1) = sdUnlearnedPair(:,1)/sqrt(13);
semUnlearnedPair(:,2) = sdUnlearnedPair(:,2)/sqrt(13);

FHLearnedPair(:,1) = nanmean(cellfun(@nanmean,LPairFHalf(Stim,:)),1);
FHLearnedPair(:,2) = nanmean(cellfun(@nanmean,LPairFHalf(~Stim,:)),1);
LHLearnedPair(:,1) = nanmean(cellfun(@nanmean,LPairLHalf(Stim,:)),1);
LHLearnedPair(:,2) = nanmean(cellfun(@nanmean,LPairLHalf(~Stim,:)),1);
sdFHLearnedPair(:,1) = nanstd(cellfun(@nanmean,LPairFHalf(Stim,:)),1);
sdFHLearnedPair(:,2) = nanstd(cellfun(@nanmean,LPairFHalf(~Stim,:)),1);
sdLHLearnedPair(:,1) = nanstd(cellfun(@nanmean,LPairLHalf(Stim,:)),1);
sdLHLearnedPair(:,2) = nanstd(cellfun(@nanmean,LPairLHalf(~Stim,:)),1);
semFHLearnedPair(:,1) = sdFHLearnedPair(:,1)/sqrt(13);
semFHLearnedPair(:,2) = sdFHLearnedPair(:,2)/sqrt(13);
semLHLearnedPair(:,1) = sdLHLearnedPair(:,1)/sqrt(13);
semLHLearnedPair(:,2) = sdLHLearnedPair(:,2)/sqrt(13);


FHUnlearnedPair(:,1) = nanmean(cellfun(@nanmean,UPairFHalf(Stim,:)),1);
FHUnlearnedPair(:,2) = nanmean(cellfun(@nanmean,UPairFHalf(~Stim,:)),1);
LHUnlearnedPair(:,1) = nanmean(cellfun(@nanmean,UPairLHalf(Stim,:)),1);
LHUnlearnedPair(:,2) = nanmean(cellfun(@nanmean,UPairLHalf(~Stim,:)),1);
sdFHUnlearnedPair(:,1) = nanstd(cellfun(@nanmean,UPairFHalf(Stim,:)),1);
sdFHUnlearnedPair(:,2) = nanstd(cellfun(@nanmean,UPairFHalf(~Stim,:)),1);
sdLHUnlearnedPair(:,1) = nanstd(cellfun(@nanmean,UPairLHalf(Stim,:)),1);
sdLHUnlearnedPair(:,2) = nanstd(cellfun(@nanmean,UPairLHalf(~Stim,:)),1);
semFHUnlearnedPair(:,1) = sdFHUnlearnedPair(:,1)/sqrt(13);
semFHUnlearnedPair(:,2) = sdFHUnlearnedPair(:,2)/sqrt(13);
semLHUnlearnedPair(:,1) = sdLHUnlearnedPair(:,1)/sqrt(13);
semLHUnlearnedPair(:,2) = sdLHUnlearnedPair(:,2)/sqrt(13);
% A = [];
% X = [];
%  for j = [1]
%      A(1,:) = cellfun(@nanmean,LPairFHalf(Stim,j));
%      A(2,:) = cellfun(@nanmean,LPairLHalf(Stim,j));
%      A(1,:) = cellfun(@nanmean,UPairFHalf(Stim,j));
%      A(2,:) = cellfun(@nanmean,UPairFHalf(~Stim,j));
%      A'
%  end
% % 
%  for x = [1]
%      X(1,:) = cellfun(@nanmean,LPairFHalf(~Stim,x));
%      X(2,:) = cellfun(@nanmean,LPairLHalf(~Stim,x));
%      X(1,:) = cellfun(@nanmean,UPairLHalf(Stim,x));
%      X(2,:) = cellfun(@nanmean,UPairLHalf(~Stim,x));
%      X'
%  end

titles = {'Adjustment in C-C Pairs' 'Adjustment in C-E Pairs' 'Adjustment in C-L Pairs' 'Adjustment in E-C Pairs' 'Adjustment in E-E Pairs' 'Adjustment in E-L Pairs' 'Adjustment in L-C Pairs' 'Adjustment in L-E Pairs' 'Adjustment in L-L Pairs'};
for j = 1:size(FHLearnedPair)
   h = figure;
   hold on
   
   CurrentPair = [FHLearnedPair(j,:) FHUnlearnedPair(j,:); LHLearnedPair(j,:) LHUnlearnedPair(j,:)];
   CurrentSEM = [semFHLearnedPair(j,:) semLHLearnedPair(j,:) ; semFHUnlearnedPair(j,:) semLHUnlearnedPair(j,:)];
   bar(CurrentPair)
   ngroups = size(CurrentPair,1);
   nbars = size(CurrentPair,2);
   groupwidth = min(0.8, nbars/(nbars+1.5));
for r = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*r-1) * groupwidth/(2*nbars);
    errorbar(x, CurrentPair(:,r), CurrentSEM(:,r), '.k');
end
   xticks([1 2])
   xticklabels({'Learned' 'Unlearned'})
   legend({'Stim FH' 'Sham FH' 'Stim LH' 'Sham LH'})
   title(titles(j));
end


TotalPairs(:,1) = nanmean(cellfun(@nanmean,TotalPair),1);

%Size of LNum and UNum is the same as LearnedPair/UnlearnedPair for S.E.M.
%Find mean, SD, SEM of LNum 
LCor(:,1) = nanmean(cellfun(@numel,LPair(Stim,1)),1);
LCor(:,2) = nanmean(cellfun(@numel,LPair(~Stim,1)),1);
LErr(:,1) = (nanmean(cellfun(@numel,LPair(Stim,4)),1)) + (nanmean(cellfun(@numel,LPair(Stim,7)),1));
LErr(:,2) = (nanmean(cellfun(@numel,LPair(~Stim,4)),1)) + (nanmean(cellfun(@numel,LPair(~Stim,7)),1));
sdLCor(:,1) = nanstd(cellfun(@numel,LPair(Stim,1)),1);
sdLCor(:,2) = nanstd(cellfun(@numel,LPair(~Stim,1)),1);
sdLErr(:,1) = (nanstd(cellfun(@numel,LPair(Stim,4)),1)) + (nanstd(cellfun(@numel,LPair(Stim,7)),1));
sdLErr(:,2) = (nanstd(cellfun(@numel,LPair(~Stim,4)),1)) + (nanstd(cellfun(@numel,LPair(~Stim,7)),1));
semLCor(:,1) = sdLCor(:,1)/sqrt(13);
semLCor(:,2) = sdLCor(:,2)/sqrt(13);
semLErr(:,1) = sdLErr(:,1)/sqrt(13);
semLErr(:,2) = sdLErr(:,2)/sqrt(13);

UCor(:,1) = nanmean(cellfun(@numel,ULPair(Stim,1)),1);
UCor(:,2) = nanmean(cellfun(@numel,ULPair(~Stim,1)),1);
UErr(:,1) = (nanmean(cellfun(@numel,ULPair(Stim,4)),1)) + (nanmean(cellfun(@numel,ULPair(Stim,7)),1));
UErr(:,2) = (nanmean(cellfun(@numel,ULPair(~Stim,4)),1)) + (nanmean(cellfun(@numel,ULPair(~Stim,7)),1));
sdUCor(:,1) = nanstd(cellfun(@numel,ULPair(Stim,1)),1);
sdUCor(:,2) = nanstd(cellfun(@numel,ULPair(~Stim,1)),1);
sdUErr(:,1) = (nanstd(cellfun(@numel,ULPair(Stim,4)),1)) + (nanstd(cellfun(@numel,ULPair(Stim,7)),1));
sdUErr(:,2) = (nanstd(cellfun(@numel,ULPair(~Stim,4)),1)) + (nanstd(cellfun(@numel,ULPair(~Stim,7)),1));
semUCor(:,1) = sdUCor(:,1)/sqrt(13);
semUCor(:,2) = sdUCor(:,2)/sqrt(13);
semUErr(:,1) = sdUErr(:,1)/sqrt(13);
semUErr(:,2) = sdUErr(:,2)/sqrt(13);

CorErr = [LCor LErr ; UCor UErr];
semCorErr = [semLCor semLErr ; semUCor semUErr];

LCor1(:,1) = cellfun(@numel,LPair(Stim,1));
LCor1(:,2) = cellfun(@numel,LPair(~Stim,1));
LErr2(:,1) = (cellfun(@numel,LPair(Stim,4))) + (cellfun(@numel,LPair(Stim,7)));
LErr2(:,2) = (cellfun(@numel,LPair(~Stim,4))) + (cellfun(@numel,LPair(~Stim,7)));

UCor1(:,1) = cellfun(@numel,ULPair(Stim,1));
UCor1(:,2) = cellfun(@numel,ULPair(~Stim,1));
UErr2(:,1) = (cellfun(@numel,ULPair(Stim,4))) + (cellfun(@numel,ULPair(Stim,7)));
UErr2(:,2) = (cellfun(@numel,ULPair(~Stim,4))) + (cellfun(@numel,ULPair(~Stim,7)));

% [h,p,ci,stats]=ttest2(LCor1(:,1), LCor1(:,2))
% [h,p,ci,stats]=ttest2(LErr2(:,1), LErr2(:,2))
% 
% [h,p,ci,stats]=ttest(LCor1(:,1), LErr2(:,1))
% [h,p,ci,stats]=ttest(LCor1(:,2), LErr2(:,2))
% 
% [h,p,ci,stats]=ttest2(UCor1(:,1), UCor1(:,2))
% [h,p,ci,stats]=ttest2(UErr2(:,1), UErr2(:,2))
% 
% [h,p,ci,stats]=ttest(UCor1(:,1), UErr2(:,1))
% [h,p,ci,stats]=ttest(UCor1(:,2), UErr2(:,2))

LNum(:,1) = nanmean(cellfun(@numel,LPair(Stim,:)),1);
LNum(:,2)= nanmean(cellfun(@numel,LPair(~Stim,:)),1);
sdLNum(:,1) = nanstd(cellfun(@numel,LPair(Stim,:)),1);
sdLNum(:,2)= nanstd(cellfun(@numel,LPair(~Stim,:)),1);
semLNum(:,1) = sdLNum(:,1)/sqrt(13);
semLNum(:,2) = sdLNum(:,2)/sqrt(13);
%Find mean, SD, SEM of UNum 
UNum(:,1) = nanmean(cellfun(@numel,ULPair(Stim,:)),1);
UNum(:,2) = nanmean(cellfun(@numel,ULPair(~Stim,:)),1);
sdUNum(:,1) = nanstd(cellfun(@numel,ULPair(Stim,:)),1);
sdUNum(:,2) = nanstd(cellfun(@numel,ULPair(~Stim,:)),1);
semUNum(:,1) = sdUNum(:,1)/sqrt(13);
semUNum(:,2) = sdUNum(:,2)/sqrt(13);

% STLPair  = nanmean(SubLPair(Stim,:),1);
% STULPair = nanmean(SubULPair(Stim,:),1);
% SHLPair  = nanmean(SubLPair(~Stim,:),1);
% SHULPair = nanmean(SubULPair(~Stim,:),1);

%Find size of LearnedPair and UnlearnedPair for S.E.M.
% sizeSTLB = sum(cellfun(@(c) size(c,1),LBins(Stim))); 
% sizeSHLB = sum(cellfun(@(c) size(c,1),LBins(~Stim)));
% sizeSTUB = sum(cellfun(@(c) size(c,1),UBins(Stim)));
% sizeSHUB = sum(cellfun(@(c) size(c,1),UBins(~Stim)));
%Find mean, SD, SEM of LBins and Ubins 
STLBins = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),LBins(Stim)','un',0)),1));
SHLBins = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),LBins(~Stim)','un',0)),1));
STUBins = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),UBins(Stim)','un',0)),1));
SHUBins = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),UBins(~Stim)','un',0)),1));

STLBins1 = cell2mat(cellfun(@(c) nanmean(c,1),LBins(Stim)','un',0));
SHLBins2 = cell2mat(cellfun(@(c) nanmean(c,1),LBins(~Stim)','un',0));

STLBins1 = cell2mat(cellfun(@(c) nanmean(c,1),LBins(Stim)','un',0));
SHLBins2 = cell2mat(cellfun(@(c) nanmean(c,1),LBins(~Stim)','un',0));
STUBins3 = cell2mat(cellfun(@(c) nanmean(c,1),UBins(Stim)','un',0));
SHUBins4 = cell2mat(cellfun(@(c) nanmean(c,1),UBins(~Stim)','un',0));


gLSTF = cell2mat(cellfun(@(c) nanmean(c,1),LBinsFHalf(Stim)','un',0));
gLSHF = cell2mat(cellfun(@(c) nanmean(c,1),LBinsFHalf(~Stim)','un',0));
gUSTF = cell2mat(cellfun(@(c) nanmean(c,1),UBinsFHalf(Stim)','un',0));
gUSHF = cell2mat(cellfun(@(c) nanmean(c,1),UBinsFHalf(~Stim)','un',0));

gLSTL = cell2mat(cellfun(@(c) nanmean(c,1),LBinsLHalf(Stim)','un',0));
gLSHL = cell2mat(cellfun(@(c) nanmean(c,1),LBinsLHalf(~Stim)','un',0));
gUSTL = cell2mat(cellfun(@(c) nanmean(c,1),UBinsLHalf(Stim)','un',0));
gUSHL = cell2mat(cellfun(@(c) nanmean(c,1),UBinsLHalf(~Stim)','un',0));

STLFBins = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),LBinsFHalf(Stim)','un',0)),1));
STUFBins = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),UBinsFHalf(Stim)','un',0)),1));

STLFBins1 = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),LBinsFHalf(Stim)','un',0)),1));
STUFBins2 = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),UBinsFHalf(Stim)','un',0)),1));

STLLBins = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),LBinsLHalf(Stim)','un',0)),1));
STULBins = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),UBinsLHalf(Stim)','un',0)),1));
SHLFBins = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),LBinsFHalf(~Stim)','un',0)),1));
SHUFBins = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),UBinsFHalf(~Stim)','un',0)),1));
SHLLBins = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),LBinsLHalf(~Stim)','un',0)),1));
SHULBins = squeeze(mean(cell2mat(cellfun(@(c) nanmean(c,1),UBinsLHalf(~Stim)','un',0)),1));
sdSTLFBins = nanstd(STLFBins);
sdSTUFBins = nanstd(STUFBins);
sdSTLLBins = nanstd(STLLBins);
sdSTULBins = nanstd(STULBins);
sdSHLFBins = nanstd(SHLFBins);
sdSHUFBins = nanstd(SHUFBins);
sdSHLLBins = nanstd(SHLLBins);
sdSHULBins = nanstd(SHULBins);
semSTLFBins = sdSTLFBins/sqrt(13);
semSTUFBins = sdSTUFBins/sqrt(13);
semSTLLBins = sdSTLLBins/sqrt(13);
semSTULBins = sdSTULBins/sqrt(13);
semSHLFBins = sdSHLFBins/sqrt(13);
semSHUFBins = sdSHUFBins/sqrt(13);
semSHLLBins = sdSHLLBins/sqrt(13);
semSHULBins = sdSHULBins/sqrt(13);

sdSTLBins = nanstd(STLBins);
sdSHLBins = nanstd(SHLBins);
sdSTUBins = nanstd(STUBins);
sdSHUBins = nanstd(SHUBins);
semSTLBins = sdSTLBins/sqrt(13);
semSHLBins = sdSHLBins/sqrt(13);
semSTUBins = sdSTUBins/sqrt(13);
semSHUBins = sdSHUBins/sqrt(13);

for j = 1:size(LBlockwDev,1)
    LFHalf(j) = nanmean(LBlockwDev(j,1:ceil(length(find(~isnan(LBlockwDev(j,:))))/2)));
    LLHalf(j) = nanmean(LBlockwDev(j,ceil(length(find(~isnan(LBlockwDev(j,:))))/2)+1:length(find(~isnan(LBlockwDev(j,:))))));
end

vStLFH = (LFHalf(Stim));
vShLFH = (LFHalf(~Stim));
vStLLH = (LLHalf(Stim));
vShLLH = (LLHalf(~Stim));

StLFH = nanmean(LFHalf(Stim));
ShLFH = nanmean(LFHalf(~Stim));
StLLH = nanmean(LLHalf(Stim));
ShLLH = nanmean(LLHalf(~Stim));

StLFHsd = nanstd(LFHalf(Stim));
ShLFHsd = nanstd(LFHalf(~Stim));
StLLHsd = nanstd(LLHalf(Stim));
ShLLHsd = nanstd(LLHalf(~Stim));


for t = 1:size(ULBlockwDev,1)
    ULFHalf(t) = nanmean(ULBlockwDev(t,1:ceil(length(find(~isnan(ULBlockwDev(t,:))))/2)));
    ULLHalf(t) = nanmean(ULBlockwDev(t,ceil(length(find(~isnan(ULBlockwDev(t,:))))/2)+1:length(find(~isnan(ULBlockwDev(t,:))))));
end

vStULFH = (ULFHalf(Stim));
vShULFH = (ULFHalf(~Stim));
vStULLH = (ULLHalf(Stim));
vShULLH = (ULLHalf(~Stim));

StULFH = nanmean(ULFHalf(Stim));
ShULFH = nanmean(ULFHalf(~Stim));
StULLH = nanmean(ULLHalf(Stim));
ShULLH = nanmean(ULLHalf(~Stim));

StULFHsd = nanstd(ULFHalf(Stim));
ShULFHsd = nanstd(ULFHalf(~Stim));
StULLHsd = nanstd(ULLHalf(Stim));
ShULLHsd = nanstd(ULLHalf(~Stim));

z = [StLFH StLLH; ShLFH ShLLH; StULFH StULLH; ShULFH ShULLH];
zstd = [StLFHsd StLLHsd; ShLFHsd ShLLHsd; StULFHsd StULLHsd; ShULFHsd ShULLHsd];


% %Two sample ttest 
% [h,p,ci,stats]=ttest2(LCor(:,1), LCor(:,2))
% [h,p,ci,stats]=ttest2(gLSTF(:,2,2), gLSHF(:,2,2))
% [h,p,ci,stats]=ttest2(gLSTF(:,3,2), gLSHF(:,3,2))
% [h,p,ci,stats]=ttest2(gLSTF(:,4,2), gLSHF(:,4,2))
% [h,p,ci,stats]=ttest2(gUSTF(:,1,2), gUSHF(:,1,2))
% [h,p,ci,stats]=ttest2(gUSTF(:,2,2), gUSHF(:,2,2))
% [h,p,ci,stats]=ttest2(gUSTF(:,3,2), gUSHF(:,3,2))
% [h,p,ci,stats]=ttest2(gUSTF(:,4,2), gUSHF(:,4,2))

% [h,p,ci,stats]=ttest2(STLBins1(:,1,2), SHLBins2(:,1,2))
% [h,p,ci,stats]=ttest2(STLBins1(:,2,2), SHLBins2(:,2,2))
% [h,p,ci,stats]=ttest2(STLBins1(:,3,2), SHLBins2(:,3,2))
% [h,p,ci,stats]=ttest2(STLBins1(:,4,2), SHLBins2(:,4,2))
% [h,p,ci,stats]=ttest2(STUBins3(:,1,2), SHUBins4(:,1,2))
% [h,p,ci,stats]=ttest2(STUBins3(:,2,2), SHUBins4(:,2,2))
% [h,p,ci,stats]=ttest2(STUBins3(:,3,2), SHUBins4(:,3,2))
% [h,p,ci,stats]=ttest2(STUBins3(:,4,2), SHUBins4(:,4,2))

% [h,p,ci,stats]=ttest(STLBins1(:,1,2), STLBins1(:,2,2))
% [h,p,ci,stats]=ttest(STLBins1(:,1,2), STLBins1(:,3,2))
% [h,p,ci,stats]=ttest(STLBins1(:,1,2), STLBins1(:,4,2))
% [h,p,ci,stats]=ttest(STLBins1(:,2,2), STLBins1(:,3,2))
% [h,p,ci,stats]=ttest(STLBins1(:,2,2), STLBins1(:,4,2))
% [h,p,ci,stats]=ttest(STLBins1(:,3,2), STLBins1(:,4,2))

% b = [repmat({'Stim'},[13,1]);repmat({'Sham'},[13,1])];
% t = table(b, [gLSTF(:,1,2); gLSHF(:,1,2)],[gLSTF(:,2,2) ; gLSHF(:,2,2)],[gLSTF(:,3,2);gLSHF(:,3,2)], [gLSTF(:,4,2);gLSHF(:,4,2)],...
%     'VariableNames', {'StimCondition', 'Bin1', 'Bin2', 'Bin3', 'Bin4'});
% Bins = table([1 2 3 4]', 'VariableNames', {'Measurements'});
% 
% rm = fitrm(t, 'Bin1-Bin4~StimCondition','WithinDesign',Bins)

% c = [repmat({'Stim'},[13,1]);repmat({'Sham'},[13,1])];
% s = table(c, [gUSTF(:,1,2); gUSHF(:,1,2)],[gUSTF(:,2,2) ; gUSHF(:,2,2)],[gUSTF(:,3,2);gUSHF(:,3,2)], [gUSTF(:,4,2);gUSHF(:,4,2)],...
%     'VariableNames', {'StimCondition', 'Bin1', 'Bin2', 'Bin3', 'Bin4'});
% Bins = table([1 2 3 4]', 'VariableNames', {'Measurements'});
% 
% rm = fitrm(s, 'Bin1-Bin4~StimCondition','WithinDesign',Bins)

% d = [repmat({'Stim'},[13,1]);repmat({'Sham'},[13,1])];
% f = table(d, [gLSTL(:,1,2); gLSHL(:,1,2)],[gLSTL(:,2,2) ; gLSHL(:,2,2)],[gLSTL(:,3,2);gLSHL(:,3,2)], [gLSTL(:,4,2);gLSHL(:,4,2)],...
%     'VariableNames', {'StimCondition', 'Bin1', 'Bin2', 'Bin3', 'Bin4'});
% Bins = table([1 2 3 4]', 'VariableNames', {'Measurements'});
% 
% rm = fitrm(f, 'Bin1-Bin4~StimCondition','WithinDesign',Bins)

% e = [repmat({'Stim'},[13,1]);repmat({'Sham'},[13,1])];
% g = table(e, [gUSTL(:,1,2); gUSHL(:,1,2)],[gUSTL(:,2,2) ; gUSHL(:,2,2)],[gUSTL(:,3,2);gUSHL(:,3,2)], [gUSTL(:,4,2);gUSHL(:,4,2)],...
%     'VariableNames', {'StimCondition', 'Bin1', 'Bin2', 'Bin3', 'Bin4'});
% Bins = table([1 2 3 4]', 'VariableNames', {'Measurements'});
% 
% rm = fitrm(g, 'Bin1-Bin4~StimCondition','WithinDesign',Bins)

d = [repmat({'Stim'},[13,1]);repmat({'Sham'},[13,1])];
f = table(d, [STLBins1(:,1,2); SHLBins2(:,1,2)],[STLBins1(:,2,2) ; SHLBins2(:,2,2)],[STLBins1(:,3,2);SHLBins2(:,3,2)], [STLBins1(:,4,2);SHLBins2(:,4,2)],...
    'VariableNames', {'StimCondition', 'Bin1', 'Bin2', 'Bin3', 'Bin4'});
Bins = table([1 2 3 4]', 'VariableNames', {'Measurements'});

rm = fitrm(f, 'Bin1-Bin4~StimCondition','WithinDesign',Bins)
% 
% e = [repmat({'Stim'},[13,1]);repmat({'Sham'},[13,1])];
% g = table(e, [STUBins3(:,1,2); SHUBins4(:,1,2)],[STUBins3(:,2,2) ; SHUBins4(:,2,2)],[STUBins3(:,3,2);SHUBins4(:,3,2)], [STUBins3(:,4,2);SHUBins4(:,4,2)],...
%     'VariableNames', {'StimCondition', 'Bin1', 'Bin2', 'Bin3', 'Bin4'});
% Bins = table([1 2 3 4]', 'VariableNames', {'Measurements'});
% 
% rm = fitrm(g, 'Bin1-Bin4~StimCondition','WithinDesign',Bins)
e = [repmat({'Stim'},[13,1]);repmat({'Sham'},[13,1])];
g = table(e, [STLBins1(:,1,2); SHLBins2(:,1,2)],[STLBins1(:,2,2) ; SHLBins2(:,2,2)],[STLBins1(:,3,2);SHLBins2(:,3,2)], [STLBins1(:,4,2);SHLBins2(:,4,2)],...
    'VariableNames', {'StimCondition', 'Bin1', 'Bin2', 'Bin3', 'Bin4'});
Bins = table([1 2 3 4]', 'VariableNames', {'Measurements'});

rm = fitrm(g, 'Bin1-Bin4~StimCondition','WithinDesign',Bins)

i =[STLBins1(:,1,2) SHLBins2(:,1,2)];
%%

figure(1); clf; hold on;
plot(squeeze(nanmean(nanmean(abs(UTrialDev(Stim,:,:)),1),2)));
plot(squeeze(nanmean(nanmean(abs(UTrialDev(~Stim,:,:)),1),2)));
plot(squeeze(nanmean(nanmean(abs(LTrialDev(Stim,:,:)),1),2)));
plot(squeeze(nanmean(nanmean(abs(LTrialDev(~Stim,:,:)),1),2)));
xlim([1 15]);
%title('Deviation from Target Window');
%legend({'UnLearned Stim' 'Unlearned Sham' 'Learned Stim' 'Learned Sham'});
%ylabel('ms');
set(get(gca,'ylabel'),'rotation',0);


figure(2); clf; hold on;
plot([nanmean(LBlockwDev(Stim,:),1) - nanmean(LBlockwDev(~Stim,:),1)]);
plot([nanmean(ULBlockwDev(Stim,:),1) - nanmean(ULBlockwDev(~Stim,:),1)]);
xlim([1 20])
%title('Learn/Unlearn(Stim Deviation - Sham Deviation from Target Window)');
%legend({'Learned(Stim-Sham)' 'Unlearned(Stim-Sham)'})
% xlabel('Block');
% ylabel('ms');
% set(get(gca,'ylabel'),'rotation',0);


figure(3); clf; hold on;
bar(z)
ngroups = size(z,1);
nbars = size(z,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for b = 1:nbars
    a = (1:ngroups) - groupwidth/2 + (2*b-1) * groupwidth/(2*nbars);
    errorbar(a, z(:,b), zstd(:,b), '.k');
end
xticks(1:4)
set (gca, 'XTickLabel', {'Stim Learned', 'Sham Learned', 'Stim UnLearned', 'Sham UnLearned'});
% title('Deviation from Target Window in Early vs. Late Learned and Unlearned Blocks');


figure(4); clf; hold on;
% plot(STLBins(:,2));
% plot(SHLBins(:,2));
% plot(STUBins(:,2));
% plot(SHUBins(:,2));
errorbar(STLBins(:,2),semSTLBins);
errorbar(SHLBins(:,2),semSHLBins);
%errorbar(STUBins(:,2),semSTUBins);
%errorbar(SHUBins(:,2),semSHUBins);
title('Absolute Deviations from Target Window')
legend({'Learned Stim' 'Learned Sham' 'UnLearned Stim' 'UnLearned Sham'});
xticks(1:4);
ylabel('ms')
xticklabels({'Bin 1' 'Bin 2' 'Bin 3' 'Bin 4'});


figure(5); clf; hold on;
bar(LearnedPair)
ngroups = size(LearnedPair,1);
nbars = size(LearnedPair,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for r = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*r-1) * groupwidth/(2*nbars);
    errorbar(x, LearnedPair(:,r), semLearnedPair(:,r), '.k');
end
xticks(1:9)
xticklabels({'C-C', 'C-E' 'C-L', 'E-C', 'E-E', 'E-L', 'L-C', 'L-E', 'L-L'});
% legend({'Stim' 'Sham'})
% title('Adjustment by Pairing in Learned Blocks')
% ylabel('ms')
%xlabel('Trial Pairings')
set(get(gca,'ylabel'),'rotation',0);
ylim([-1000 900])
% errorbar(LearnedPair, sdLearnedPair);
% Llabels = reshape(LNum',[],1);


figure(6); clf; hold on;
bar(UnlearnedPair)
ngroups = size(UnlearnedPair,1);
nbars = size(UnlearnedPair,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for r = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*r-1) * groupwidth/(2*nbars);
    errorbar(x, UnlearnedPair(:,r), semUnlearnedPair(:,r), '.k');
end
xticks(1:9)
xticklabels({'C-C', 'C-E' 'C-L', 'E-C', 'E-E', 'E-L', 'L-C', 'L-E', 'L-L'});
% legend({'Stim' 'Sham'})
% title('Adjustment by Pairing in Unlearned Blocks')
% ylabel('ms')
% xlabel('Trial Pairings')
set(get(gca,'ylabel'),'rotation',0);
ylim([-1100 900])
% Ulabels = reshape(UNum',[],1);


figure(7); clf; hold on;
bar(TotalPairs)
xticks(1:9)
xticklabels({'C-C', 'C-E' 'C-L', 'E-C', 'E-E', 'E-L', 'L-C', 'L-E', 'L-L'});
% title('Adjustment by Pairing in All Blocks')
% ylabel('ms')
% xlabel('Trial Pairings')
set(get(gca,'ylabel'),'rotation',0);
ylim([-1000 1000])


figure(8); clf; hold on;
bar(LNum)
ngroups = size(LNum,1);
nbars = size(LNum,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for r = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*r-1) * groupwidth/(2*nbars);
    errorbar(x, LNum(:,r), semLNum(:,r), '.k');
end
xticks(1:9)
xticklabels({'C-C', 'C-E' 'C-L', 'E-C', 'E-E', 'E-L', 'L-C', 'L-E', 'L-L'});
% legend({'Stim' 'Sham'})
% title('Number of Trials by Pairing in Learned Blocks')
% ylabel('No. of Trial Pairs')
% xlabel('Trial Pairings')


figure(9); clf; hold on;
bar(UNum)
ngroups = size(UNum,1);
nbars = size(UNum,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for r = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*r-1) * groupwidth/(2*nbars);
    errorbar(x, UNum(:,r), semUNum(:,r), '.k');
end
xticks(1:9)
xticklabels({'C-C', 'C-E' 'C-L', 'E-C', 'E-E', 'E-L', 'L-C', 'L-E', 'L-L'});
% legend({'Stim' 'Sham'})
% title('Number of Trials by Pairing in Unlearned Blocks')
% ylabel('No. of Trial Pairs')
% xlabel('Trial Pairings')
ylim([0 60])

% figure(10); clf; hold on;
% errorbar(STLFBins(:,2),semSTLFBins);
% errorbar(STUFBins(:,2),semSTUFBins);
% errorbar(STLLBins(:,2),semSTLLBins);
% errorbar(STULBins(:,2),semSTULBins);
% errorbar(SHLFBins(:,2),semSHLFBins);
% errorbar(SHUFBins(:,2),semSHUFBins);
% errorbar(SHLLBins(:,2),semSHLLBins);
% errorbar(SHULBins(:,2),semSHULBins);
% title('Absolute Deviations from Target Window')
% legend({'FH Learned Stim' 'FH Unlearned Stim' 'LH Learned Stim' 'LH Unlearned Stim' 'FH Learned Sham' 'FH Unlearned Sham' 'LH Learned Sham' 'LH Unlearned Sham'});
% xticks(1:4);
% ylabel('ms')
% xticklabels({'Bin 1' 'Bin 2' 'Bin 3' 'Bin 4'});

figure(11); clf; hold on;
errorbar(STLFBins(:,2),semSTLFBins);
errorbar(STUFBins(:,2),semSTUFBins);
errorbar(SHLFBins(:,2),semSHLFBins);
errorbar(SHUFBins(:,2),semSHUFBins);
title('Absolute Deviations from Target Window')
legend({'FH Learned Stim' 'FH Unlearned Stim' 'FH Learned Sham' 'FH Unlearned Sham'});
xticks(1:4);
ylim([-10 650])
ylabel('ms')
xticklabels({'Bin 1' 'Bin 2' 'Bin 3' 'Bin 4'});

figure(12); clf; hold on;
errorbar(STLLBins(:,2),semSTLLBins);
errorbar(STULBins(:,2),semSTULBins);
errorbar(SHLLBins(:,2),semSHLLBins);
errorbar(SHULBins(:,2),semSHULBins);
title('Absolute Deviations from Target Window')
legend({'LH Learned Stim' 'LH Unlearned Stim' 'LH Learned Sham' 'LH Unlearned Sham'});
xticks(1:4);
ylim([-10 650])
ylabel('ms')
xticklabels({'Bin 1' 'Bin 2' 'Bin 3' 'Bin 4'});

figure(13); clf; hold on;
bar(FHLearnedPair)
ngroups = size(FHLearnedPair,1);
nbars = size(FHLearnedPair,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for r = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*r-1) * groupwidth/(2*nbars);
    errorbar(x, FHLearnedPair(:,r), semFHLearnedPair(:,r), '.k');
end
xticks(1:9)
xticklabels({'C-C', 'C-E' 'C-L', 'E-C', 'E-E', 'E-L', 'L-C', 'L-E', 'L-L'});
legend({'Stim' 'Sham'})
title('FH Learned Pair')
ylabel('ms')
xlabel('Trial Pairings')
set(get(gca,'ylabel'),'rotation',0);
ylim([-1100 900])

figure(14); clf; hold on;
bar(LHLearnedPair)
ngroups = size(LHLearnedPair,1);
nbars = size(LHLearnedPair,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for r = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*r-1) * groupwidth/(2*nbars);
    errorbar(x, LHLearnedPair(:,r), semLHLearnedPair(:,r), '.k');
end
xticks(1:9)
xticklabels({'C-C', 'C-E' 'C-L', 'E-C', 'E-E', 'E-L', 'L-C', 'L-E', 'L-L'});
legend({'Stim' 'Sham'})
title('LH Learned Pair')
ylabel('ms')
xlabel('Trial Pairings')
set(get(gca,'ylabel'),'rotation',0);
ylim([-1100 900])

figure(15); clf; hold on;
bar(FHUnlearnedPair)
ngroups = size(FHUnlearnedPair,1);
nbars = size(FHUnlearnedPair,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for r = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*r-1) * groupwidth/(2*nbars);
    errorbar(x, FHUnlearnedPair(:,r), semFHUnlearnedPair(:,r), '.k');
end
xticks(1:9)
xticklabels({'C-C', 'C-E' 'C-L', 'E-C', 'E-E', 'E-L', 'L-C', 'L-E', 'L-L'});
legend({'Stim' 'Sham'})
title('FH Unlearned Pair')
ylabel('ms')
xlabel('Trial Pairings')
set(get(gca,'ylabel'),'rotation',0);
ylim([-1100 900])

figure(16); clf; hold on;
bar(LHUnlearnedPair)
ngroups = size(LHUnlearnedPair,1);
nbars = size(LHUnlearnedPair,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for r = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*r-1) * groupwidth/(2*nbars);
    errorbar(x, LHUnlearnedPair(:,r), semLHUnlearnedPair(:,r), '.k');
end
xticks(1:9)
xticklabels({'C-C', 'C-E' 'C-L', 'E-C', 'E-E', 'E-L', 'L-C', 'L-E', 'L-L'});
legend({'Stim' 'Sham'})
title('LH Unlearned Pair')
ylabel('ms')
xlabel('Trial Pairings')
set(get(gca,'ylabel'),'rotation',0);
ylim([-1100 900])

figure(17); clf; hold on;
bar(CorErr)
ngroups = size(CorErr,1);
nbars = size(CorErr,2);
groupwidth = min(0.8, nbars/(nbars+1.5));
for r = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*r-1) * groupwidth/(2*nbars);
    errorbar(x, CorErr(:,r), semCorErr(:,r), '.k');
end
xticks(1:2)
xticklabels({'Learned Blocks', 'Unlearned Blocks'});
legend({'Stim Corr - Corr' 'Sham Corr - Corr' 'Stim Err - Corr' 'Sham Err - Corr'})
% title('Number of Trials by Pairing in Learned Blocks')
ylabel('No. of Trial Pairs')
% xlabel('Trial Pairings')