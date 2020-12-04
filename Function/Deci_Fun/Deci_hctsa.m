function Deci_hctsa(Deci,info,data,params)

%% Seperate by Condition

for Cond = 1:length(Deci.Analysis.Conditions)
    maxt = length(find(cellfun(@(c) any(ismember(Deci.Analysis.Conditions{Cond},c)), Deci.DT.Markers)));
    info.alltrials = find(sum(ismember(data.events,Deci.Analysis.Conditions{Cond}),2) == maxt);

        %% ignore all locks with missing nans
    if Deci.Analysis.IgnoreNanLocks
        minamountofnans = min(mean(isnan(data.locks(info.alltrials,:)),2));
        info.nanlocks = mean(isnan(data.locks(info.alltrials,:)),2) ~= minamountofnans;
        
        if any(info.nanlocks)
            display(['ignoring ' num2str(length(find(info.nanlocks))) ' trials with missing locks'])
        end
    else
        info.nanlocks = logical(size(info.alltrials));
    end
    
    %% Reject Arts
    ccfg.trials =  info.alltrials(~info.nanlocks);
    
    dataplaceholder{Cond} = ft_selectdata(ccfg,data);
end
%% Do hctsa

%Deci.Analysis.CondTitle     = {'GG Correct'20 'GG Incorrect'10  'G0 Correct'10 'G0 Incorrect'0 ...
%                               'N0 Correct'0 'N0 Incorrect'-10 'NN Correct'-10 'NN Incorrect'20};  


%ismember(dataplaceholder{2}.label,{'FCz'})
allchan_cond2=dataplaceholder{2}.trial;
[~,size_cond2]=size(allchan_cond2);
allchan_cond5=dataplaceholder{5}.trial;
[~,size_cond5]=size(allchan_cond5);
allchan_cond3=dataplaceholder{3}.trial;
[~,size_cond3]=size(allchan_cond3);
allchan_cond7=dataplaceholder{7}.trial;
[~,size_cond7]=size(allchan_cond7);

lo2=zeros(size_cond2,1);
for i=1:size_cond2
    [~,lo2(i)]=size(allchan_cond2{i});
end
lo5=zeros(size_cond5,1);
for i=1:size_cond5
    [~,lo5(i)]=size(allchan_cond5{i});
end
lo3=zeros(size_cond3,1);
for i=1:size_cond3
    [~,lo3(i)]=size(allchan_cond3{i});
end
lo7=zeros(size_cond7,1);
for i=1:size_cond7
    [~,lo7(i)]=size(allchan_cond7{i});
end

lo=min([min(lo2) min(lo5) min(lo3) min(lo7)]);


FCz_cond2_trials=zeros(size_cond2,min(lo));
for i=1:size_cond2
    FCz_cond2_trials(i,:)=allchan_cond2{i}(59,1:min(lo)); %change 59 to ismember notation later; also NOTE this has been normalized
end
FCz_cond5_trials=zeros(size_cond5,min(lo));
for i=1:size_cond5
    FCz_cond5_trials(i,:)=allchan_cond5{i}(59,1:min(lo)); %change 59 to ismember notation later; also NOTE this has been normalized
end
FCz_cond3_trials=zeros(size_cond3,min(lo));
for i=1:size_cond3
    FCz_cond3_trials(i,:)=allchan_cond3{i}(59,1:min(lo)); %change 59 to ismember notation later; also NOTE this has been normalized
end
FCz_cond7_trials=zeros(size_cond7,min(lo));
for i=1:size_cond7
    FCz_cond7_trials(i,:)=allchan_cond7{i}(59,1:min(lo)); %change 59 to ismember notation later; also NOTE this has been normalized
end

compare_case=2;
if compare_case==1
    FCz_cor_trials=[FCz_cond3_trials; FCz_cond7_trials];
    FCz_inc_trials=[FCz_cond2_trials; FCz_cond5_trials];
    
    [~,sizen2]=size(FCz_cor_trials);
    [~,sizen1]=size(FCz_inc_trials);
    allconditions=[FCz_cor_trials; FCz_inc_trials]; %add second condition once code confirmed
    timeSeriesData=allconditions;
    bothkeys=[repmat({'Correct'},1,size(FCz_cor_trials,1)) repmat({'Incorrect'},1,size(FCz_inc_trials,1))];
    %save bothkeys
    uniquelabel=cell(1,size(timeSeriesData,1));
    keywords=bothkeys;
elseif compare_case==2
    [~,sizen2]=size(FCz_cond2_trials);
    [~,sizen3]=size(FCz_cond3_trials);
    [~,sizen5]=size(FCz_cond5_trials);
    [~,sizen7]=size(FCz_cond7_trials);
    allconditions=[FCz_cond2_trials; FCz_cond3_trials; FCz_cond5_trials; FCz_cond7_trials];
    timeSeriesData=allconditions;
    allkeys=[repmat({'Condition 2'},1,size(FCz_cond2_trials,1)) repmat({'Condition 3'},1,size(FCz_cond3_trials,1)) repmat({'Condition 5'},1,size(FCz_cond5_trials,1)) repmat({'Condition 7'},1,size(FCz_cond7_trials,1))];
    %save allkeys
    keywords=allkeys;
end
    %save timeSeriesData
    for i=1:size(timeSeriesData,1)
        count=num2str(i);
        uniquelabel{1,i}=['Sample_' count];
    end
    labels=uniquelabel;
    save('INP_test.mat', 'timeSeriesData','labels','keywords')
    TS_Init('INP_test.mat')

% TS_Compute();
% TS_Normalize('mixedSigmoid',[0.4,1.0]);
% TS_LabelGroups('norm')
% TS_PlotTimeSeries('norm')
% TS_PlotDataMatrix('norm')
% TS_Cluster()
% TS_PlotDataMatrix('norm')
% TS_PlotLowDim('norm','pca')
% TS_PlotLowDim('norm','tsne')
% TS_Classify('HCTSA_N.mat')
% cfnParams=GiveMeDefaultClassificationParams('HCTSA_N.mat');
% numNulls=100;
% %TS_Classify('HCTSA_N.mat',cfnParams,num Nulls,'doParallel',true)
% featuresets= {'notLocationDependent','locationDependent','notLengthDependent','lengthDependent','notSpreadDependent','spreadDependent'};
% TS_CompareFeatureSets('norm',cfnParams,featuresets)
% TS_ClassifyLowDim('norm')
% TS_TopFeatures('norm', 'classification')

end