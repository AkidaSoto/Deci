function Deci_hctsa(Deci,params)

%% Deci Checks
Dims = {'Topo' 'Square' 'Wire' 'Bar' 'MTopo'};
[Deci, info] = dc_plotcheck(Deci,Dims);
info.isfreq = true;
info.isconn = false;

%% Load
disp('----------------------');
display(' ')

display(['Plotting ' Deci.Plot.Freq.Type]);

info.extension = ['Freq_' Deci.Plot.Freq.Type];
info.parameter = 'powspctrm';
info.variable = 'freq';

[Subjects,info] =  dc_plotload(Deci,info);

%% Baseline Correction
display(' ');
display(['Using Lock: ' Deci.Plot.Lock]); 
display(['Using Ref: ' Deci.Plot.BslRef ' at times ' strrep(regexprep(num2str(Deci.Plot.Bsl),' +',' '),' ','-')]);

Subjects = dc_plotbsl(Deci,Subjects,info);

%% Math
if ~isempty(Deci.Plot.Math)
     [Deci, Subjects] = dc_math(Deci,Subjects,info);
end

%% Data Management
if size(Subjects,1) == 1
    Deci.Plot.GrandAverage = false;
end

if Deci.Plot.GrandAverage
    if any(~isnan(info.trllen))
        info.trlstd = nanstd(info.trllen,[],1);
        info.trllen = nanmean(info.trllen,1);
    end
    
    %if any(~isnan(info.lockers))
        info.lockersstd = nanstd(info.lockers,[],1);
        info.lockers = nanmean(info.lockers,1);
    %end
end

% Mean across freqs

cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';

%% DataSegmentation

for conds = 1:size(Subjects,2)
    
    if Deci.Plot.GrandAverage
        facfg.parameter =  info.parameter;
        facfg.type = 'mean';
        facfg.keepindividual = 'yes';
        
        if info.isfreq
            evalc('AvgData{1,conds} = ft_freqgrandaverage(facfg,Subjects{:,conds});');
        else
            evalc('AvgData{1,conds} = ft_timelockgrandaverage(facfg,Subjects{:,conds});');
        end
        AvgData{1,conds} = rmfield(AvgData{1,conds},'cfg');
        
        if Deci.Plot.GroupLevel
            AvgData2{2,conds} = AvgData{1,conds};
            AvgData2{1,conds} = AvgData{1,conds};
            AvgData2{2,conds}.powspctrm = AvgData{1,conds}.powspctrm(Deci.Plot.Groups{2},:,:,:);
            AvgData2{1,conds}.powspctrm = AvgData{1,conds}.powspctrm(Deci.Plot.Groups{1},:,:,:);
        end
        
    else
        Deci.Plot.Stat.do = false;
        AvgData(:,conds) = Subjects(:,conds);
    end
    
    for subj = 1:size(AvgData(:,conds),1)
        tcfg = [];
        tcfg.nanmean = Deci.Plot.nanmean;
        
        tcfg.latency = Deci.Plot.Wire.Toi;
        tcfg.frequency = Deci.Plot.Wire.Foi;
        tcfg.channel = Deci.Plot.Wire.Channel;
        
        Segdata{subj,conds} = ft_selectdata(tcfg,AvgData{subj,conds});
        
        %powspctrm [chan_freq_time] 1_1_1001
    end
end

save Segdata
%% Do hctsa

%Deci.Analysis.CondTitle     = {'GG Correct'20 'GG Incorrect'10  'G0 Correct'10 'G0 Incorrect'0 ...
%                               'N0 Correct'0 'N0 Incorrect'-10 'NN Correct'-10 'NN Incorrect'20};  



% compare_case=2;
% if compare_case==1
%     FCz_cor_trials=[FCz_cond3_trials; FCz_cond7_trials];
%     FCz_inc_trials=[FCz_cond2_trials; FCz_cond5_trials];
%     
%     [~,sizen2]=size(FCz_cor_trials);
%     [~,sizen1]=size(FCz_inc_trials);
%     allconditions=[FCz_cor_trials; FCz_inc_trials]; %add second condition once code confirmed
%     timeSeriesData=allconditions;
%     bothkeys=[repmat({'Correct'},1,size(FCz_cor_trials,1)) repmat({'Incorrect'},1,size(FCz_inc_trials,1))];
%     %save bothkeys
%     uniquelabel=cell(1,size(timeSeriesData,1));
%     keywords=bothkeys;
% elseif compare_case==2
%     [~,sizen2]=size(FCz_cond2_trials);
%     [~,sizen3]=size(FCz_cond3_trials);
%     [~,sizen5]=size(FCz_cond5_trials);
%     [~,sizen7]=size(FCz_cond7_trials);
%     allconditions=[FCz_cond2_trials; FCz_cond3_trials; FCz_cond5_trials; FCz_cond7_trials];
%     timeSeriesData=allconditions;
%     allkeys=[repmat({'Condition 2'},1,size(FCz_cond2_trials,1)) repmat({'Condition 3'},1,size(FCz_cond3_trials,1)) repmat({'Condition 5'},1,size(FCz_cond5_trials,1)) repmat({'Condition 7'},1,size(FCz_cond7_trials,1))];
%     %save allkeys
%     keywords=allkeys;
% end
%     %save timeSeriesData
%     for i=1:size(timeSeriesData,1)
%         count=num2str(i);
%         uniquelabel{1,i}=['Sample_' count];
%     end
%     labels=uniquelabel;
%     save('INP_test.mat', 'timeSeriesData','labels','keywords')
%     TS_Init('INP_test.mat')

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