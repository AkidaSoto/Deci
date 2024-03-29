function Deci_hctsa(Deci,info,data,params)


%% Seperate by Condition

for newloop=1:length(params.channel)
locks = data.locks;
events = data.events;
trlnum = data.trlnum;
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
    
    dataplaceholder = ft_selectdata(ccfg,data);
    
    for Lock = 1:length(params.Lock)
        display(' ')
        display(['---Starting Lock #' num2str(Lock) ': ' Deci.Analysis.LocksTitle{Lock} '---'])
        display(' ')
        info.Lock = Lock;
        
        cfg.offset = locks(ccfg.trials,Deci.Analysis.Locks(Lock));
        
        if all(isnan(cfg.offset)) || isempty(cfg.offset)
            continue
        end
        
        cfg.toilim = Deci.Analysis.Toilim;
        evalc('dat = ft_datashift2(cfg,dataplaceholder)');
        
        for lockstd = 1:size(data.trialinfo,2)
            lockers(lockstd)  =  mean(data.trialinfo(:,Lock) - data.trialinfo(:,lockstd));
        end
        info.lockers = lockers;
        
        
        %% Do Freq Analyses
        
        if ~strcmp(Deci.Analysis.Freq.method,'hilbert')
            fcfg = Deci.Analysis.Freq;
            fcfg.output='fourier';
            fcfg.pad = 'maxperlen';
            fcfg.scc = 0;
            fcfg.keeptapers = 'yes';
            fcfg.keeptrials = 'yes';
            fcfg.toi = Deci.Analysis.Toi(1):round(diff([data.time{1}(1) data.time{1}(2)]),5):Deci.Analysis.Toi(2);
            fcfg.gpu = Deci.GCom;
            fcfg.cpu = Deci.DCom;
            
            fcfg.channel = params.channel(newloop);
            
            Fourier = dc_freqanalysis(fcfg, dat);
            trllength = size(Fourier.fourierspctrm,1);
        else
            display('Applying Hilbert Transformation')
            fcfg = Deci.Analysis.Freq;
            fcfg.channel = params.channel(newloop);
            nyquist = data.fsample/2;
            
            freqs = params.foi;
            
            tempfreq = [];
            
            for foi = 1:length(freqs)
                
                hcfg = [];
                hcfg.bpfilter2 = 'yes';  %Modified implementation to work with MikexCohen's formula
                hcfg.bpfreq =[freqs(foi)-fcfg.width(foi) freqs(foi)+fcfg.width(foi)];
                hcfg.bpfiltord = round(fcfg.order*(data.fsample/hcfg.bpfreq(1)));
                hcfg.bpfilttype = 'firls';
                hcfg.transition_width = fcfg.transition_width;
                hcfg.hilbert = 'complex';
                
                evalc('hil = ft_preprocessing(hcfg,dat)');
                
                rcfg.latency = [params.Tfoi];
                Fo = ft_selectdata(rcfg,hil);
                
                tempfreq{foi}.fourierspctrm = permute(cell2mat(permute(Fo.trial,[3 1 2])),[3 1 4 2]);
                tempfreq{foi}.label = Fo.label;
                tempfreq{foi}.freq = freqs(foi);
                tempfreq{foi}.trialinfo = Fo.trialinfo;
                tempfreq{foi}.time = Fo.time{1}';
                tempfreq{foi}.dimord = 'rpt_chan_freq_time';
                
            end
            
            acfg.parameter = 'fourierspctrm';
            acfg.appenddim = 'freq';
            
            Fourier(Cond) = rmfield(ft_appendfreq(acfg,tempfreq{:}),'cfg');
            Fourier(Cond).dimord = 'rpt_chan_freq_time';
            trllength(Cond) = size(Fourier(Cond).fourierspctrm,1);
            
            
            
            
        end
        
        
        
        
    end
    
    
    %% Do hctsa UNPOOLED
    
    %Fourier(Cond) corresponds to the fourier of each condition
    
    %Deci.Analysis.CondTitle     = {'GG Correct'20 'GG Incorrect'10  'G0 Correct'10 'G0 Incorrect'0 ...
    %                               'N0 Correct'0 'N0 Incorrect'-10 'NN Correct'-10 'NN Incorrect'20};
    
    
    % Apply the function below to change the complex fourier into total power
    variable=Fourier;
    Fpowspctrm=abs(variable.fourierspctrm).^2;
    toi = variable.time >= round(params.bsltoi(1),4) & variable.time <= round(params.bsltoi(2),4);
    Bsl = nanmean(Fpowspctrm(:,:,:,toi),4);
    Bsl = repmat(Bsl,[1 1 1 size(Fpowspctrm,4)]);
    
    
    switch params.bsltype
        case 'none'
        case 'absolute'
            Fpowspctrm =  Fpowspctrm - Bsl;
        case 'relative'
            Fpowspctrm=  Fpowspctrm ./ Bsl;
        case 'relchange'
            Fpowspctrm = (Fpowspctrm - Bsl) ./ Bsl;
        case 'db'
            Fpowspctrm= 10*log10( Fpowspctrm ./ Bsl);
    end
    
    %Get Freqs and Time out of Fourier here using code below
    toi2 = round(variable.time,4) >= params.toi(1) &round(variable.time,4) <= params.toi(2);
    foi=round(variable.freq,4)>=params.foi(1) & round(variable.freq,4) <= params.foi(2);
    Fpowspctrm_int(:,:,:,:)=Fpowspctrm(:,:,foi,toi2);
    save newtestvar
    if  Cond==1
        Fpowspctrm_int1=Fpowspctrm_int;
    elseif Cond==2
        Fpowspctrm_int2=Fpowspctrm_int;
    end
    clear Fpowspctrm_int
end
%hctsa can only run after all conditions have been loaded
clear variable Fpowspctrm
fvar1_1=squeeze(mean(shiftdim(squeeze(Fpowspctrm_int1),1)));
fvar1_2=squeeze(mean(shiftdim(squeeze(Fpowspctrm_int2),1)));
corr=fvar1_1';
incorr=fvar1_2';

if params.pooled==2
    corr_ts=corr;
    incorr_ts=incorr;
    FCz_cor_trials=[corr_ts];
    FCz_inc_trials=[incorr_ts];
    allconditions=[FCz_cor_trials; FCz_inc_trials]; %add second condition once code confirmed
    timeSeriesData=allconditions;
    bothkeys=[repmat({'Correct'},1,size(FCz_cor_trials,1)) repmat({'Incorrect'},1,size(FCz_inc_trials,1))];
    uniquelabel=cell(1,size(timeSeriesData,1));
    keywords=bothkeys;
    
    for i=1:size(timeSeriesData,1)
        count=num2str(i);
        uniquelabel{1,i}=['Sample_' count];
    end
    
    labels=uniquelabel;
    heading0=['INP_' Deci.SubjectList{info.subject_list} 'CI.mat'];
    heading=['HCTSA_' Deci.SubjectList{info.subject_list} 'CI.mat'];
    save(heading0, 'timeSeriesData','labels','keywords')
    TS_Init(heading0,'INP_mops.txt','INP_ops.txt',[true,false,false],heading)
    TS_Compute(true,[],17:18,'error',heading);
    %heading_norm=TS_Normalize('mixedSigmoid',[0.4,1.0],heading);
    headingnew=['Subset_' heading];
    TS_Subset(heading,[],17:18,1,headingnew);
    heading_norm=headingnew;
    TS_LabelGroups(heading_norm)
    cfnParams=GiveMeDefaultClassificationParams(heading_norm);
    cfnParams.whatClassifier='svm_linear';
    all_tenfold_SVM_acc=TS_Classify(heading_norm);
    cfnParamsLDA=cfnParams;
    cfnParamsLDA.whatClassifier='linear';
    cfnParamsLDA.classifierText='linear classifier';
    all_tenfold_LDA_acc=TS_Classify(heading_norm,cfnParamsLDA);
    [a1,b1,~,~]=TS_TopFeatures(heading_norm,'classification',cfnParams);
    topfeat_SVM_mat=[a1 b1(a1)];
    [a2,b2,~,~]=TS_TopFeatures(heading_norm,'classification',cfnParamsLDA);
    topfeat_LDA_mat=[a2 b2(a2)];
    hctsa_norm_mat=load(heading_norm);
    
    
    mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'SWM_Del1Rsp2_UnpooledSubj_Gamma_Burstiness' ])
    save([Deci.Folder.Analysis filesep 'Extra' filesep 'SWM_Del1Rsp2_UnpooledSubj_Gamma_Burstiness' filesep Deci.SubjectList{info.subject_list}],'topfeat_SVM_mat','topfeat_LDA_mat','all_tenfold_SVM_acc','all_tenfold_LDA_acc','hctsa_norm_mat');
    
    
    
elseif params.pooled==1
    
    if info.subject_list==1
        corr_ts=corr;
        incorr_ts=incorr;
    elseif info.subject_list>1
        AllTr=load('TwoCondPooled_Tr');
        corr_ts=AllTr.corr_ts;
        incorr_ts=AllTr.incorr_ts;
        corr_ts(end+1:end+size(corr,1),:)=corr;
        incorr_ts(end+1:end+size(incorr,1),:)=incorr;
    end
    save('TwoCondPooled_Tr','corr_ts','incorr_ts')
    
    lastsub=length(Deci.SubjectList);
    if info.subject_list==lastsub
        FCz_cor_trials=[corr_ts];
        FCz_inc_trials=[incorr_ts];
        allconditions=[FCz_cor_trials; FCz_inc_trials];
        timeSeriesData=allconditions;
        bothkeys=[repmat({'Correct'},1,size(FCz_cor_trials,1)) repmat({'Incorrect'},1,size(FCz_inc_trials,1))];
        uniquelabel=cell(1,size(timeSeriesData,1));
        keywords=bothkeys;
        
        for i=1:size(timeSeriesData,1)
            count=num2str(i);
            uniquelabel{1,i}=['Sample_' count];
        end
        
        labels=uniquelabel;
        heading0='INP_Pooled_CI.mat';
        heading='HCTSA_Pooled_CI.mat';
        save(heading0, 'timeSeriesData','labels','keywords')
        TS_Init(heading0,'INP_mops.txt','INP_ops.txt',[true,false,false],heading)
        TS_Compute(true,[],[17 18 780],'missing',heading);
        heading_norm=TS_Normalize('mixedSigmoid',[0.4,1.0],heading);
         TS_Subset(heading_norm,[],[17 18],headingnew);
    heading_norm=headingnew;
        TS_LabelGroups(heading_norm)
        cfnParams=GiveMeDefaultClassificationParams(heading_norm);
        cfnParams.whatClassifier='svm_linear';
        all_tenfold_SVM_acc=TS_Classify(heading_norm);
        cfnParamsLDA=cfnParams;
        cfnParamsLDA.whatClassifier='linear';
        cfnParamsLDA.classifierText='linear classifier';
        all_tenfold_LDA_acc=TS_Classify(heading_norm,cfnParamsLDA);
        [a1,b1,~,~]=TS_TopFeatures(heading_norm,'classification',cfnParams);
        topfeat_SVM_mat=[a1 b1(a1)];
        [a2,b2,~,~]=TS_TopFeatures(heading_norm,'classification',cfnParamsLDA);
        topfeat_LDA_mat=[a2 b2(a2)];
        hctsa_norm_mat=load(heading_norm);
        
        
        mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'SWM_Del1Rsp2_PooledSubj_Gamma_Burstiness' ])
        save([Deci.Folder.Analysis filesep 'Extra' filesep 'SWM_Del1Rsp2_PooledSubj_Gamma_Burstiness' filesep 'Del3Imp2Del4_C6_CorrIncorr_PooledSubj'],'topfeat_SVM_mat','topfeat_LDA_mat','all_tenfold_SVM_acc','all_tenfold_LDA_acc','hctsa_norm_mat');
        
    end

    clear variable Fpowspctrm_int
end
savedchanfreq=[Deci.SubjectList{info.subject_list} 'Sing_chanfreq_burst.mat'];
if newloop~=1
saved_chanfreq_burst=load(savedchanfreq);
sing_chanfreq_burst=saved_chanfreq_burst.sing_chanfreq_burst;
end
sing_chanfreq_burst(:,newloop)=hctsa_norm_mat.TS_DataMat(:,1); %concatenate the channels
save(savedchanfreq)

clearvars -except Deci info data params
end

% saved_chanfreq_burst=load('sing_chanfreq_burst.mat');
% sing_chanfreq_burst=saved_chanfreq_burst.sing_chanfreq_burst;
% heading_norm=saved_chanfreq_burst.heading_norm;
% hctsa_norm_mat=load(heading_norm);
% sing_chanfreq_labels=hctsa_norm_mat.TimeSeries(:,3).Keywords;
% for ii1=1:length(sing_chanfreq_labels)
% if isequal(sing_chanfreq_labels{ii1},'Incorrect') && ~exist('theline')
%     theline=ii1;
% end
% end

% figure()
% %undersampling
% lowersize=size(sing_chanfreq_labels(1:theline-1),1);
% sing_chanfreq_burstnew(1:lowersize,:)=sing_chanfreq_burst(1:lowersize,:);
% sing_chanfreq_burstnew((lowersize+1):(lowersize*2+1),:)=sing_chanfreq_burst((theline:theline+lowersize),:);
% sing_chanfreq_labelsnew(1:lowersize,:)=sing_chanfreq_labels(1:lowersize,:);
% sing_chanfreq_labelsnew((lowersize+1):(lowersize*2+1),:)=sing_chanfreq_labels((theline:theline+lowersize),:);
% 
% %sing_chanfreq_burstnew_amp=sing_chanfreq_burstnew.^5; %amplifying
% for iii=1:length(saved_chanfreq_burst.params.channel)
% hold on
% set1=sing_chanfreq_burstnew((lowersize+1):end,iii);
% set2=sing_chanfreq_burstnew(1:lowersize,iii);
% outlier1=isoutlier(set1,'mean');
% outlier2=isoutlier(set2,'mean');
% set1(outlier1)=[];
% set2(outlier2)=[];
% errorbar(iii,mean(set1),(std(set1)/sqrt(length(set1))),'r*');
% errorbar(iii,mean(set2),(std(set2)/sqrt(length(set2))),'b*');
% 
% end
% mdl_linear=fitcdiscr(sing_chanfreq_burstnew,sing_chanfreq_labelsnew);
% %[mdl_auto,FitInfo,HyperparameterOptimizationResults]=fitcdiscr(sing_chanfreq_burst,sing_chanfreq_labels,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('Holdout',0.3,'AcquisitionFunctionName','expected-improvement-plus'));
% K_coeff=mdl_linear.Coeffs(1,2).Const;
% L_coeff=mdl_linear.Coeffs(1,2).Linear;
% % fcustomLDA= @(x1,x2) K_coeff + L_coeff(1)*x1 + L_coeff(2)*x2;
% %  boundline=fimplicit(fcustomLDA);
% %  boundline.Color='b';
% %  boundline.LineWidth=1.5;
% tester=predict(mdl_linear,sing_chanfreq_burstnew);
% simivec=zeros(length(tester),1);
% for ii2=1:length(sing_chanfreq_labelsnew)
%     if isequal(sing_chanfreq_labelsnew{ii2},tester{ii2})
%         simivec(ii2)=1;
%     elseif ~isequal(sing_chanfreq_labelsnew{ii2}, tester{ii2})
%         simivec(ii2)=0;
%     end
% end
% LDAaccSimple = 100*mean(simivec);
% LDAerror=loss(mdl_linear,sing_chanfreq_burstnew,sing_chanfreq_labelsnew);
% xlim([0 size(sing_chanfreq_burstnew,2)])
% hold off
end



