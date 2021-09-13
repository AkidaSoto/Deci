function Deci_hctsa(Deci,info,data,params)


%% Seperate by Condition
if params.channel{1}=='all'
    params.channel=data.label; %if 'all' is specified, then all channels that the data has are used
end
for freqbandloop=1:size(params.foi,1) %frequency band loop
    for chanloop=1:length(params.channel) %channel loop nested within frequency band loop
        locks = data.locks;
        events = data.events;
        trlnum = data.trlnum;
        for Cond = 1:length(Deci.Analysis.Conditions) %condition loop nested within channel loop
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
                    
                    fcfg.channel = params.channel(chanloop);
                    
                    Fourier = dc_freqanalysis(fcfg, dat);
                    trllength = size(Fourier.fourierspctrm,1);
                else
                    display('Applying Hilbert Transformation')
                    fcfg = Deci.Analysis.Freq;
                    fcfg.channel = params.channel(chanloop);
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
            
            
            %% Do hctsa Unpooled (subject by subject)
            
            %Fourier(Cond) corresponds to the fourier of each condition
            % Applying the function below to change the complex fourier into total power
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
            foi=round(variable.freq,4)>=params.foi(freqbandloop,1) & round(variable.freq,4) <= params.foi(freqbandloop,2);
            Fpowspctrm_int(:,:,:,:)=Fpowspctrm(:,:,foi,toi2);
            
            if  Cond==1    %using if statement to generate 2 different variables for 2 conditions of SWM working memory task
                Fpowspctrm_int1=Fpowspctrm_int;
            elseif Cond==2
                Fpowspctrm_int2=Fpowspctrm_int;
            end
            
            clear Fpowspctrm_int
        end %condition loop ends here
        %hctsa can only run after all conditions have been loaded
        clear variable Fpowspctrm
        fvar1_1=squeeze(mean(shiftdim(squeeze(Fpowspctrm_int1),1)));
        fvar1_2=squeeze(mean(shiftdim(squeeze(Fpowspctrm_int2),1)));
        corr=fvar1_1';
        incorr=fvar1_2';
        
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
        %heading_norm=TS_Normalize('mixedSigmoid',[0.4,1.0],heading); %we
        %stopped normalizing in order to use the raw burstiness values
        headingnew=['Subset_' heading];
        TS_Subset(heading,[],17:18,1,headingnew)
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
        
        burstiness_per_trial=hctsa_norm_mat.TS_DataMat(:,1);
        sing_chanfreq_labels=hctsa_norm_mat.TimeSeries(:,3).Keywords;
        counter1=0;
        counter2=0;
        for i2=1:length(sing_chanfreq_labels)
            if isequal(sing_chanfreq_labels{i2},'Correct')
                counter1=counter1+1;
                newmp_correct(counter1)=burstiness_per_trial(i2);
            elseif isequal(sing_chanfreq_labels{i2},'Incorrect')
                counter2=counter2+1;
                newmp_incorrect(counter2)=burstiness_per_trial(i2);
            end
        end
        newmp_mcorr(info.subject_list,chanloop)=mean(newmp_correct);
        newmp_minco(info.subject_list,chanloop)=mean(newmp_incorrect);
        
        
        mkdir([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{info.subject_list} filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))]  filesep 'Burstiness per trial'])
        
        save([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{info.subject_list} filesep 'SWM_Del1Rsp2_Gamma_Burstiness'  filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'Burstiness per trial' filesep params.channel{chanloop} ],'burstiness_per_trial','sing_chanfreq_labels');
        vali.powspctrm=newmp_mcorr(info.subject_list,chanloop);
        vali.label=params.channel{chanloop};
        vali.dimord='chan';
        mkdir([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{info.subject_list} filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))]  filesep 'Burstiness per channel' filesep 'correct'])
        save([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{info.subject_list} filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))]  filesep 'Burstiness per channel' filesep 'correct' filesep params.channel{chanloop}], 'vali')
        clear vali
        vali.powspctrm=newmp_minco(info.subject_list,chanloop);
        vali.label=params.channel{chanloop};
        vali.dimord='chan';
        mkdir([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{info.subject_list} filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))]  filesep 'Burstiness per channel' filesep 'incorrect'])
        save([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{info.subject_list} filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))]  filesep 'Burstiness per channel' filesep 'incorrect' filesep params.channel{chanloop}], 'vali')
        clear vali
        vali.powspctrm=newmp_minco(info.subject_list,chanloop)-newmp_mcorr(info.subject_list,chanloop);
        vali.label=params.channel{chanloop};
        vali.dimord='chan';
        mkdir([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{info.subject_list} filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))]  filesep 'Burstiness per channel' filesep 'incorrect_minus_correct'])
        save([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{info.subject_list} filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))]  filesep 'Burstiness per channel' filesep 'incorrect_minus_correct' filesep params.channel{chanloop}], 'vali')
        
        
        
        savedchanfreq=[Deci.SubjectList{info.subject_list} 'Sing_chanfreq_burst.mat'];
        if chanloop~=1
            saved_chanfreq_burst=load(savedchanfreq);
            sing_chanfreq_burst=saved_chanfreq_burst.sing_chanfreq_burst;
        end
        sing_chanfreq_burst(:,chanloop)=hctsa_norm_mat.TS_DataMat(:,1); %concatenate the channels
        save(savedchanfreq)
        
        clearvars -except Deci info data params freqbandloop
    end
end
%% LDA Classifier
if Deci.SubjectList{info.subject_list}==Deci.SubjectList{end}
    for si=1:length(Deci.SubjectList)
        heading_norm=['Subset_HCTSA_' Deci.SubjectList{info.subject_list} 'CI.mat'];
        hctsa_norm_mat=load(heading_norm);
        sing_chanfreq_labels=hctsa_norm_mat.TimeSeries(:,3).Keywords;
        savedchanfreq=[Deci.SubjectList{info.subject_list} 'Sing_chanfreq_burst.mat'];
        saved_chanfreq_burst=load(savedchanfreq);
        sing_chanfreq_burst=saved_chanfreq_burst.sing_chanfreq_burst;
        for ii1=1:length(sing_chanfreq_labels)
            if isequal(sing_chanfreq_labels{ii1},'Incorrect') && ~exist('theline')
                theline=ii1;
            end
        end
        figure()
        %undersampling as there are more incorrects than correct
        if theline<(0.5*size(sing_chanfreq_labels,1))
            lowersize=size(sing_chanfreq_labels(1:theline-1),1);
            sing_chanfreq_burstnew(1:lowersize,:)=sing_chanfreq_burst(1:lowersize,:);
            sing_chanfreq_burstnew((lowersize+1):(lowersize*2+1),:)=sing_chanfreq_burst((theline:theline+lowersize),:);
            sing_chanfreq_labelsnew(1:lowersize,:)=sing_chanfreq_labels(1:lowersize,:);
            sing_chanfreq_labelsnew((lowersize+1):(lowersize*2+1),:)=sing_chanfreq_labels((theline:theline+lowersize),:);
        elseif theline>=(0.5*size(sing_chanfreq_labels,1))
            lowersize=size(sing_chanfreq_labels(theline:end),1);
            sing_chanfreq_burstnew(1:lowersize,:)=sing_chanfreq_burst(1:lowersize,:);
            sing_chanfreq_burstnew((lowersize+1):(lowersize*2+1),:)=sing_chanfreq_burst((theline:theline+lowersize),:);
            sing_chanfreq_labelsnew(1:lowersize,:)=sing_chanfreq_labels(1:lowersize,:);
            sing_chanfreq_labelsnew((lowersize+1):(lowersize*2+1),:)=sing_chanfreq_labels((theline:theline+lowersize),:);
            
        end
        
        
        
        for iii=1:size(sing_chanfreq_burstnew,2)
            set1=sing_chanfreq_burstnew((lowersize+1):end,iii); %INCORRECT
            set2=sing_chanfreq_burstnew(1:lowersize,iii); %CORRECT
            outlier1=isoutlier(set1,'mean');
            outlier2=isoutlier(set2,'mean');
            set1(outlier1)=[];
            set2(outlier2)=[];
        end
        
        sing_chanfreq_burstnewtrain=sing_chanfreq_burstnew;
        sing_chanfreq_labelsnewtrain=sing_chanfreq_labelsnew;
        
        mdl_linear=fitcdiscr(sing_chanfreq_burstnewtrain,sing_chanfreq_labelsnewtrain,'DiscrimType','linear');
        K_coeff(:,si)=mdl_linear.Coeffs(1,2).Const;
        L_coeff(:,si)=mdl_linear.Coeffs(1,2).Linear;
        mdl_cross=crossval(mdl_linear);
        tester=predict(mdl_linear,sing_chanfreq_burstnew);
        simivec=zeros(length(tester),1);
        for ii2=1:length(sing_chanfreq_labelsnew)
            if isequal(sing_chanfreq_labelsnew{ii2},tester{ii2})
                simivec(ii2)=1;
            elseif ~isequal(sing_chanfreq_labelsnew{ii2}, tester{ii2})
                simivec(ii2)=0;
            end
        end
        LDAaccSimple(:) = 100*mean(simivec);
        LDAerror(:)=loss(mdl_linear,sing_chanfreq_burstnew,sing_chanfreq_labelsnew);
        LDAerror_cross(:,si)=kfoldLoss(mdl_cross);
        
        %for si=1:length(Deci.SubjectList)
        for freqbandloop=1:size(params.foi,1)
            for chanloop= 1:length(params.channel)
                vali.powspctrm=L_coeff(chanloop,si);
                vali.label=params.channel{chanloop};
                vali.dimord='chan';
                mkdir([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{info.subject_list} filesep 'SWM_Del1Rsp2_Gamma_Burstiness'  filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'LDA_coefficient_per_channel']);
                save([Deci.Folder.Analysis filesep 'Extra' filesep Deci.SubjectList{info.subject_list} filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'LDA_coefficient_per_channel' filesep params.channel{chanloop} ], 'vali');
            end
        end
        
        clear sing_chanfreq_burstnew sing_chanfreq_labelsnew sing_chanfreqlabels sing_chanfreq_burst theline simivec mdl_linear tester set1 set2 outlier1 outlier2
    end
    %% Cluster Analysis
    LDA_cross_acc=(1-LDAerror_cross);
    select_cross=L_coeff(:,LDA_cross_acc>0.5);
    select_cross_lab=Deci.SubjectList(LDA_cross_acc>0.5);%finding the appropriate subject names
    varplacehold=LDA_cross_acc(LDA_cross_acc>0.5);
    for i=1:length(varplacehold)
        select_cross_with_labels{i,1}=varplacehold(i);
        select_cross_with_labels{i,2}=select_cross_lab{i};
        newmpacc_mcorr=newmp_mcorr(i,:);
        newmpacc_minco=newmp_minco(i,:);
    end
    
    select_coeff=select_cross';
    [ids2,c2]=kmeans(select_coeff,2); %kmeans clustering specified for 2 clusters
    for subclusters=1:max(ids2) %clusters loop to produce and save plottable data for each cluster
        subsel=L_coeff(:,ids2==subclusters); %choosing subjects within this cluster
        newmpacc_mcorr=newmp_mcorr(ids2==subclusters,:);
        newmpacc_minco=newmp_minco(ids2==subclusters,:);
        newmpcorr_totmean=mean(newmpacc_mcorr,1);
        newmpincorr_totmean=mean(newmpacc_minco,1);
        
        LDA_coeff_significant_ini=mean(subsel,2);
        for freqbandloop=1:size(params.foi,1) %frequency loop nested within clusters loop
            for chanloop= 1:length(params.channel) %channels loop nested within frequency loop
                val.powspctrm=LDA_coeff_significant_ini(chanloop);
                val.label=params.channel{chanloop};
                val.dimord='chan';
                mkdir([Deci.Folder.Analysis filesep 'Extra' filesep ['Cluster' num2str(subclusters)] filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'Stats'])
                save([Deci.Folder.Analysis filesep 'Extra' filesep ['Cluster' num2str(subclusters)] filesep 'SWM_Del1Rsp2_Gamma_Burstiness'  filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'Stats' filesep 'Stat_Accuracy' ], 'select_cross_with_labels');
                mkdir([Deci.Folder.Analysis filesep 'Extra' filesep ['Cluster' num2str(subclusters)] filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'LDA_Coefficient'])
                save([Deci.Folder.Analysis filesep 'Extra' filesep ['Cluster' num2str(subclusters)] filesep 'SWM_Del1Rsp2_Gamma_Burstiness'  filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'LDA_Coefficient' filesep params.channel{chanloop} ], 'val');
                
                
                vali.powspctrm=newmpcorr_totmean(chanloop);
                vali.label=params.channel{chanloop};
                vali.dimord='chan';
                mkdir([Deci.Folder.Analysis filesep 'Extra' filesep ['Cluster' num2str(subclusters)] filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'burstiness_per_channel' filesep 'correct'])
                save([Deci.Folder.Analysis filesep 'Extra' filesep ['Cluster' num2str(subclusters)] filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'burstiness_per_channel' filesep 'correct' filesep params.channel{chanloop}], 'vali')
                clear vali
                vali.powspctrm=newmpincorr_totmean(chanloop);
                vali.label=params.channel{chanloop};
                vali.dimord='chan';
                mkdir([Deci.Folder.Analysis filesep 'Extra' filesep ['Cluster' num2str(subclusters)] filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'burstiness_per_channel' filesep 'incorrect'])
                save([Deci.Folder.Analysis filesep 'Extra' filesep ['Cluster' num2str(subclusters)] filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'burstiness_per_channel' filesep 'incorrect' filesep params.channel{chanloop}], 'vali')
                clear vali
                vali.powspctrm=newmpincorr_totmean(chanloop)-newmpcorr_totmean(chanloop);
                vali.label=params.channel{chanloop};
                vali.dimord='chan';
                mkdir([Deci.Folder.Analysis filesep 'Extra' filesep ['Cluster' num2str(subclusters)] filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'burstiness_per_channel' filesep 'incorrect_minus_correct'])
                save([Deci.Folder.Analysis filesep 'Extra' filesep ['Cluster' num2str(subclusters)] filesep 'SWM_Del1Rsp2_Gamma_Burstiness' filesep ['Freq- ' num2str(params.foi(freqbandloop,:))] filesep 'burstiness_per_channel' filesep 'incorrect_minus_correct' filesep params.channel{chanloop}], 'vali')
                
            end
        end
    end
end
end



