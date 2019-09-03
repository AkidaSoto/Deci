function out =  Deci_Corr(Deci,info,freq,params)

All_Trials = info.alltrials & info.allnonnans;

Magnitude = abs(freq.fourierspctrm.^2);
toi = round(freq.time,4) >= Deci.Analysis.Extra.Corr.Bsl(1) & round(freq.time,4) <= Deci.Analysis.Extra.Corr.Bsl(2);
Magnitude = 10*log10( Magnitude ./ mean(Magnitude(:,:,:,toi),4));

Magnitude = zscore(Magnitude);

Phase = angle(freq.fourierspctrm);


latency = reshape([freq.time [nan(-rem(length(freq.time),Deci.Analysis.Extra.Corr.timebin)+Deci.Analysis.Extra.Corr.timebin,1)]],[ceil([length(freq.time)]/Deci.Analysis.Extra.Corr.timebin) Deci.Analysis.Extra.Corr.timebin]);
latency = [min(latency,[],1);max(latency,[],1)];


for Var = 1:length(params.Variable)
    
    Variable = load([Deci.Folder.Analysis filesep 'Extra' filesep params.Variable{Var} filesep Deci.SubjectList{info.subject_list}  filesep filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}],params.Variable{Var});
    
    Type = fieldnames(Variable);
    
    if length(Variable.(Type{1})) ~= size(freq.trialinfo,1)
        
        corrs = 1:length(Variable.(Type{1}));
        
        if iscell(Variable.(Type{1}))
            
            if length(Variable.(Type{1}){1}) ~= size(freq.trialinfo,1)
                error('mismatch in trial count with Extra parameter')
            end
            
        end
        
        if isfield(Deci.Analysis.Extra.Corr,'Varnums')
            
            corrs = Deci.Analysis.Extra.Corr.Varnums;
        end
        
    else
        
        corrs = 1;
    end
    
    
    
    for num = corrs
        
        if iscell(Variable.(Type{1}))
            parameter = Variable.(Type{1}){num};
        else
            parameter = Variable.(Type{1});
        end
        
        parameter = zscore(parameter);
        
        for param = 1:length(params.Freq)
            
            R = [];
            P = [];
            Beta = [];
            
            for fois = 1:length(Deci.Analysis.Extra.Corr.freqbin)
                
                foi = freq.freq >= Deci.Analysis.Extra.Corr.freqbin{fois}(1) & freq.freq <= Deci.Analysis.Extra.Corr.freqbin{fois}(2);
                
                for tois = 1:Deci.Analysis.Extra.Corr.timebin
                    

                    toi = freq.time >= latency(1,tois) & freq.time <= latency(2,tois);
                    
                    switch params.Freq{param}
                        
                        case 'Magnitude'
                            
                            [r,p] = corrcoef(squeeze(nanmean(nanmean(Magnitude(:,:,foi,toi),3),4))',parameter);
                            
                            R(ismember(corrs,num),1,fois,tois) = r(1,2);
                            P(ismember(corrs,num),1,fois,tois) = p(1,2);
                            
                            
                            if Deci.Analysis.Extra.Corr.Regression
                                LM = fitlm(parameter,squeeze(nanmean(nanmean(Magnitude(:,:,foi,toi),3),4))');
                                Beta(ismember(corrs,num),1,fois,tois) = LM.Coefficients.Estimate(2);
                            end
                            
                        case 'Phase'
                            
                            [R(ismember(corrs,num),1,fois,tois),P(ismember(corrs,num),1,fois,tois)] =  circ_corrcl(squeeze(circ_mean(circ_mean(Phase(:,:,foi,toi),[],3),[],4)), parameter);
                            
                            if Deci.Analysis.Extra.Corr.Regression
                                LM = CircularRegression(parameter,squeeze(circ_mean(circ_mean(Phase(:,:,foi,toi),[],3),[],4)));
                                Beta(ismember(corrs,num),1,fois,tois) = LM(2);
                            end
                            
                            
                            
                    end
                    
                    if R(ismember(corrs,num),1,fois,tois) > 1
                        k = 0;
                    end
                    
                end
            end
            
            extracorr.label = freq.label;
            extracorr.freq = mean(cat(1,Deci.Analysis.Extra.Corr.freqbin{:}),2)';
            extracorr.time = mean(latency,1);
            extracorr.trialinfo = [corrs];
            extracorr.dimord =  freq.dimord;
            
            extracorr.powspctrm = R;
            R = extracorr;
            
            extracorr.powspctrm = P;
            P = extracorr;
            
            if Deci.Analysis.Extra.Corr.Regression
                extracorr.powspctrm = Beta;
                Beta = extracorr;
            end
            
            
            mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Corr' filesep params.Freq{param} '_' Type{1} filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}])
            save([Deci.Folder.Analysis filesep 'Extra' filesep 'Corr' filesep params.Freq{param} '_' Type{1} filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep info.Channels{info.ChanNum}],'R','P','Beta');
            
        end
    end
end
end