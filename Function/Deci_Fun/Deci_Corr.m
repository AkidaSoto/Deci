function out =  Deci_Corr(Deci,info,freq,params)

All_Trials = info.alltrials & info.allnonnans;

[toi] = downsample(freq.time,round(1/diff([freq.time(1),freq.time(2)]))/Deci.Analysis.Extra.Corr.Downsample);
toi = ismember(freq.time,toi);
freq.fourierspctrm = freq.fourierspctrm(:,:,:,toi);
freq.time = freq.time(toi);

Magnitude = abs(freq.fourierspctrm.^2);
toi = round(freq.time,4) >= Deci.Analysis.Extra.Corr.Bsl(1) & round(freq.time,4) <= Deci.Analysis.Extra.Corr.Bsl(2);
Magnitude = 10*log10( Magnitude ./ mean(Magnitude(:,:,:,toi),4));

Magnitude = zscore(Magnitude);

Phase = circ_rad2ang(angle(freq.fourierspctrm));


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
            
            for foi = 1:length(freq.freq)
                
                for toi = 1:length(freq.time)
                    
                    switch params.Freq{param}
                        
                        case 'Magnitude'
                            [r,p] = corrcoef(squeeze(Magnitude(:,:,foi,toi)),parameter);
                            
                            R(num,1,foi,toi) = r(1,2);
                            P(num,1,foi,toi) = p(1,2);
                            
                            
                            if Deci.Analysis.Extra.Corr.Regression
                                LM = fitlm(parameter,squeeze(Magnitude(:,:,foi,toi)));
                                Beta(num,1,foi,toi) = LM.Coefficients.Estimate(2);
                            end
                            
                        case 'Phase'
                            
                            [R(num,1,foi,toi),P(num,1,foi,toi)] =  circ_corrcl(squeeze(Phase(:,:,foi,toi)), parameter);
                            
                            if Deci.Analysis.Extra.Corr.Regression
                                LM = CircularRegression(parameter,angle(freq.fourierspctrm(:,:,foi,toi)));
                                Beta(num,1,foi,toi) = LM(2);
                            end
                    end
                    
                    
                end
            end
            
            extracorr.label = freq.label;
            extracorr.freq = freq.freq;
            extracorr.time = freq.time;
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