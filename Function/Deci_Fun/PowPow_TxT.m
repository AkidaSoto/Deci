function PowPow_TxT(Deci,info,Fourier,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by CTGill
% last edited 1/8/20
%
% PowPowTxT.m computes the power time series for two electrodes, and then
% computes the correlation coefficient between time varying power of the
% two electrodes in sliding time segments within a trial, and then for each
% trial.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ismember(Deci.Analysis.CondTitle(info.Cond),Deci.Analysis.Connectivity.cond)
    
    if Deci.Analysis.Laplace
        mkdir([Deci.Folder.Analysis filesep 'PowPow_TxT'  filesep 'Laplacian' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep [Deci.Analysis.Connectivity.chan1{1} '_' Deci.Analysis.Connectivity.chan2{1}]]);
    else
        mkdir([Deci.Folder.Analysis filesep 'PowPow_TxT' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep [Deci.Analysis.Connectivity.chan1{1} '_' Deci.Analysis.Connectivity.chan2{1}]]);
    end
    
    times2save = -.5:.02:1.5;
    time_window_options = linspace(1.5,3.5,length(Deci.Analysis.Freq.foi));
    
    %time in indicies
    times2save_idx = dsearchn(Fourier.time', times2save');
    
    
    chan_idx = zeros(1,2);
    chan_idx(1) = find(strcmpi(Deci.Analysis.Connectivity.chan1,Fourier.label));
    chan_idx(2) = find(strcmpi(Deci.Analysis.Connectivity.chan2,Fourier.label));
    
    time_window = zeros(length(Deci.Analysis.Connectivity.freq),1);
    
    for fh = 1:length(Deci.Analysis.Connectivity.freq)
        switch Deci.Analysis.Connectivity.freq{fh}
            case 'theta'
                HF{fh} = [4 8];
                time_window(fh) = time_window_options(dsearchn(Fourier.freq',HF{fh}(1)));
            case 'beta'
                HF{fh} = [12.5 30];
                time_window(fh) = time_window_options(dsearchn(Fourier.freq',HF{fh}(1)));
            case 'alpha'
                HF{fh} =[8 12.5];
                time_window(fh) = time_window_options(dsearchn(Fourier.freq',HF{fh}(1)));
            case 'lowgamma'
                HF{fh} =[30 60];
                time_window(fh) = time_window_options(dsearchn(Fourier.freq',HF{fh}(1)));
        end
    end
    
    
    
    rho = zeros(size(times2save,2),size(Fourier.fourierspctrm,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% left off here
    for fi = 1:length(Deci.Analysis.Connectivity.freq)
        
        cfg = [];
        cfg.frequency = HF{fi};
        
        cfg.channel = Fourier.label{chan_idx(1)};
        data1 = ft_selectdata(cfg,Fourier);
        data1_temp_powspctrm = squeeze(data1.fourierspctrm);
        data1_temp_powspctrm = abs(data1_temp_powspctrm).^2 ;
        data1.FreqBand_Avg_powspctrm = squeeze(mean(data1_temp_powspctrm,2));
        clear data1_temp_powspctrm
        data1 = rmfield(data1,'fourierspctrm');
        data1 = rmfield(data1,'cumtapcnt');
        data1 = rmfield(data1,'time');
        data1 = rmfield(data1,'trialinfo');
        data1 = rmfield(data1,'condinfo');
        data1 = rmfield(data1,'cfg');
        
        cfg.channel = Fourier.label{chan_idx(2)};
        data2 = ft_selectdata(cfg,Fourier);
        data2_temp_powspctrm = squeeze(data2.fourierspctrm);
        data2_temp_powspctrm = abs(data2_temp_powspctrm).^2 ;
        data2.FreqBand_Avg_powspctrm = squeeze(mean(data2_temp_powspctrm,2));
        clear data2_temp_powspctrm
        data2 = rmfield(data2,'fourierspctrm');
        data2 = rmfield(data2,'cumtapcnt');
        data2 = rmfield(data2,'time');
        data2 = rmfield(data2,'trialinfo');
        data2 = rmfield(data2,'condinfo');
        data2 = rmfield(data2,'cfg');
        

        %compute time window in indicies for this freq
        time_window_idx = round((1000/HF{fi}(1))*time_window(fi)/(1000/Deci.Analysis.DownSample));
        
        for ti = 1:length(times2save)
            
            for trl = 1:size(data1.FreqBand_Avg_powspctrm,1)
                
                %compute Spearman correlation coeff
                rho(ti,trl) = corr(data1.FreqBand_Avg_powspctrm(trl,times2save_idx(ti)-time_window_idx:times2save_idx(ti)+time_window_idx)',data2.FreqBand_Avg_powspctrm(trl,times2save_idx(ti)-time_window_idx:times2save_idx(ti)+time_window_idx)','Type','Spearman');
            end
        end
        clear data1
        clear data2
        
        PowPow_corr_data.rho = rho;
        PowPow_corr_data.dimord = 'time_trial';
        PowPow_corr_data.chan1 = Deci.Analysis.Connectivity.chan1;
        PowPow_corr_data.chan2 = Deci.Analysis.Connectivity.chan2;
        PowPow_corr_data.time = times2save;
        PowPow_corr_data.freq = Fourier.freq;
        PowPow_corr_data.label = Fourier.label;
        
        if Deci.Analysis.Laplace
            save([Deci.Folder.Analysis filesep 'PowPow_TxT' filesep 'Laplacian' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep [Deci.Analysis.Connectivity.chan1{1} '_' Deci.Analysis.Connectivity.chan2{1}] filesep Deci.Analysis.Connectivity.freq{fi}],'PowPow_corr_data','-v7.3');
        else
            save([Deci.Folder.Analysis filesep 'PowPow_TxT' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep [Deci.Analysis.Connectivity.chan1{1} '_' Deci.Analysis.Connectivity.chan2{1}] filesep Deci.Analysis.Connectivity.freq{fi}],'PowPow_corr_data','-v7.3');
        end
        
        clear PowPow_corr_data rho 
        
    end
end
