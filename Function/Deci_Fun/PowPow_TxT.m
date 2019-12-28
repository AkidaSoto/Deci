function PowPow_TxT(Deci,info,Fourier,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by CTGill
% last edited 12/27/19
%
% PowPowTxT.m computes the power time series for two electrodes, averages
% power across specified frequency bins, and then computes the correlation
% coefficient between time varying power of the two electrodes in sliding
% time segments within a trial, and then for each trial.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ismember(Deci.Analysis.CondTitle(info.Cond),Deci.Analysis.Connectivity.cond)
    
    mkdir([Deci.Folder.Analysis filesep 'PowPow_TxT' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
    
    times2save = -.5:.02:1.5;
    time_window = linspace(1.5,3.5,length(Deci.Analysis.Freq.foi));
    
    %time in indicies
    times2save_idx = dsearchn(Fourier.time', times2save');
    
    
    chan_idx = zeros(1,2);
    chan_idx(1) = find(strcmpi(Deci.Analysis.Connectivity.chan1,Fourier.label));
    chan_idx(2) = find(strcmpi(Deci.Analysis.Connectivity.chan2,Fourier.label));
    
    
    %     %average power across each specified frequency band
    %     for fh = 1:length(Deci.Analysis.Connectivity.freq)
    %         switch Deci.Analysis.Connectivity.freq{fh}
    %             case 'theta'
    %                 HF{fh} = [4 8];
    %             case 'beta'
    %                 HF{fh} = [12.5 30];
    %             case 'alpha'
    %                 HF{fh} =[8 12.5];
    %             case 'lowgamma'
    %                 HF{fh} =[30 60];
    %         end
    %     end
    %
    %     for f =  1:length(Deci.Analysis.Connectivity.freq)
    cfg = [];
    %         cfg.frequency = [HF{f}];
    %
    cfg.channel = Fourier.label{chan_idx(1)};
    data1 = ft_selectdata(cfg,Fourier);
    data1_temp_powspctrm = squeeze(data1.fourierspctrm);
    data1.powspctrm = abs(data1_temp_powspctrm).^2 ;
    data1 = rmfield(data1,'fourierspctrm');
    
    cfg.channel = Fourier.label{chan_idx(2)};
    data2 = ft_selectdata(cfg,Fourier);
    data2_temp_powspctrm = squeeze(data2.fourierspctrm);
    data2.powspctrm = abs(data2_temp_powspctrm).^2 ;
    data2 = rmfield(data2,'fourierspctrm');
    
    rho = zeros(size(Fourier.freq,2),size(times2save,2),size(data1.powspctrm,1));
    
    for fi = 1:length(Fourier.freq)
        
        %compute time window in indicies for this freq
        time_window_idx = round((1000/Fourier.freq(fi))*time_window(fi)/(1000/Deci.Analysis.DownSample));
        
        for ti = 1:length(times2save)
            
            for trl = 1:size(data1.powspctrm,1)
                
                %compute Spearman correlation coeff
                rho(fi,ti,trl) = corr(squeeze(data1.powspctrm(trl,fi,times2save_idx(ti)-time_window_idx:times2save_idx(ti)+time_window_idx)),squeeze(data2.powspctrm(trl,fi,times2save_idx(ti)-time_window_idx:times2save_idx(ti)+time_window_idx)),'Type','Spearman');
            end
        end
    end
    
    
    PowPow_corr_data.rho = rho;
    PowPow_corr_data.dimord = 'freq_time_trial';
    PowPow_corr_data.chan1 = Deci.Analysis.Connectivity.chan1;
    PowPow_corr_data.chan2 = Deci.Analysis.Connectivity.chan2;
    PowPow_corr_data.time = times2save;
    PowPow_corr_data.freq = Fourier.freq;
    PowPow_corr_data.label = Fourier.label;
    
    save([Deci.Folder.Analysis filesep 'PowPow_TxT' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep [Deci.Analysis.Connectivity.chan1{1} '_' Deci.Analysis.Connectivity.chan2{1}]],'PowPow_corr_data','-v7.3');
    
end
%end