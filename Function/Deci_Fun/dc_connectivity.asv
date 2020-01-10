function dc_connectivity(Deci,info,Fourier,params)

for Conn = 1:length(params.List)
    
    Current = params.List{Conn};
    chanl = Current{1};
    chanh = Current{2};
    freqlow = Current{3};
    freqhigh = Current{4};
    conne = Current{5};
    
    chanlow = find(ismember(Fourier.label,chanl));
    chanhigh = find(ismember(Fourier.label,chanh));
    
    
    switch freqlow
        case 'theta'
            LF = [4 8];
        case 'lowbeta'
            LF = [12.5 21];
        case 'highbeta'
            LF = [21 30];
        case 'alpha'
            LF =[8 12.5];
        case 'lowgamma'
            LF =[30 55];
        case 'highgamma'
            LF = [55 80];
    end
    
    switch freqhigh
        case 'theta'
            HF= [4 8];
        case 'lowbeta'
            HF = [12.5 21];
        case 'highbeta'
            HF = [21 30];
        case 'alpha'
            HF =[8 12.5];
        case 'lowgamma'
            HF =[30 55];
        case 'highgamma'
            HF = [55 80];
    end
    
    %lcfg.latency = params.toi;
    lcfg.frequency = LF;
    lcfg.channel = chanl;
    evalc('datalow = ft_selectdata(lcfg,Fourier)');
    freqs = datalow.freq;
    
    time_window = params.window; %linspace(1.5,3.5,length(Deci.Analysis.Freq.foi));
    time_window = time_window(ismember(Fourier.freq,datalow.freq));
    
    cfg.resamplefs = params.DownSample;
    evalc('datalow = ft_resampledata(cfg,datalow)');
    datalow.fourierspctrm = permute(cell2mat(permute(datalow.trial,[3 1 2])),[3 4 1 2]);
    datalow.freq = freqs;
    datalow.time = datalow.time{1};
    datalow.dimord = 'rpt_chan_freq_time';
    
    %hcfg.latency = params.toi;
    hcfg.frequency = HF;
    hcfg.channel = chanh;
    evalc('datahigh = ft_selectdata(hcfg,Fourier)');
    freqs = datahigh.freq;
    
    cfg.resamplefs = params.DownSample;
    evalc('datahigh = ft_resampledata(cfg,datalow)');
    datahigh.fourierspctrm = permute(cell2mat(permute(datahigh.trial,[3 1 2])),[3 4 1 2]);
    datahigh.freq = freqs;
    datahigh.time = datahigh.time{1};
    datahigh.dimord = 'rpt_chan_freq_time';
    
    toi = find(datalow.time >= round(params.toi(1),4) & datalow.time <= round(params.toi(2),4));
    
    switch conne
        case 'ispc'
            
            if chanlow == chanhigh
                continue
            end
            
            for foi = 1:size(datalow.fourierspctrm,3)
                
                phase_low = angle(datalow.fourierspctrm(:,:,foi,:));
                phase_high = angle(datahigh.fourierspctrm(:,:,foi,:));
                %display('ispc only uses freqlow')
                
                %phase angle differences
                phase_angle_diffs = phase_low - phase_high;
                
                %compute time window in indicies for this freq
                time_window_idx = round((1000/datalow.freq(foi))*time_window(foi)/(1000*mean(diff(datalow.time))));
                
                for ti = 1:length(toi)
                    %compute phase snychronization
                    ispc(:,foi,ti) = abs(mean(exp(1i*phase_angle_diffs(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx)),4));
                end
            end
            
            conn.param(:,1,:) = nanmean(ispc,2);
            clear ispc phase_angle_diffs phase_low phase_high
            
            conn.dimord = 'rpt_freq_time';
            conn.chanlow = Fourier.label(chanlow);
            conn.chanhigh = Fourier.label(chanhigh);
            conn.time = datalow.time(toi);
            conn.freqlow = mean(LF,1);
            
            if params.rmvtrls
                conn.dimord = 'freq_time';
                conn.param = permute(nanmean(conn.param,1),[2 3 1]);
            end
            
            conn.lockers = info.lockers;
            conn.trllen = info.trllen;
            
            mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
            save([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep strjoin(Current,'_')],'conn','-v7.3');
            clear conn
            
            
        case 'plv'
            
            if isequal(HF,LF) && chanlow == chanhigh
                continue;
            else
                
                for foi_low = 1:size(datalow.fourierspctrm,3)
                    time_window_idx = round((1000/datalow.freq(foi_low))*time_window(foi_low)/(1000*mean(diff(datalow.time))));
                    
                    for foi_high = 1:size(datahigh.fourierspctrm,3)
                        
                        phase_low = angle(Fourier.fourierspctrm(:,chanlow,foi_low,:));
                        phase_high    = abs(Fourier.fourierspctrm(:,chanhigh,foi_high,:));
                        phase_high(isnan(phase_high(:))) = 0;
                        phase_high   = angle(hilbert(phase_high));
                        
                        phase_angle_diffs = phase_low - phase_high;
                        
                        for ti = 1:length(toi)
                            plv(:,foi_low,foi_high,ti) = abs(mean(exp(1i*phase_angle_diffs(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx)),4));
                        end
                    end
                end
            end
            
            conn.param(:,1,1,:) = permute(nanmean(nanmean(plv,2),3),[1 4 2 3]);
            clear plv phase_low phase_high phase_angle_diffs
            
            conn.dimord = 'rpt_freqlow_freqhigh_time';
            conn.chanlow = Fourier.label(chanlow);
            conn.chanhigh = Fourier.label(chanhigh);
            conn.time = datalow.time(toi);
            conn.freqlow = mean(LF,1);
            conn.freqhigh = mean(HF,1);
            if params.rmvtrls
                conn.dimord = 'freqlow_freqhigh_time';
                conn.param = permute(nanmean(conn.param,1),[2 3 4 1]);
            end
            
            conn.lockers = info.lockers;
            conn.trllen = info.trllen;
            mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
            save([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep strjoin(Current,'_')],'conn','-v7.3');
            clear conn
            
        case 'mvl'
            
            for foi_low = 1:size(datalow.fourierspctrm,3)
                time_window_idx = round((1000/datalow.freq(foi_low))*time_window(foi_low)/(1000*mean(diff(datalow.time))));
                
                for foi_high = 1:size(datahigh.fourierspctrm,3)
                    
                    phaselow = angle(Fourier.fourierspctrm(:,chanlow,foi_low,:));
                    amphigh    = abs(Fourier.fourierspctrm(:,chanhigh,foi_high,:));
                    
                    for ti = 1:length(toi)
                        
                        mvl(:,foi_low,foi_high,ti) = abs(mean(amphigh(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx).*exp(1i*phaselow(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx)),4));
                        
                    end
                end
            end
            
            conn.param(:,1,1,:) = permute(nanmean(nanmean(mvl,2),3),[1 4 2 3]);
            clear mvl phaselow amphigh
            
            conn.dimord = 'rpt_freqlow_freqhigh_time';
            conn.chanlow = Fourier.label(chanlow);
            conn.chanhigh = Fourier.label(chanhigh);
            conn.time = datalow.time(toi);
            conn.freqlow = mean(LF,1);
            conn.freqhigh = mean(HF,1);
            if params.rmvtrls
                conn.dimord = 'freqlow_freqhigh_time';
                conn.param = permute(nanmean(conn.param,1),[2 3 4 1]);
            end
            
            conn.lockers = info.lockers;
            conn.trllen = info.trllen;
            mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
            save([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep strjoin(Current,'_')],'conn','-v7.3');
            clear conn
            
        case 'pac'
            
            nbin = 21;
            
            for foi_low = 1:size(datalow.fourierspctrm,3)
                time_window_idx = round((1000/datalow.freq(foi_low))*time_window(foi_low)/(1000*mean(diff(datalow.time))));
                
                for foi_high = 1:size(datahigh.fourierspctrm,3)
                    
                    phaselow = angle(Fourier.fourierspctrm(:,chanlow,foi_low,:));
                    amphigh    = abs(Fourier.fourierspctrm(:,chanhigh,foi_high,:));
                    
                    for ti = 1:length(toi)
                        
                        [~,bin] = histc(phaselow(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx), linspace(-pi,pi,nbin));  % binned low frequency phase
                        binamp = zeros(size(amphigh(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx),1),nbin);      % binned amplitude
                        
                        for k = 1:nbin-1
                            idx = bin == k ;
                            pacdata(k) = squeeze(mean(mean(amphigh(idx),4),1));
                        end
                        
                        Q =ones(nbin-1,1)/[nbin-1];
                        P = pacdata/ nansum(pacdata);
                        
                        pac(foi_low,foi_high,ti) = nansum(P.* log2(P./Q'))./log2(nbin-1);
                    end
                end
            end
            conn.param(1,1,:) = permute(nanmean(nanmean(pac,1),2),[3 2 1]);
            clear pac phaselow amphigh
            
            conn.dimord = 'freqlow_freqhigh_time';
            conn.chanlow = Fourier.label(chanlow);
            conn.chanhigh = Fourier.label(chanhigh);
            conn.time = datalow.time(toi);
            conn.freqlow = mean(LF,1);
            conn.freqhigh = mean(HF,1);
            
            conn.lockers = info.lockers;
            conn.trllen = info.trllen;
            mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
            save([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep strjoin(Current,'_')],'conn','-v7.3');
            clear conn
            
        case 'cs_cl'
            
            for foi_low = 1:size(datalow.fourierspctrm,3)
                time_window_idx = round((1000/datalow.freq(foi_low))*time_window(foi_low)/(1000*mean(diff(datalow.time))));
                
                for foi_high = 1:size(datahigh.fourierspctrm,3)
                    
                    phaselow = angle(Fourier.fourierspctrm(:,chanlow,foi_low,:));
                    amphigh    = abs(Fourier.fourierspctrm(:,chanhigh,foi_high,:));
                    
                    
                    for ti = 1:length(toi)
                        
                        pha = circ_ang2rad(phaselow(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx));
                        amp =  amphigh(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx);
                        cs_cl(foi_low,foi_high,ti) = circ_corrcl(pha(:),amp(:));
                        
                    end
                end
            end
            
            conn.param(1,1,:) = permute(nanmean(nanmean(cs_cl,1),2),[3 1 2]);
            clear cs_cl phaselow amphigh
            
            conn.dimord = 'freqlow_freqhigh_time';
            conn.chanlow = Fourier.label(chanlow);
            conn.chanhigh = Fourier.label(chanhigh);
            conn.time = datalow.time(toi);
            conn.freqlow = mean(LF,1);
            conn.freqhigh = mean(HF,1);
            
            conn.lockers = info.lockers;
            conn.trllen = info.trllen;
            mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
            save([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep strjoin(Current,'_')],'conn','-v7.3');
            clear conn
        case 'cs_cc'
            
            if isequal(HF,LF) && chanlow == chanhigh
                continue;
            else
                
                
                for foi_low = 1:size(datalow.fourierspctrm,3)
                    time_window_idx = round((1000/datalow.freq(foi_low))*time_window(foi_low)/(1000*mean(diff(datalow.time))));
                    
                    for foi_high = 1:size(datahigh.fourierspctrm,3)
                        
                        phaselow = angle(Fourier.fourierspctrm(:,chanlow,foi_low,:));
                        phasehigh = angle(Fourier.fourierspctrm(:,chanhigh,foi_high,:));
                        
                        for ti = 1:length(toi)
                            
                            phalow = circ_ang2rad(phaselow(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx));
                            phahigh = circ_ang2rad(phasehigh(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx));
                            
                            cs_cl(foi_low,foi_high,ti) = circ_corrcc(phalow(:),phahigh(:));
                            
                        end
                    end
                end
            end
            
            conn.param(1,1,:) = permute(nanmean(nanmean(cs_cl,1),2),[3 1 2]);
            clear cs_cl phaselow amphigh
            
            conn.dimord = 'rpt_freqlow_freqhigh';
            conn.chanlow = Fourier.label(chanlow);
            conn.chanhigh = Fourier.label(chanhigh);
            conn.time = datalow.time(toi);
            conn.freqlow = mean(LF,1);
            conn.freqhigh = mean(HF,1);
            
            conn.lockers = info.lockers;
            conn.trllen = info.trllen;
            mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
            save([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep strjoin(Current,'_')],'conn','-v7.3');
            clear conn
    end
    
    
    
    
end

end