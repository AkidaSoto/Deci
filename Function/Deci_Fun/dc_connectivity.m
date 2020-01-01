function dc_connectivity(Deci,info,Fourier,params)

%% Freqs
for fl = 1:length(params.freqlow)
    switch params.freqlow{fl}
        case 'theta'
            LF(:,fl) = [4 8];
        case 'lowbeta'
            LF(:,fl) = [12.5 21];
        case 'highbeta'
            LF(:,fl) = [21 30];
        case 'alpha'
            LF(:,fl) =[8 12.5];
        case 'lowgamma'
            LF(:,fl) =[30 55];
        case 'highgamma'
            LF(:,fl) = [55 80];
    end
end

for fh = 1:length(params.freqhigh)
    switch params.freqhigh{fh}
        case 'theta'
            HF(:,fh) = [4 8];
        case 'lowbeta'
            HF(:,fh) = [12.5 21];
        case 'highbeta'
            HF(:,fh) = [21 30];
        case 'alpha'
            HF(:,fh) =[8 12.5];
        case 'lowgamma'
            HF(:,fh) =[30 55];
        case 'highgamma'
            HF(:,fh) = [55 80];
    end
end


%% Channels
chanlow = find(ismember(Fourier.label,params.chanlow));
chanhigh = find(ismember(Fourier.label,params.chanhigh));

%% Time
toi = find(Fourier.time >= round(params.toi(1),4) & Fourier.time <= round(params.toi(2),4));
time_window = params.window; %linspace(1.5,3.5,length(Deci.Analysis.Freq.foi));

%% fourier spectrum for Chan/Freq
for conne = 1:length(params.type)
    
    for cl = 1:length(chanlow)
        for ch = 1:length(chanhigh)
            
            switch params.type{conne}
                case 'ispc'
                    
                    if chanlow(cl) == chanhigh(ch)
                        continue
                    end
                    
                    for fi = 1:size(LF,2)
                        
                        floi = find(Fourier.freq >= round(LF(1,fi),4) & Fourier.freq <= round(LF(2,fi),4));
                        
                        for foi = 1:length(floi)
                            
                            phase_low = angle(Fourier.fourierspctrm(:,chanlow(cl),floi(foi),:));
                            phase_high = angle(Fourier.fourierspctrm(:,chanhigh(ch),floi(foi),:));
                            %display('ispc only uses freqlow')
                            
                            %phase angle differences
                            phase_angle_diffs = phase_low - phase_high;
                            
                            %compute time window in indicies for this freq
                            time_window_idx = round((1000/Fourier.freq(floi(foi)))*time_window(floi(foi))/(1000*mean(diff(Fourier.time))));
                            
                            for ti = 1:length(toi)
                                %compute phase snychronization
                                ispc(:,foi,ti) = abs(mean(exp(1i*phase_angle_diffs(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx)),4));
                            end
                        end
                        
                        conn.param(:,fi,:) = nanmean(ispc,2);
                        clear ispc phase_angle_diffs phase_low phase_high
                    end
                    
                    conn.dimord = 'rpt_freq_time';
                    conn.chanlow = Fourier.label(chanlow(cl));
                    conn.chanhigh = Fourier.label(chanhigh(ch));
                    conn.time = Fourier.time(toi);
                    conn.freqlow = mean(LF,1);
                    
                    if params.rmvtrls
                        conn.dimord = 'freq_time';
                        conn.param = permute(nanmean(conn.param,1),[2 3 1]);
                    end
                    
                    mkdir([Deci.Folder.Analysis filesep params.type{conne} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
                    save([Deci.Folder.Analysis filesep params.type{conne} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep [Fourier.label(chanlow(cl)) '_' Fourier.label(chanhigh(ch))]],'conn','-v7.3');
                    clear conn
                    
                    
                case 'plv'
                    
                    for fli = 1:size(LF,2)
                        floi = find(Fourier.freq >= round(LF(1,fli),4) & Fourier.freq <= round(LF(2,fli),4));
                        for fhi = 1:size(HF,2)
                            fhoi = find(Fourier.freq >= round(HF(1,fhi),4) & Fourier.freq <= round(HF(2,fhi),4));
                            if isequal(HF(:,fhi),LF(:,fli)) && chanlow(cl) == chanhigh(ch)
                                plv = nan([size(Fourier.fourierspctrm,1) 1 1 size(toi,2)]);
                            else
                                for foi_low = 1:length(floi)
                                    time_window_idx = round((1000/Fourier.freq(floi(foi_low)))*time_window(floi(foi_low))/(1000*mean(diff(Fourier.time))));
                                    for foi_high = 1:length(fhoi)
                                        
                                        phase_low = angle(Fourier.fourierspctrm(:,chanlow(cl),floi(foi_low),:));
                                        phase_high    = abs(Fourier.fourierspctrm(:,chanhigh(ch),fhoi(foi_high),:));
                                        phase_high(isnan(phase_high(:))) = 0;
                                        phase_high   = angle(hilbert(phase_high));
                                        
                                        phase_angle_diffs = phase_low - phase_high;
                                        
                                        for ti = 1:length(toi)
                                            plv(:,foi_low,foi_high,ti) = abs(mean(exp(1i*phase_angle_diffs(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx)),4));
                                        end
                                    end
                                end
                            end
                            
                            conn.param(:,fli,fhi,:) = permute(nanmean(nanmean(plv,2),3),[1 4 2 3]);
                            clear plv phase_low phase_high phase_angle_diffs
                        end
                    end
                    
                    conn.dimord = 'rpt_freqlow_freqhigh_time';
                    conn.chanlow = Fourier.label(chanlow(cl));
                    conn.chanhigh = Fourier.label(chanhigh(ch));
                    conn.time = Fourier.time(toi);
                    conn.freqlow = mean(LF,1);
                    conn.freqhigh = mean(HF,1);
                    if params.rmvtrls
                        conn.dimord = 'freqlow_freqhigh_time';
                        conn.param = permute(nanmean(conn.param,1),[2 3 4 1]);
                    end
                    
                    mkdir([Deci.Folder.Analysis filesep 'Extra' filesep params.type{conne} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
                    save([Deci.Folder.Analysis filesep 'Extra' filesep params.type{conne} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep [Fourier.label{chanlow(cl)} '_' Fourier.label{chanhigh(ch)}]],'conn','-v7.3');
                    clear conn
                    
                case 'mvl'
                    
                    
                    for fli = 1:size(LF,2)
                        floi = find(Fourier.freq >= round(LF(1,fli),4) & Fourier.freq <= round(LF(2,fli),4));
                        for fhi = 1:size(HF,2)
                            fhoi = find(Fourier.freq >= round(HF(1,fhi),4) & Fourier.freq <= round(HF(2,fhi),4));
                            %                             if isequal(HF(:,fhi),LF(:,fli)) && chanlow(cl) == chanhigh(ch)
                            %                                 mvl = nan([size(Fourier.fourierspctrm,1) 1 1 size(toi,2)]);
                            %                             else
                            for foi_low = 1:length(floi)
                                time_window_idx = round((1000/Fourier.freq(floi(foi_low)))*time_window(floi(foi_low))/(1000*mean(diff(Fourier.time))));
                                for foi_high = 1:length(fhoi)
                                    
                                    phaselow = angle(Fourier.fourierspctrm(:,chanlow(cl),floi(foi_low),:));
                                    amphigh    = abs(Fourier.fourierspctrm(:,chanhigh(ch),fhoi(foi_high),:));
                                    
                                    for ti = 1:length(toi)
                                        
                                        mvl(:,foi_low,foi_high,ti) = abs(mean(amphigh(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx).*exp(1i*phaselow(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx)),4));
                                        
                                    end
                                end
                            end
                            %                             end
                            conn.param(:,fli,fhi,:) = permute(nanmean(nanmean(mvl,2),3),[1 4 2 3]);
                            clear mvl phaselow amphigh
                        end
                    end
                    
                    conn.dimord = 'rpt_freqlow_freqhigh_time';
                    conn.chanlow = Fourier.label(chanlow(cl));
                    conn.chanhigh = Fourier.label(chanhigh(ch));
                    conn.time = Fourier.time(toi);
                    conn.freqlow = mean(LF,1);
                    conn.freqhigh = mean(HF,1);
                    if params.rmvtrls
                        conn.dimord = 'freqlow_freqhigh_time';
                        conn.param = permute(nanmean(conn.param,1),[2 3 4 1]);
                    end
                    
                    mkdir([Deci.Folder.Analysis filesep 'Extra' filesep params.type{conne} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
                    save([Deci.Folder.Analysis filesep 'Extra' filesep params.type{conne} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep [Fourier.label{chanlow(cl)} '_' Fourier.label{chanhigh(ch)}]],'conn','-v7.3');
                    clear conn
                    
                case 'pac'
                    
                    nbin = 21;
                    
                    
                    for fli = 1:size(LF,2)
                        floi = find(Fourier.freq >= round(LF(1,fli),4) & Fourier.freq <= round(LF(2,fli),4));
                        for fhi = 1:size(HF,2)
                            fhoi = find(Fourier.freq >= round(HF(1,fhi),4) & Fourier.freq <= round(HF(2,fhi),4));
                            %                             if isequal(HF(:,fhi),LF(:,fli)) && chanlow(cl) == chanhigh(ch)
                            %                                 pac = nan([1 1 size(toi,2)]);
                            %                             else
                            for foi_low = 1:length(floi)
                                time_window_idx = round((1000/Fourier.freq(floi(foi_low)))*time_window(floi(foi_low))/(1000*mean(diff(Fourier.time))));
                                for foi_high = 1:length(fhoi)
                                    
                                    phaselow = angle(Fourier.fourierspctrm(:,chanlow(cl),floi(foi_low),:));
                                    amphigh    = abs(Fourier.fourierspctrm(:,chanhigh(ch),fhoi(foi_high),:));
                                    
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
                            %                             end
                            
                            conn.param(fli,fhi,:) = permute(nanmean(nanmean(pac,1),2),[3 2 1]);
                            clear pac phaselow amphigh
                        end
                    end
                    
                    conn.dimord = 'freqlow_freqhigh_time';
                    conn.chanlow = Fourier.label(chanlow(cl));
                    conn.chanhigh = Fourier.label(chanhigh(ch));
                    conn.time = Fourier.time(toi);
                    conn.freqlow = mean(LF,1);
                    conn.freqhigh = mean(HF,1);
                    
                    mkdir([Deci.Folder.Analysis filesep 'Extra' filesep params.type{conne} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
                    save([Deci.Folder.Analysis filesep 'Extra' filesep params.type{conne} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep [Fourier.label{chanlow(cl)} '_' Fourier.label{chanhigh(ch)}]],'conn','-v7.3');
                    clear conn
                    
                case 'cs_cl'
                    
                    for fli = 1:size(LF,2)
                        floi = find(Fourier.freq >= round(LF(1,fli),4) & Fourier.freq <= round(LF(2,fli),4));
                        for fhi = 1:size(HF,2)
                            fhoi = find(Fourier.freq >= round(HF(1,fhi),4) & Fourier.freq <= round(HF(2,fhi),4));
                            %                             if isequal(HF(:,fhi),LF(:,fli)) && chanlow(cl) == chanhigh(ch)
                            %                                 cl = nan([size(Fourier.fourierspctrm,1) 1 1 size(toi,2)]);
                            %                             else
                            for foi_low = 1:length(floi)
                                time_window_idx = round((1000/Fourier.freq(floi(foi_low)))*time_window(floi(foi_low))/(1000*mean(diff(Fourier.time))));
                                for foi_high = 1:length(fhoi)
                                    
                                    phaselow = angle(Fourier.fourierspctrm(:,chanlow(cl),floi(foi_low),:));
                                    amphigh    = abs(Fourier.fourierspctrm(:,chanhigh(ch),fhoi(foi_high),:));
                                    
                                    
                                    for ti = 1:length(toi)
                                        
                                        pha = circ_ang2rad(phaselow(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx));
                                        amp =  amphigh(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx);
                                        cs_cl(fli,fhi,ti) = circ_corrcl(pha(:),amp(:));
                                        
                                    end
                                end
                            end
                            %                             end
                            conn.param(fli,fhi,:) = permute(nanmean(nanmean(cs_cl,1),2),[3 1 2]);
                            clear cs_cl phaselow amphigh
                        end
                    end
                    conn.dimord = 'rpt_freqlow_freqhigh';
                    conn.chanlow = Fourier.label(chanlow(cl));
                    conn.chanhigh = Fourier.label(chanhigh(ch));
                    conn.time = Fourier.time(toi);
                    conn.freqlow = Fourier.freq(floi);
                    
                    mkdir([Deci.Folder.Analysis filesep 'Extra' filesep params.type{conne} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
                    save([Deci.Folder.Analysis filesep 'Extra' filesep params.type{conne} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep [Fourier.label{chanlow(cl)} '_' Fourier.label{chanhigh(ch)}]],'conn','-v7.3');
                    clear conn
                case 'cs_cc'
                    
                    for fli = 1:size(LF,2)
                        floi = find(Fourier.freq >= round(LF(1,fli),4) & Fourier.freq <= round(LF(2,fli),4));
                        for fhi = 1:size(HF,2)
                            fhoi = find(Fourier.freq >= round(HF(1,fhi),4) & Fourier.freq <= round(HF(2,fhi),4));
                            if isequal(HF(:,fhi),LF(:,fli)) && chanlow(cl) == chanhigh(ch)
                                cl = nan([size(Fourier.fourierspctrm,1) 1 1 size(toi,2)]);
                            else
                                for foi_low = 1:length(floi)
                                    time_window_idx = round((1000/Fourier.freq(floi(foi_low)))*time_window(floi(foi_low))/(1000*mean(diff(Fourier.time))));
                                    for foi_high = 1:length(fhoi)
                                        
                                        phaselow = angle(Fourier.fourierspctrm(:,chanlow(cl),floi(foi_low),:));
                                        phasehigh = angle(Fourier.fourierspctrm(:,chanhigh(ch),fhoi(foi_high),:));
                                        
                                        
                                        for ti = 1:length(toi)
                                            
                                            phalow = circ_ang2rad(phaselow(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx));
                                            phahigh = circ_ang2rad(phasehigh(:,:,:,toi(ti)-time_window_idx:toi(ti)+time_window_idx));
                                            
                                            cs_cl(fli,fhi,ti) = circ_corrcc(phalow(:),phahigh(:));
                                            
                                        end
                                    end
                                end
                            end
                            conn.param(fli,fhi,:) = permute(nanmean(nanmean(cs_cl,1),2),[3 1 2]);
                            clear cs_cl phaselow amphigh
                        end
                    end
                    conn.dimord = 'rpt_freqlow_freqhigh';
                    conn.chanlow = Fourier.label(chanlow(cl));
                    conn.chanhigh = Fourier.label(chanhigh(ch));
                    conn.time = Fourier.time(toi);
                    conn.freqlow = Fourier.freq(floi);
                    
                    mkdir([Deci.Folder.Analysis filesep 'Extra' filesep params.type{conne} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}]);
                    save([Deci.Folder.Analysis filesep 'Extra' filesep params.type{conne} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep [Fourier.label{chanlow(cl)} '_' Fourier.label{chanhigh(ch)}]],'conn','-v7.3');
                    clear conn
            end
            
        end
    end
    
end

end