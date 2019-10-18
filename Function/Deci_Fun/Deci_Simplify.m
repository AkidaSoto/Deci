function Deci_Simplify(Deci,info,freq,params)

Magnitude = abs(freq.fourierspctrm.^2);
Phase = angle(freq.fourierspctrm);

HF = [];

for fh = 1:length(params.freqbin)
    switch params.freqbin{fh}
        case 'theta'
            HF{fh} = [4 8];
        case 'highbeta'
            HF{fh} = [21 30];
        case 'lowbeta'
            HF{fh} = [13 20];
        case 'alpha'
            HF{fh} =[9 12];
        case 'lowgamma'
            HF{fh} =[30 55];
        case 'highgamma'
            HF{fh} = [55 80];
    end
end
params.freqbin = HF;

latency = reshape([freq.time [nan(-rem(length(freq.time),params.timebin)+params.timebin,1)]'],[ceil([length(freq.time)]/params.timebin) params.timebin]);
latency = [min(latency,[],1);max(latency,[],1)];

for fois = 1:length(params.freqbin)
    
    foi = freq.freq >= params.freqbin{fois}(1) & freq.freq <= params.freqbin{fois}(2);
    
    for tois = 1:params.timebin
        
        toi = freq.time >= latency(1,tois) & freq.time <= latency(2,tois);
        
        Mag(:,:,fois,tois) = squeeze(nanmean(nanmean(Magnitude(:,:,foi,toi),3),4));
        Pha(:,:,fois,tois) = squeeze(circ_mean(circ_mean(Phase(:,:,foi,toi),[],3),[],4));
        
    end
end

extra = freq;
extra.freq = mean(cat(1,params.freqbin{:}),2)';
extra.time = mean(latency,1);
extra = rmfield(extra,'fourierspctrm');
extra.powspctrm = Mag;
Mag = extra;

extra.powspctrm = Pha;
Pha = extra;

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Simp' filesep 'Mag'  filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'Simp' filesep 'Mag' filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep info.Channels{info.ChanNum}],'Mag');

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Simp' filesep 'Pha'  filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'Simp' filesep 'Pha' filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep info.Channels{info.ChanNum}],'Pha');

end