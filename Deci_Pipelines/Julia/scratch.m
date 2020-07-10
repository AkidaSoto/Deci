freq.cfg.method = 'wavelet';
freq.cfg.output      = 'pow';
freq.cfg.keeptrials  = 'yes' ;
phase_freqs = params.phase_freqs;
freq.cfg.foilim        = [phase_freqs(1)-1 phase_freqs(1)+1];
freq.cfg.toi = 'all';

[freq2] = ft_freqanalysis(freq.cfg, dat);
[freq2] = ft_freqanalysis(cfg.cfg, cfg.dat);
filt = squeeze(sum(freq2.powspctrm,3));

        [filt] = ft_preproc_bandpassfilter(data, Fs,...
            [phase_freqs(phase)-1 phase_freqs(phase)+1],...
            filt_order, 'but', 'twopass', 'no');
[wavelet,cycles,freqresol,timeresol] = dftfilt3([phase_freqs(phase)-1:phase_freqs(phase)+1], 3, Fs, 'winsize', size(data,2));
ws = wavelet_spectrogram(data, Fs, [phase_freqs(phase)-1:phase_freqs(phase)+1], 3, 0, '');


for t = 0:1
    params.trialflag = t;
    for c = 1:3
        for d = 1:3
            if d ~= c
                params.PAC_cond1 = c;
                params.PAC_cond2 = d;
                [stat] = Deci_get_PAC_stats(Deci,params);
            end
        end
    end
end