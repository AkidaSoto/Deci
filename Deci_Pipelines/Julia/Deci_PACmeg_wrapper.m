function Deci_PACmeg_wrapper(Deci,info,freq,dat,params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
num_trials = size(dat.trial,2);
time_len = size(cell2mat(dat.trial(1)),2);
data_table = zeros(num_trials,time_len);
for t = 1:num_trials
    data_table(t,:) = cell2mat(dat.trial(t));
end

% cfg.Fs            = Sampling frequency (in Hz)
% cfg.phase_freqs   = Phase Frequencies in Hz (e.g. [8:1:13])
% cfg.amp_freqs     = Amplitude Frequencies in Hz (e.g. [40:2:100])
% cfg.filt_order    = Filter order used by ft_preproc_bandpassfilter
% amp_bandw_method  = Method for calculating bandwidth
% cfg.method        = Method for PAC Computation
% cfg.surr_method   = Method to compute surrogates
% cfg.surr_N        = Number of iterations used for surrogate analysis
% cfg.mask          = filters ALL data but masks between times [a b]
% cfg.avg_PAC       = Average PAC over trials ('yes' or 'no')
params.Fs = 1/(dat.time{1,1}(2)-dat.time{1,1}(1));
[MI_matrix_raw,MI_matrix_surr] = PACmeg(params,data_table);
MI_full = freq;
MI_full.MI_matrix_raw = MI_matrix_raw;
MI_full.MI_matrix_surr = MI_matrix_surr;
MI_full.params = params;
MI_full.info = info;
%MI_full = rmfield(MI_full,'freq');
mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Corr' filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'Corr' filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep info.Channels{info.ChanNum}],'MI_full','-v7.3');
end

