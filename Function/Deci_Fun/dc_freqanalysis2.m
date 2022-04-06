function [freq] = ft_freqanalysis2(cfg, data)

EEG.data = cell2mat(data.trial');
EEG.data = permute(EEG.data,[3 2 1]);

% Setup Wavelet Params
num_freqs=50;
frex=logspace(.01,1.7,num_freqs);
s=logspace(log10(3),log10(10),num_freqs)./(2*pi*frex);
t=data.time{1};

% Definte Convolution Parameters
dims = size(EEG.data);
n_wavelet = length(t);
n_data = dims(2)*dims(3);
n_convolution = n_wavelet+n_data-1;
n_conv_pow2 = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

tx=t;
T1=nearest(tx,min(cfg.toi));  T2=nearest(tx,max(cfg.toi));

% Get FFT of data

    EEG_fft = fft(reshape(EEG.data,1,n_data),n_conv_pow2);

    for fi=1:num_freqs

        wavelet = fft( exp(2*1i*pi*frex(fi).*t) .* exp(-t.^2./(2*(s(fi)^2))) , n_conv_pow2 );  % sqrt(1/(s(fi)*sqrt(pi))) *

        % convolution
        EEG_conv = ifft(wavelet.*EEG_fft);
        EEG_conv = EEG_conv(1:n_convolution);
        EEG_conv = EEG_conv(half_of_wavelet_size+1:end-half_of_wavelet_size);
        EEG_conv = reshape(EEG_conv,dims(2),dims(3));

        Fourier(1,fi,:,:) = EEG_conv(T1:T2,:);

        % Get power by condi (different time windows)
        %temp_POWER1 = mean(abs(EEG_conv(T1:T2,:)).^2,2);
        %temp_POWER2 = mean(abs(EEG_conv(T1:T2,:)).^2,2);

    end

    freq.fourierspctrm = permute(Fourier,[4 1 2 3]);
    freq.time = tx(T1:T2);
    freq.label = data.label;
    freq.dimord = 'rpttap_chan_freq_time';
    freq.trialinfo = data.trialinfo;
    freq.freq = frex;
end