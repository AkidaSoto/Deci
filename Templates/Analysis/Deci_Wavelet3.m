Deci.Analysis.Freq.method        = 'wavelet';                                                  % Currently only uses 'wavelet' and 'hilbert'
Deci.Analysis.Freq.foi           = 4:.5:8;                           % Frequency of Interest
Deci.Analysis.Freq.width         = 5;                                                         % Width
Deci.Analysis.Freq.gwidth        = 7;                                                          % Gwidth
Deci.Analysis.Toi           = [-.5 1.5];                                                  % Time Range to save
Deci.Analysis.Toilim        = [-2 3];

%pls use dc_checkmyfreqs(exp(linspace(log(1),log(60),30)),[2500],500)
% to check validity of the freq replace [2500] with the number of samples 
% remember that default is 500samples/sec (last parameter) with downsampling. 