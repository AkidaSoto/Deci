Deci.Analysis.Freq.method        = 'wavelet';                                                  % Currently only uses 'wavelet' and 'hilbert'
Deci.Analysis.Freq.foi           = exp(linspace(log(2),log(40),30));                           % Frequency of Interest
Deci.Analysis.Freq.width         = exp(linspace(log(3),log(13),30));                                                         % Width
Deci.Analysis.Freq.gwidth        = 1;                                                          % Gwidth
Deci.Analysis.Toi           = [-.5 1.5];                                                  % Time Range to save
Deci.Analysis.Toilim        = [-2 3];

%pls use dc_checkmyfreqs(exp(linspace(log(1),log(60),30)),[2500],500)
% to check validity of the freq replace [2500] with the number of samples 
% remember that default is 500samples/sec (last parameter) with downsampling. 