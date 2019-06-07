Deci.Analysis.Freq.method        = 'wavelet';                                                  % Currently only uses 'wavelet' and 'hilbert'
Deci.Analysis.Freq.foi           = exp(linspace(log(4),log(40),20));                           % Frequency of Interest
Deci.Analysis.Freq.width         = 7 ;                                                         % Width
Deci.Analysis.Freq.gwidth        = 3;                                                          % Gwidth
Deci.Analysis.Freq.Toi           = [-.5 1.5];                                                  % Time Range to save
Deci.Analysis.Freq.Toilim        = [-2 3];