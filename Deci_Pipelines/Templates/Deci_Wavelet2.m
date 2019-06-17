Deci.Analysis.Freq.method        = 'wavelet';                                                  % Currently only uses 'wavelet' and 'hilbert'
Deci.Analysis.Freq.foi           = exp(linspace(log(1),log(50),50));                           % Frequency of Interest
Deci.Analysis.Freq.width         = 4;                                                         % Width
Deci.Analysis.Freq.gwidth        = 10;                                                          % Gwidth
Deci.Analysis.Freq.Toi           = [-.5 1.5];                                                  % Time Range to save
Deci.Analysis.Freq.Toilim        = [-2 3];