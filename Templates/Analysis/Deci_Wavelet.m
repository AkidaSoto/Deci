
Deci.Analysis.Freq.method        = 'wavelet';                                                  % Currently only uses 'wavelet' and 'hilbert'
Deci.Analysis.Freq.foi           = exp(linspace(log(2),log(40),40));                           % Frequency of Interest
Deci.Analysis.Freq.width         = exp(linspace(log(3),log(13),40));                                                         % Width
Deci.Analysis.Freq.gwidth        = 4;                                                          % Gwidth
Deci.Analysis.Freq.Toi           = [-.5 1.5];                                                  % Time Range to save
Deci.Analysis.Freq.Toilim        = [-2 3];

Deci.Analysis.ERP.Toi = Deci.Analysis.Freq.Toi;