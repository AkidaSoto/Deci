Deci.Analysis.Freq.method        = 'hilbert';                                                  % Currently only uses Wavelet Decomp
Deci.Analysis.Freq.foi           = exp(linspace(log(4),log(80),40)) ;                           % Frequency of Interest
Deci.Analysis.Freq.width         = log(Deci.Analysis.Freq.foi) ;                                                         % Width
Deci.Analysis.Freq.transition_width = .2;
Deci.Analysis.Freq.order         = 3;
Deci.Analysis.Freq.Toi           = [-.5 1.5];                                                  % Time Range to save
Deci.Analysis.Freq.Toilim        = [-2 3];                                                     % Time Range to do freq analysis on