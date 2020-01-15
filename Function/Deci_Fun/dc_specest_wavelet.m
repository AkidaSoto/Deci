function [spectrum,freqoi,timeoi] = ft_specest_wavelet(dat, time, varargin)

freqoi    = ft_getopt(varargin, 'freqoi', 'all');
timeoi    = ft_getopt(varargin, 'timeoi', 'all');
width     = ft_getopt(varargin, 'width', 7);
gwidth    = ft_getopt(varargin, 'gwidth', 3);
pad       = ft_getopt(varargin, 'pad');
padtype   = ft_getopt(varargin, 'padtype', 'zero');
polyorder = ft_getopt(varargin, 'polyorder', 0);

% Set n's
[nchan,ndatsample] = size(dat);

% Remove polynomial fit from the data -> default is demeaning
if polyorder >= 0
  dat = ft_preproc_polyremoval(dat, polyorder, 1, ndatsample);
end

% Determine fsample and set total time-length of data
fsample = 1./mean(diff(time));
dattime = ndatsample / fsample; % total time in seconds of input data

% Zero padding
if round(pad * fsample) < ndatsample
  error('the padding that you specified is shorter than the data');
end

if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end

postpad    = round((pad - dattime) * fsample);
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;            % total time in seconds of padded data

% Set freqboi and freqoi
freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1; % is equivalent to: round(freqoi .* endtime) + 1;
freqboi   = unique(freqboi);
freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output

nfreqoi  = length(freqoi);

% Set timeboi and timeoi
offset = round(time(1)*fsample);

timeoi   = unique(round(timeoi .* fsample) ./ fsample);
timeboi  = round(timeoi .* fsample - offset) + 1;
ntimeboi = length(timeboi);

% Compute fft
spectrum = complex(nan(nchan,nfreqoi,ntimeboi),nan(nchan,nfreqoi,ntimeboi));
datspectrum = fft(ft_preproc_padding(dat, padtype, 0, postpad), [], 2);

% Creating wavelets
% expand width to array if constant width
if numel(width) == 1
  width = ones(1,nfreqoi) * width;
end

for ifreqoi = 1:nfreqoi

  dt = 1/fsample;
  sf = freqoi(ifreqoi) / width(ifreqoi);
  st = 1/(2*pi*sf);
  toi2 = -gwidth*st:dt:gwidth*st;
  A = 1/sqrt(st*sqrt(pi));
  tap = (A*exp(-toi2.^2/(2*st^2)))';
  acttapnumsmp = size(tap,1);
  ins = ceil(endnsample./2) - floor(acttapnumsmp./2);
  prezer = zeros(ins,1);
  pstzer = zeros(endnsample - ((ins-1) + acttapnumsmp)-1,1);
  
  % produce angle with convention: cos must always be 1  and sin must always be centered in upgoing flank, so the centre of the wavelet (untapered) has angle = 0
  ind  = (-(acttapnumsmp-1)/2 : (acttapnumsmp-1)/2)'   .*  ((2.*pi./fsample) .* freqoi(ifreqoi));
  
  % create wavelet and fft it
  % compute indices that will be used to extracted the requested fft output
  reqtimeboiind    = find((timeboi >=  (acttapnumsmp ./ 2)) & (timeboi < (ndatsample - (acttapnumsmp ./2))));
  reqtimeboi       = timeboi(reqtimeboiind);
  
  % compute datspectrum*wavelet, if there are reqtimeboi's that have data
  if ~isempty(reqtimeboi)
          dum = [fftshift(ifft(datspectrum .* repmat(fft(complex(vertcat(prezer,tap.*cos(ind),pstzer), vertcat(prezer,tap.*sin(ind),pstzer)),[],1)',[nchan 1]), [], 2),2)] .* sqrt(2 ./ fsample);
          spectrum(:,ifreqoi,reqtimeboiind) = dum(:,reqtimeboi);
  end
end

