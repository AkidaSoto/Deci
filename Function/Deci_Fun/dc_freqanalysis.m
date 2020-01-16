function [freq] = ft_freqanalysis(cfg, data)

% ensure that the required options are present
cfg.trials      = ft_getopt(cfg, 'trials',     'all', 1);
cfg.channel     = ft_getopt(cfg, 'channel',    'all');

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw', 'raw+comp', 'mvar'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'channel', 'showcallinfo'});

data = ft_selectdata(tmpcfg, data);

% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

% some proper error handling
if isfield(data, 'trial') && numel(data.trial)==0
    ft_error('no trials were selected'); % this does not apply for MVAR data
end

if numel(data.label)==0
    ft_error('no channels were selected');
end

% switch over method and do some of the method specfic checks and defaulting
cfg.width  = ft_getopt(cfg, 'width',  7);
cfg.gwidth = ft_getopt(cfg, 'gwidth', 3);

% set all the defaults
cfg.pad       = ft_getopt(cfg, 'pad',       []);
if isempty(cfg.pad)
    ft_notice('Default cfg.pad=''maxperlen'' can run slowly. Consider using cfg.pad=''nextpow2'' for more efficient FFT computation.')
    cfg.pad = 'maxperlen';
end

cfg.padtype   = ft_getopt(cfg, 'padtype',   'zero');
cfg.channel   = ft_getopt(cfg, 'channel',   'all');
cfg.precision = ft_getopt(cfg, 'precision', 'double');
cfg.foi       = ft_getopt(cfg, 'foi',       []);

cfg.polyremoval      = ft_getopt(cfg, 'polyremoval', 0);

% determine trial characteristics
ntrials = numel(data.trial);
trllength = cellfun(@(c) size(c,2),data.trial);

if any(logical(diff(trllength)))
   error('dc_freq is assumed to work only for trials will all equal length') 
end

if strcmp(cfg.pad, 'maxperlen')
    padding = max(trllength);
    cfg.pad = padding/data.fsample;
elseif strcmp(cfg.pad, 'nextpow2')
    padding = 2^nextpow2(max(trllength));
    cfg.pad = padding/data.fsample;
else
    padding = cfg.pad*data.fsample;
    if padding<max(trllength)
        ft_error('the specified padding is too short');
    end
end

oldfoi = cfg.foi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main loop over trials, inside fourierspectra are obtained and transformed into the appropriate outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is done on trial basis to save memory


% determine the corresponding indices of all channels
chanind    = match_str(data.label, cfg.channel);
nchan      = size(chanind,1);


nfoi   = numel(cfg.foi);
ntoi = numel(cfg.toi);

dimord    = 'rpttap_chan_freq_time';

freqoi    = cfg.foi;
timeoi    = cfg.toi;
width     = cfg.width;
gwidth    = cfg.gwidth;

polyorder = cfg.polyremoval;
[nchan,ndatsample] = size(data.trial{1});
    
% Determine fsample and set total time-length of data
fsample = 1./mean(diff(data.time{1}(1:2)));

endnsample = round(cfg.pad * fsample);

% Set freqboi and freqoi
freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1; % is equivalent to: round(freqoi .* endtime) + 1;
freqboi   = unique(freqboi);
freqoi    = (freqboi-1) ./ cfg.pad; % boi - 1 because 0 Hz is included in fourier output
nfreqoi  = length(freqoi);

if numel(width) == 1
    width = ones(1,nfreqoi) * width;
end

% Set timeboi and timeoi
offset = round(data.time{1}(1)*fsample);
timeoi   = unique(round(timeoi .* fsample) ./ fsample);
timeboi  = round(timeoi .* fsample - offset) + 1;

if polyorder >= 0
    data.trial = cellfun(@(c) ft_preproc_polyremoval(c,polyorder,1,ndatsample), data.trial,'un',0);
end



%data.trial = permute(cell2mat(permute(data.trial,[3 1 2])),[3 1 2]);

cfg = Exist(cfg,'gpu',0);
cfg = Exist(cfg,'cpu',0);

if cfg.cpu
     fourierspctrm = cell([ntrials 1]);
     
    parfor itrial = 1:ntrials
       fourierspctrm{itrial} = complex(zeros(1,nchan,nfoi,ntoi,cfg.precision));
        
        % Compute fft
        if cfg.gpu
            datspectrum = gpuArray(data.trial{itrial});
        else
            datspectrum = data.trial{itrial};
        end
        datspectrum = fft(datspectrum, [], 2);
        
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
                fourierspctrm{itrial}(1,:,ifreqoi,reqtimeboiind) = gather(dum(:,reqtimeboi));
            end
        end
    end
    
    fourierspctrm = cell2mat(fourierspctrm);
else
    
    fourierspctrm = complex(zeros(ntrials,nchan,nfoi,ntoi,cfg.precision));
    
    for itrial = 1:ntrials
        
        % Compute fft
        if cfg.gpu
            datspectrum = gpuArray(data.trial{itrial});
        else
            datspectrum = data.trial{itrial};
        end
        datspectrum = fft(datspectrum, [], 2);
        
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
                fourierspctrm(itrial,:,ifreqoi,reqtimeboiind) = gather(dum(:,reqtimeboi));
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END: Main loop over trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set output variables
freq        = [];
freq.label  = data.label;
freq.dimord = dimord;
freq.freq   = freqoi;
hasdc       = find(freqoi==0);
hasnyq      = find(freqoi==data.fsample./2);
hasdc_nyq   = [hasdc hasnyq];
freq.time = timeoi;
freq.cumtapcnt = ones([ntrials nfoi]);

% correct the 0 Hz or Nyqist bin if present, scaling with a factor of 2 is only appropriate for ~0 Hz
if ~isempty(hasdc_nyq)
    fourierspctrm(:,:,hasdc_nyq,:) = fourierspctrm(:,:,hasdc_nyq,:)./sqrt(2);
end
freq.fourierspctrm = fourierspctrm;

% some fields from the input should always be copied over in the output
freq = copyfields(data, freq, {'grad', 'elec', 'opto', 'topo', 'topolabel', 'unmixing'});

if isfield(data, 'trialinfo') && strcmp(cfg.keeptrials, 'yes')
    % copy the trialinfo into the output, but not the sampleinfo
    freq.trialinfo = data.trialinfo;
end

if length(oldfoi) == length(freq.freq)
    freq.freq = oldfoi;
end
