function dc_connectivitystatistics(cfg,varargin)

%% Data Management
% All Varargins are concatenated into subject-wise matrix dat

dimord = getdimord(varargin{1}, cfg.parameter);
dimtok = tokenize(dimord, '_');
dimsiz = getdimsiz(varargin{1}, cfg.parameter);
dimsiz(end+1:length(dimtok)) = 1; % there can be additional trailing singleton dimensions
rptdim = find( strcmp(dimtok, 'subj') |  strcmp(dimtok, 'rpt') |  strcmp(dimtok, 'rpttap'));
datdim = find(~strcmp(dimtok, 'subj') & ~strcmp(dimtok, 'rpt') & ~strcmp(dimtok, 'rpttap'));
datsiz = dimsiz(datdim);

% pass these fields to the low-level functions, they should be removed further down
cfg.dimord = sprintf('%s_', dimtok{datdim});
cfg.dimord = cfg.dimord(1:end-1); % remove trailing _
cfg.dim    = dimsiz(datdim);

if isempty(rptdim)
    % repetitions are across multiple inputs
    dat = nan(prod(dimsiz), length(varargin));
    for i=1:length(varargin)
        tmp = varargin{i}.(cfg.parameter);
        dat(:,i) = tmp(:);
    end
else
    % repetitions are within inputs
    dat = cell(size(varargin));
    for i=1:length(varargin)
        tmp = varargin{i}.(cfg.parameter);
        if rptdim~=1
            % move the repetitions to the first dimension
            tmp = permute(tmp, [rptdim datdim]);
        end
        dat{i} = reshape(tmp, size(tmp,1), []);
    end
    dat = cat(1, dat{:});   % repetitions along 1st dimension
    dat = dat';             % repetitions along 2nd dimension
end

%% Design Matrix is computed
% resampling design for permuting is also computed

if size(cfg.design,2)~=size(dat,2)
    cfg.design = transpose(cfg.design);
end

design = cfg.design;

resampled_design = dc_resampledesign(cfg, design);

%% If simple Stat

if ~strcmp(cfg.method,'montecarlo')
    
    % check if fieldtrip is flexible
    [stat, cfg] = ft_statistics_analytic(cfg, dat, resampled_design);

else
    [stat, cfg] = ft_statistics_montecarlo(cfg, dat, resampled_design);
end
end