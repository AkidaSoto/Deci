function stat = dc_connectivitystatistics(cfg,varargin)


%% Design Matrix is computed
% resampling design for permuting is also computed


%% Data Management
% All Varargins are concatenated into subject-wise matrix dat


dimord = varargin{1}.dimord;
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


if size(cfg.design,2)~=size(dat,2)
    cfg.design = transpose(cfg.design);
end

design = cfg.design;


if ~strcmp(cfg.method,'montecarlo')
    %% If simple Stat
    
    % check if fieldtrip is flexible
    [stat, cfg] = ft_statistics_analytic(cfg, dat, design);
    
    
else
    cfg.resampling   = ft_getopt(cfg, 'resampling', 'permutation');
    resampled_design = dc_resampledesign(cfg, design);
    Nrand = size(resampled_design,1);
    
    % fetch function handle to the low-level statistics function
    statfun = str2func(['ft_statfun_' cfg.statistic]);
    if isempty(statfun)
        error('could not locate the appropriate statistics function');
    else
        fprintf('using "%s" for the single-sample statistics\n', func2str(statfun));
    end

    % set the defaults for clustering
    cfg.clusterstatistic = ft_getopt(cfg, 'clusterstatistic', 'maxsum');
    cfg.clusterthreshold = ft_getopt(cfg, 'clusterthreshold', 'parametric');
    cfg.clusteralpha     = ft_getopt(cfg, 'clusteralpha',     0.05);
    cfg.clustercritval   = ft_getopt(cfg, 'clustercritval',   []);
    cfg.clustertail      = ft_getopt(cfg, 'clustertail',      cfg.tail);
    cfg.connectivity     = ft_getopt(cfg, 'connectivity',     []); % the default is dealt with below
    cfg.ivar         = ft_getopt(cfg, 'ivar',       'all');
    cfg.uvar         = ft_getopt(cfg, 'uvar',       []);
    cfg.cvar         = ft_getopt(cfg, 'cvar',       []);
    cfg.wvar         = ft_getopt(cfg, 'wvar',       []);
    cfg.feedback     = ft_getopt(cfg, 'feedback',   'text');
    cfg.precondition = ft_getopt(cfg, 'precondition', []);

    % determine the critical value for cluster thresholding
    
    fprintf('computing a parametric threshold for clustering\n');
    tmpcfg = [];
    tmpcfg.dimord         = cfg.dimord;
    tmpcfg.alpha          = cfg.clusteralpha;
    tmpcfg.tail           = cfg.clustertail;
    tmpcfg.ivar           = cfg.ivar;
    tmpcfg.uvar           = cfg.uvar;
    tmpcfg.cvar           = cfg.cvar;
    tmpcfg.wvar           = cfg.wvar;
    if isfield(cfg, 'contrastcoefs'), tmpcfg.contrastcoefs = cfg.contrastcoefs; end % needed for Erics F-test statfun
    tmpcfg.computecritval = 'yes';  % explicitly request the computation of the crtitical value
    tmpcfg.computestat    = 'no';   % skip the computation of the statistic
    
    try
        cfg.clustercritval    = getfield(statfun(tmpcfg, dat, design), 'critval');
    catch
        disp(lasterr);
        error('could not determine the parametric critical value for clustering');
    end

    % compute the statistic for the observed data
    ft_progress('init', cfg.feedback, 'computing statistic');
    % get an estimate of the time required per evaluation of the statfun
    time_pre = cputime;
    
    % both the statistic and the (updated) configuration are returned
    [statobs, cfg] = statfun(cfg, dat, design);

    
    if isstruct(statobs)
        % remember all details for later reference, continue to work with the statistic
        statfull = statobs;
        statobs  = getfield(statfull, 'stat');
    else
        % remember the statistic for later reference, continue to work with the statistic
        statfull.stat = statobs;
    end
    
    time_eval = cputime - time_pre;
    fprintf('estimated time per randomization is %.2f seconds\n', time_eval);
    
    % pre-allocate some memory
    statrand = zeros(size(statobs,1), size(resampled_design,1));

    
    % compute the statistic for the randomized data and count the outliers
    for i=1:Nrand
        ft_progress(i/Nrand, 'computing statistic %d from %d\n', i, Nrand);
        
        tmpdesign = design(:,resampled_design(i,:));     % the columns in the design matrix are reshufled by means of permutation
        tmpdat    = dat;                        % the data itself is not shuffled
        if size(tmpdesign,1)==size(tmpdat,2)
            tmpdesign = transpose(tmpdesign);
        end
        
        % keep each randomization in memory for cluster postprocessing
        dum = statfun(cfg, tmpdat, tmpdesign);
        if isstruct(dum)
            statrand(:,i) = dum.stat;
        else
            statrand(:,i) = dum;
        end
        
    end
    
    ft_progress('close');

    %% clustering stat
    
        cfg.neighbours   = ft_getopt(cfg, 'neighbours', []);
    
    tempcfg = cfg;
    tempcfg.channel =unique(varargin{1}.labelcmb(:,1));
    connectivitylow = channelconnectivity(tempcfg);
    
    if ~isfield(cfg, 'inside')
        cfg.inside = true(cfg.dim);
    end % cfg.inside is set in ft_sourcestatistics, but is also needed for timelock and freq
    
    tempcfg = cfg;
    tempcfg.channel =unique(varargin{1}.labelcmb(:,2));
    connectivityhigh = channelconnectivity(tempcfg);
    
    cfg.connectivity = connectivitylow;
    [stat, cfg] = clusterstat(cfg, statrand, statobs);
    %%
    
    if strcmp(cfg.correcttail, 'prob') && cfg.tail==0
        stat.prob = stat.prob .* 2;
        stat.prob(stat.prob>1) = 1; % clip at p=1
        % also correct the probabilities in the pos/negcluster fields
        if isfield(stat, 'posclusters')
            for i=1:length(stat.posclusters)
                stat.posclusters(i).prob = stat.posclusters(i).prob*2;
                if stat.posclusters(i).prob>1; stat.posclusters(i).prob = 1; end
            end
        end
        if isfield(stat, 'negclusters')
            for i=1:length(stat.negclusters)
                stat.negclusters(i).prob = stat.negclusters(i).prob*2;
                if stat.negclusters(i).prob>1; stat.negclusters(i).prob = 1; end
            end
        end
    elseif strcmp(cfg.correcttail, 'alpha') && cfg.tail==0
        cfg.alpha = cfg.alpha / 2;
    end
    
    % compute range of confidence interval p ? 1.96(sqrt(var(p))), with var(p) = var(x/n) = p*(1-p)/N
    stddev = sqrt(stat.prob.*(1-stat.prob)/Nrand);
    stat.cirange = 1.96*stddev;
    
    if isfield(stat, 'posclusters')
        for i=1:length(stat.posclusters)
            stat.posclusters(i).stddev  = sqrt(stat.posclusters(i).prob.*(1-stat.posclusters(i).prob)/Nrand);
            stat.posclusters(i).cirange =  1.96*stat.posclusters(i).stddev;
            if i==1 && stat.posclusters(i).prob<cfg.alpha && stat.posclusters(i).prob+stat.posclusters(i).cirange>=cfg.alpha
                warning('FieldTrip:posCluster_exceeds_alpha', sprintf('The p-value confidence interval of positive cluster #%i includes %.3f - consider increasing the number of permutations!', i, cfg.alpha));
            end
        end
    end
    if isfield(stat, 'negclusters')
        for i=1:length(stat.negclusters)
            stat.negclusters(i).stddev  = sqrt(stat.negclusters(i).prob.*(1-stat.negclusters(i).prob)/Nrand);
            stat.negclusters(i).cirange =  1.96*stat.negclusters(i).stddev;
            if i==1 && stat.negclusters(i).prob<cfg.alpha && stat.negclusters(i).prob+stat.negclusters(i).cirange>=cfg.alpha
                warning('FieldTrip:negCluster_exceeds_alpha', sprintf('The p-value confidence interval of negative cluster #%i includes %.3f - consider increasing the number of permutations!', i, cfg.alpha));
            end
        end
    end
    
    
end


% the statistical output contains multiple elements, e.g. F-value, beta-weights and probability
fn = fieldnames(stat);

for i=1:length(fn)
    if numel(stat.(fn{i}))==prod(datsiz)
        % reformat into the same dimensions as the input data
        stat.(fn{i}) = reshape(stat.(fn{i}), [datsiz 1]);
    end
end

end