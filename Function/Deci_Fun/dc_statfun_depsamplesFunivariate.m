function [s, cfg] = ft_statfun_depsamplesFmultivariate(cfg, dat, design)

% FT_STATFUN_DEPSAMPLESFMULTIVARIATE calculates the MANOVA dependent samples 
% F-statistic on the biological data in dat (the dependent variable), using 
% the information on the independent variable (ivar) in design.
%
% Use this function by calling one of the high-level statistics functions as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option
%   cfg.statistic = 'ft_statfun_depsamplesFmultivariate'
% see FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = ft_statfun_depsamplesFmultivariate(cfg, dat, design);
% where
%   dat    contains the biological data, Nsamples x Nreplications
%   design contains the independent variable (ivar) and the unit-of-observation (uvar) 
%          factor, Nfac x Nreplications
%
% Configuration options
%   cfg.contrastcoefs  = matrix of contrast coefficients determining the
%                        effect being tested. The number of columns of this
%                        matrix has to be equal to the number of conditions. 
%                        The default is a matrix that specifies the
%                        main effect of the independent variable. This matrix
%                        has size [(ncond-1),ncond]. 
%   cfg.computestat    = 'yes' or 'no', calculate the statistic (default='yes')
%   cfg.computecritval = 'yes' or 'no', calculate the critical values of the test statistics (default='no')
%   cfg.computeprob    = 'yes' or 'no', calculate the p-values (default='no')
%
% The following options are relevant if cfg.computecritval='yes' and/or
% cfg.computeprob='yes'.
%   cfg.alpha = critical alpha-level of the statistical test (default=0.05)
%   cfg.tail  = -1, 0, or 1, left, two-sided, or right (default=1)
%               cfg.tail in combination with cfg.computecritval='yes'
%               determines whether the critical value is computed at
%               quantile cfg.alpha (with cfg.tail=-1), at quantiles
%               cfg.alpha/2 and (1-cfg.alpha/2) (with cfg.tail=0), or at
%               quantile (1-cfg.alpha) (with cfg.tail=1).
%
% Design specification
%   cfg.ivar  = row number of the design that contains the labels of the conditions that must be 
%               compared (default=1). The labels range from 1 to the number of conditions.
%   cfg.uvar  = row number of design that contains the labels of the units-of-observation (subjects or trials)
%               (default=2). The labels are assumed to be integers ranging from 1 to 
%               the number of units-of-observation.

% Copyright (C) 2006, Eric Maris
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$


% set defaults
if ~isfield(cfg, 'computestat'),       cfg.computestat='yes';     end;
if ~isfield(cfg, 'computecritval'),    cfg.computecritval='no';   end;
if ~isfield(cfg, 'computeprob'),       cfg.computeprob='no';      end;
if ~isfield(cfg, 'alpha'),             cfg.alpha=0.05;            end;
if ~isfield(cfg, 'tail'),              cfg.tail=1;                end;

nconds=length(unique(design(cfg.ivar,:)));
if ~isfield(cfg,'contrastcoefs')
    % specify the default contrast coefficient matrix.
    ncontrasts = nconds-1;
    cfg.contrastcoefs = zeros(ncontrasts,nconds);
    cfg.contrastcoefs(:,1) = 1;
    for contrastindx=1:ncontrasts
        cfg.contrastcoefs(contrastindx,contrastindx+1)=-1;
    end;
else
    ncontrasts = size(cfg.contrastcoefs,1);
end;

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') & strcmp(cfg.computestat,'no')
    error('P-values can only be calculated if the test statistics are calculated.');
end;
if ~isfield(cfg,'uvar') || isempty(cfg.uvar)
    error('uvar must be specified for dependent samples statistics');
end

% perform some checks on the design
nuospercond=zeros(nconds,1);
for condindx=1:nconds
    nuospercond(condindx)=length(find(design(cfg.ivar,:)==condindx));
end;
if sum(nuospercond)<size(design,2) | nuospercond~=(nuospercond(1)*ones(nconds,1))
  error('Invalid specification of the design array.');
end;
nunits = max(design(cfg.uvar,:));
dfdenom = nunits - ncontrasts;
if dfdenom<1
    error('The data must contain more units-of-observation (usually subjects) than the number of contrasts.')
end;
nrepl=nunits*nconds;
if (nrepl~=sum(nuospercond)) | (nrepl~=size(dat,2))
  error('Invalid specification of the design array.');
end;
nsmpls = size(dat,1);

if strcmp(cfg.computestat,'yes')
    
    
   
   for c = 1:size(dat,1)
   
       if length(unique(design(1,:))) == 4
           tbl = rm_anova2(dat(c,:),design(2,:),double(design(1,:)<3),rem(design(1,:),2),{'Inherent' 'Contextual'});
           FValues(c,:) = [tbl{2:4,5}];
           s.prob(c,:) = [tbl{2:4,6}];
       elseif length(unique(design(1,:))) == 8
           tbl = rm_anova3(dat(c,:),design(2,:),double(ismember(design(1,:),[1 2 5 6])),double(ismember(design(1,:),[1 2 3 4])),double(ismember(design(1,:),[1 4 5 8])),{'Inherent' 'Contextual','Magnitude'});
           FValues(c,:) = [tbl{2:7,5}];
           s.prob(c,:) = [tbl{2:7,6}];
       else
           error('Cannot find fitting anova for design')
       end
       

    
   end
   
    s.stat = FValues;
    
end;

if length(unique(design(1,:))) == 4
    tbl = rm_anova2(dat(1,:),design(2,:),double(design(1,:)<3),rem(design(1,:),2),{'Inherent' 'Contextual'})  ;
    s.dfnum = [tbl{2:4,3}];
    s.dfdenom = tbl{5,3};
    
    if strcmp(cfg.computecritval,'yes')
        s.critval = arrayfun(@(c) finv(.95,c,tbl{5,3}),[tbl{2,3}]);
    end
    
elseif length(unique(design(1,:))) == 8
    tbl = rm_anova3(dat(1,:),design(2,:),double(ismember(design(1,:),[1 2 5 6])),double(ismember(design(1,:),[1 2 3 4])),double(ismember(design(1,:),[1 4 5 8])),{'Inherent' 'Contextual','Magnitude'});
    s.dfnum = [tbl{2:7,3}];
    s.dfdenom = [tbl{8,3}];
    
    if strcmp(cfg.computecritval,'yes')
        s.critval = arrayfun(@(c) finv(.95,c,tbl{8,3}),[tbl{2,3}]);
    end
end

if strcmp(cfg.computeprob,'yes')
  % also compute the p-values

  if cfg.tail==-1
      error('For a dependent samples F-statistic, it does not make sense to calculate a left tail p-value.');
  end;
  if cfg.tail==0
      error('For a dependent samples F-statistic, it does not make sense to calculate a two-sided p-value.');
  end;
  if cfg.tail==1
    s.prob = cell2mat(arrayfun(@(c) 1-fcdf(s.stat(:,c),s.dfnum(c),s.dfdenom),1:length(s.dfnum),'un',0));
  end;
end
