function StatData = dc_plotstat(Deci,SegStatdata,info)

display(' ')
display('Calculating Statistics')

if isfield(SegStatdata{1},'label')
    if length(SegStatdata{1}.label) > 1
        neighbours       = load('easycapM11_neighb');
        Deci.Plot.Stat.neighbours = neighbours.neighbours;
    end
end

Deci.Plot.Stat.ivar = 1;

% if isfield(SegStatdata{1},'cfg')
%    SegStatdata = cellfun(@(c) rmfield(c,'cfg'), SegStatdata,'un',0);
% end

Deci.Plot.Stat.randomseed = 7;

for conds = 1:length(Deci.Plot.Draw)
    design = [];
    
    switch Deci.Plot.Stat.Comp
        case 'Conditions'
            if length(Deci.Plot.Draw{conds}) ~= 1 && ~Deci.Plot.GroupLevel
                Deci.Plot.Stat.uvar = 2;
                
                
                for subcond = 1:length(Deci.Plot.Draw{conds})
                    for subj = 1:size(SegStatdata{1}.(info.parameter),1)
                        design(1,subj+size(SegStatdata{1}.(info.parameter),1)*[subcond-1]) =  subcond;
                        design(2,subj+size(SegStatdata{1}.(info.parameter),1)*[subcond-1]) = subj;
                    end
                end
                
                Deci.Plot.Stat.design = design;
                
                if length(Deci.Plot.Draw{conds}) > 2
                    Deci.Plot.Stat.tail = 1;
                    
                    Deci.Plot.Stat = Exist(Deci.Plot.Stat,'twoway',[]);
                    Deci.Plot.Stat.twoway = Exist(Deci.Plot.Stat.twoway,'do',false);
                    
                    if Deci.Plot.Stat.twoway.do && length(Deci.Plot.Draw{conds}) == 4
                        Deci.Plot.Stat.statistic = 'dc_statfun_depsamplesFunivariate';
                    else
                        Deci.Plot.Stat.statistic = 'ft_statfun_depsamplesFunivariate';
                    end
                    Deci.Plot.Stat.clustertail      = 1;
                else
                    Deci.Plot.Stat.statistic = 'depsamplesT';
                    Deci.Plot.Stat.tail = 0;
                    Deci.Plot.Stat.clustertail      = 0;
                end
                
                if info.isfreq
                    [StatData{conds}] = ft_freqstatistics(Deci.Plot.Stat, SegStatdata{:,Deci.Plot.Draw{conds}});
                elseif info.isconn
                    Deci.Plot.Stat.minnbchan = 1;
                    %Deci.Plot.Stat.neighbours = [];
                    
                    if Deci.Plot.Extra.Conn.CL.do || Deci.Plot.Extra.Conn.CH.do
                         neighbours       = load('easycapM11_neighb');
                        Deci.Plot.Stat.neighbours = neighbours.neighbours;
                        
                        if Deci.Plot.Extra.Conn.CL.do 
                            Deci.Plot.Stat.channel = SegStatdata{1}.chanlow;
                        else
                            Deci.Plot.Stat.channel = SegStatdata{1}.chanhigh;
                        end
                        
                    end
                    
                    [StatData{conds}] = dc_freqstatistics(Deci.Plot.Stat, SegStatdata{:,Deci.Plot.Draw{conds}});
                    
                    %[StatData{conds}] = dc_connectivitystatistics(Deci.Plot.Stat, SegStatdata{:,Deci.Plot.Draw{conds}});
                else
                    [StatData{conds}] = ft_timelockstatistics(Deci.Plot.Stat, SegStatdata{:,Deci.Plot.Draw{conds}});
                end
                
            elseif Deci.Plot.GroupLevel 
                design = [];
                
                design(1,:) = Deci.Plot.Groups{2}+1;
                
                Deci.Plot.Stat.design = design;
                
                if length(Deci.Plot.Draw{conds}) > 2
                    
                    error('Currently does not support 2-Way Indep Anova')
                    
                    Deci.Plot.Stat.tail = 1;
                    Deci.Plot.Stat.statistic = 'indepsamplesF';
                    Deci.Plot.Stat.clustertail      = 1;
                else
                    Deci.Plot.Stat.statistic = 'indepsamplesT';
                    Deci.Plot.Stat.tail = 0;
                    Deci.Plot.Stat.clustertail      = 0;
                end
                
                if info.isfreq
                    [StatData{conds}] = ft_freqstatistics(Deci.Plot.Stat, SegStatdata{:,Deci.Plot.Draw{conds}});
                elseif info.isconn
                    Deci.Plot.Stat.minnbchan = 1;
                    [StatData{conds}] = dc_connectivitystatistics(Deci.Plot.Stat, SegStatdata{:,Deci.Plot.Draw{conds}});
                else
                    [StatData{conds}] = ft_timelockstatistics(Deci.Plot.Stat, SegStatdata{:,Deci.Plot.Draw{conds}});
                end
                
            else
                [StatData{conds}.mask] = permute(ones(size(SegStatdata{:,Deci.Plot.Draw{conds}(1)}.(info.parameter)(1,:,:,:))),[2 3 4 1]);
                [StatData{conds}.prob] = permute(ones(size(SegStatdata{:,Deci.Plot.Draw{conds}(1)}.(info.parameter)(1,:,:,:))),[2 3 4 1]);
                [StatData{conds}.stat] = permute(ones(size(SegStatdata{:,Deci.Plot.Draw{conds}(1)}.(info.parameter)(1,:,:,:))),[2 3 4 1]);
            end
            
        case 'Bsl'
            Deci.Plot.Stat.tail = 0;
            Deci.Plot.Stat.statistic = 'indepsamplesT';
            
            Deci.Plot.Stat.design = ones([size(SegStatdata{1}.(info.parameter),1) 1]);
            
            if info.isfreq
                allSegStatdata = ft_freqgrandaverage([],SegStatdata{Deci.Plot.Draw{conds}});
                allSegStatdata.dimord = 'rpt_chan_freq_time';
                [StatData{conds}] = ft_freqstatistics(Deci.Plot.Stat, allSegStatdata);
            elseif info.isconn
                
                allSegStatdata = SegStatdata{1};
                allSegStatdata.(info.parameter) = cellfun(@(c) c.(info.parameter),SegStatdata(Deci.Plot.Draw{conds}),'UniformOutput',false);
                allSegStatdata.(info.parameter) = mean(cat(4,allSegStatdata.(info.parameter){:}),4);
                [StatData{conds}] = dc_connectivitystatistics(Deci.Plot.Stat, allSegStatdata);
            else
                allSegStatdata = ft_timelockgrandaverage([],SegStatdata{Deci.Plot.Draw{conds}});
                allSegStatdata.dimord = 'rpt_chan_time';
                [StatData{conds}] = ft_timelockstatistics(Deci.Plot.Stat, allSegStatdata);
                
            end
    end
    
    if iscell(StatData{conds})
        tmp.prob = cellfun(@(c) c.prob,StatData{conds},'UniformOutput',false);
        tmp.prob = cat(4,tmp.prob{:});
        
        tmp.stat = cellfun(@(c) c.stat,StatData{conds},'UniformOutput',false);
        tmp.stat= cat(4,tmp.stat{:});
        
        tmp.mask = cellfun(@(c) c.mask,StatData{conds},'UniformOutput',false);
        tmp.mask = cat(4,tmp.mask{:});
        %tmp.critval = cell2mat(cellfun(@(c) c.critval,StatData{conds},'UniformOutput',false));
        
    else
        
        if Deci.Plot.Stat.twoway.do
            statdim = length(size(StatData{conds}.prob));
            
            tmp.prob = permute(StatData{conds}.prob,[1:statdim-1 statdim+1:statdim+[4 - statdim] statdim]);
            tmp.stat = permute(StatData{conds}.stat,[1:statdim-1 statdim+1:statdim+[4 - statdim] statdim]);
            tmp.mask = permute(StatData{conds}.mask,[1:statdim-1 statdim+1:statdim+[4 - statdim] statdim]);
            %tmp.critval = StatData{conds}.critval;
        else
            tmp.prob = StatData{conds}.prob;
            tmp.stat = StatData{conds}.stat;
            tmp.mask = StatData{conds}.mask;
        end

        
    end
    
    StatData{conds} = tmp;
end

end