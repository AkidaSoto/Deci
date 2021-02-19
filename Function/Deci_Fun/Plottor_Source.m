function Plottor_Source(Deci,Params)

%% Load
cfg        = [];

display(['Plotting Source']);




for subject_list = 1:length(Deci.SubjectList)
    
    display(['Loading Plottor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
    
    for Conditions = 1:length(Deci.Plot.CondTitle)
        
        load([Deci.Folder.Analysis filesep 'Extra' filesep 'Source' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.CondTitle{Conditions} filesep Params.type '_' Params.foi  '.mat'],'source'); 
        source.avg.pow = mean(source.avg.pow,2);
        Subjects{subject_list,Conditions} = source;
        
    end
    
end

%% Baseline Correction


if length(Subjects{subject_list,Conditions}.time) > 1
    display(' ');
    display(['Using Lock: ' Deci.Plot.Lock]);
    display(['Using Ref: ' Deci.Plot.BslRef ' at times ' strrep(regexprep(num2str(Deci.Plot.Freq.Bsl),' +',' '),' ','-')]);
    
    if Params.Bsl
        for Conditions = 1:size(Subjects,2)
            for subject_list = 1:size(Subjects,1)
                
                if ~strcmpi(Deci.Plot.BslRef,Deci.Plot.Lock)
                    
                    load([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.BslRef filesep Deci.Plot.CondTitle{Conditions} filesep Deci.Plot.Extra.Source.List{SourceList} '.mat'],'conn');
                    
                    ccfg.latency = Deci.Plot.Freq.Bsl;
                    ccfg.avgovertime = 'yes';
                    
                    Bsl{subject_list,Conditions} = conn;
                    
                    toi1 = round(Bsl{subject_list,Conditions}.time,4) >= Deci.Plot.Freq.Bsl(1) & round(Bsl{subject_list,Conditions}.time,4) <= Deci.Plot.Freq.Bsl(2);
                    
                    Bsl{subject_list,Conditions}.param =  Bsl{subject_list,Conditions}.param(:,:,toi1);
                    Bsl{subject_list,Conditions}.time = Bsl{subject_list,Conditions}.time(toi1);
                    Bsl{subject_list,Conditions} = ft_selectdata(ccfg, Bsl{subject_list,Conditions});
                    Bsl{subject_list,Conditions}.param = repmat(Bsl{subject_list,Conditions}.param,[1 1 size(Subjects{subject_list,Conditions}.param ,3)]);
                    
                else
                    ccfg.latency = Deci.Plot.Freq.Bsl;
                    ccfg.avgovertime = 'yes';
                    
                    toi1 = round(Subjects{subject_list,Conditions}.time,4) >= Deci.Plot.Freq.Bsl(1) & round(Subjects{subject_list,Conditions}.time,4) <= Deci.Plot.Freq.Bsl(2);
                    Bsl{subject_list,Conditions} =Subjects{subject_list,Conditions};
                    Bsl{subject_list,Conditions}.param =  Bsl{subject_list,Conditions}.param(:,:,toi1);
                    Bsl{subject_list,Conditions}.time = Bsl{subject_list,Conditions}.time(toi1);
                    Bsl{subject_list,Conditions} = ft_selectdata(ccfg, Bsl{subject_list,Conditions});
                    Bsl{subject_list,Conditions}.param = repmat(Bsl{subject_list,Conditions}.param,[1 1 size(Subjects{subject_list,Conditions}.param ,3)]);
                end
                
                switch Deci.Plot.Freq.BslType
                    case 'none'
                    case 'absolute'
                        Subjects{subject_list,Conditions}.param =  Subjects{subject_list,Conditions}.param - Bsl{subject_list,Conditions}.param;
                    case 'relative'
                        Subjects{subject_list,Conditions}.param=  Subjects{subject_list,Conditions}.param ./ Bsl{subject_list,Conditions}.param;
                    case 'relchange'
                        Subjects{subject_list,Conditions}.param = ( Subjects{subject_list,Conditions}.param - Bsl{subject_list,Conditions}.param) ./ Bsl{subject_list,Conditions}.param;
                    case 'db'
                        Subjects{subject_list,Conditions}.param = 10*log10( Subjects{subject_list,Conditions}.param ./ Bsl{subject_list,Conditions}.param);
                end
                
                %Subjects{subject_list,Conditions}.time = Subjects{subject_list,Conditions}.time(toi);
                %Subjects{subject_list,Conditions}.param = Subjects{subject_list,Conditions}.param(:,:,toi);
                
            end
        end
    end
    
end

%% Math
if ~isempty(Deci.Plot.Math)
    
    display(' ')
    display(['Doing ' num2str(length(Deci.Plot.Math)) ' Maths'] )
    
    for conds = 1:length(Deci.Plot.Math)
        for subj = 1:size(Subjects,1)
            scfg.parameter = 'avg.pow';
            scfg.operation = Deci.Plot.Math{conds};
            evalc('MathData{subj} = ft_math(scfg,Subjects{subj,:})');
        end
        
        Subjects(:,size(Subjects,2)+1) = MathData;
    end
end

clear MathData
%% Data Management
if size(Subjects,1) == 1
    Deci.Plot.GrandAverage = false;
end


for conds = 1:size(Subjects,2)
    
    if Deci.Plot.GrandAverage
        facfg.parameter =  'avg.pow';
        facfg.type = 'mean';
        facfg.keepindividual = 'yes';
        evalc('SourceData{1,conds} = ft_sourcegrandaverage(facfg,Subjects{:,conds});');
        %ConnData{1,conds}.param = ConnData{1,conds}.individual;
        SourceData{1,conds} = rmfield(SourceData{1,conds},'cfg');
        
        SourcetData{1,conds} = SourceData{1,conds};
        
        SourcetData{1,conds}.pow = squeeze(SourcetData{1,conds}.('pow'));
        
        SourceData{1,conds}.pow = mean(SourceData{1,conds}.pow,1);
        SourceData{1,conds}.inside = source.inside;
        SourceData{1,conds}.dim = source.dim;
        SourceData{1,conds}.pos = source.pos;
    else
        Params.Stat = false;
        SourceData(:,conds) = Subjects(:,conds);
    end
end

s1 = size(Subjects,1) > 1;
clear Subjects
clear source
%% Stat
Params = Exist(Params,'Stat',false);

if Deci.Plot.GrandAverage && s1
    
    display(' ')
    display('Calculating Statistics')
    
    neighbours       = load('easycap_neighbours','neighbours');
    Deci.Plot.Stat.neighbours = neighbours.neighbours;
    Deci.Plot.Stat.ivar = 1;
    
    Deci.Plot.Stat.parameter = 'pow';
    
    for conds = 1:length(Deci.Plot.Draw)
        design = [];
        
        Deci.Plot.Stat.uvar = 2;
        
        for subcond = 1:length(Deci.Plot.Draw{conds})
            for subj = 1:length(Deci.SubjectList)
                design(1,subj+length(Deci.SubjectList)*[subcond-1]) =  subcond;
                design(2,subj+length(Deci.SubjectList)*[subcond-1]) = subj;
            end
        end
        
        design = design';
        
        Deci.Plot.Stat.design = design;
        
        if length(Deci.Plot.Draw{conds}) > 2
            Deci.Plot.Stat.tail = 1;
            Deci.Plot.Stat.statistic = 'depsamplesFmultivariate';
            Deci.Plot.Stat.clustertail      = 1;
        else
            Deci.Plot.Stat.statistic = 'depsamplesT';
            Deci.Plot.Stat.tail = 0;
            Deci.Plot.Stat.clustertail      = 0;
        end
        [Sourcestat{conds}] = ft_sourcestatistics(Deci.Plot.Stat, SourcetData{:,Deci.Plot.Draw{conds}});
        Sourcestat{conds}.inside = SourceData{1,conds}.inside;
        Sourcestat{conds}.dim = SourceData{1,conds}.dim;
        Sourcestat{conds}.pos = SourceData{1,conds}.pos;
    end
end

clear SourcetData
%% Plot

if Deci.Plot.GrandAverage
    SubjectList = {'Group Average'};
end

mri = ft_read_mri('standard_seg.mat');
mri = ft_convert_units(mri, 'cm');

% cfg = [];
% cfg.method = 'interactive';
% mri = ft_volumerealign(cfg,mri);

cfg  = [];
mrirs = ft_volumereslice(cfg,mri);

% cfg = [];
% ft_sourceplot(cfg,mrirs);



for cond = 1:length(Deci.Plot.Draw)
    
    for subj = 1:size(SourceData,1)
       
          
        if length(Deci.Plot.Draw{cond}) == 2
            
            cfg              = [];
            cfg.parameter    = 'pow';
            cfg.interpmethod = 'nearest';
            
            subtraction = SourceData{subj,Deci.Plot.Draw{cond}(2)};
            subtraction.pow = (subtraction.pow - SourceData{subj,Deci.Plot.Draw{cond}(1)}.pow)./ subtraction.pow;
            
            if Deci.Plot.GrandAverage && s1
                
                try
                subtraction.pow(~Sourcestat{cond}.prob.mask) = nan;
                
                if ~any(Sourcestat{cond}.prob.mask)
                    continue
                end
                
                catch
                subtraction.pow(~Sourcestat{cond}.mask) = nan;
                   
                if ~any(Sourcestat{cond}.mask)
                    continue
                end
                
                end
                

            end
            source_diff_int  = ft_sourceinterpolate(cfg, subtraction, mrirs);
            
            source_diff_int.mask = isnan(source_diff_int.pow);
            
            cfg               = [];
            cfg.method        = 'slice';
            cfg.funparameter  = 'pow';
            cfg.maskparameter  = 'mask';
            cfg.opacitymap    =  'rampdown';
            if Deci.Plot.GrandAverage && s1
                
                %cfg.maskparameter = 'stat';
                %cfg.funcolorlim   = 'zeromax';
                
                %cfg.opacitylim    = 'minzero';
            end
            source_diff_int.pow = mean(source_diff_int.pow,2);
            ft_sourceplot(cfg, source_diff_int);
            title([Deci.Plot.Subtitle{cond}{2} ' - ' Deci.Plot.Subtitle{cond}{1}]);

        end
          
          
        for subcond = 1:length(Deci.Plot.Draw{cond})
            
            
            cfg              = [];
            cfg.parameter    = 'pow';
            cfg.interpmethod = 'nearest';
            
            plotdata = SourceData{subj,Deci.Plot.Draw{cond}(subcond)};
            
            if Deci.Plot.GrandAverage && s1
                
                try
                    plotdata.pow(~Sourcestat{cond}.prob.mask) = nan;
                    
                    if ~any(Sourcestat{cond}.prob.mask)
                        continue
                    end
                    
                catch
                     plotdata.pow(~Sourcestat{cond}.mask) = nan;
                     
                     if ~any(Sourcestat{cond}.mask)
                         continue
                     end
                     
                end
                
            end
            source_diff_int  = ft_sourceinterpolate(cfg, plotdata, mrirs);
            
            cfg               = [];
            cfg.method        = 'slice';
            cfg.funparameter  = 'pow';
            %cfg.maskparameter  = 'pow';
            %cfg.opacitymap    =  'rampdown';
            if Deci.Plot.GrandAverage && s1
                
                %cfg.maskparameter = 'stat';
                %cfg.funcolorlim   = 'zeromax';
                
                %cfg.opacitylim    = 'minzero';
            end
            source_diff_int.pow = mean(source_diff_int.pow,2);
            ft_sourceplot(cfg, source_diff_int);
            title(Deci.Plot.Subtitle{cond}{subcond})
            
        end
    end
end
