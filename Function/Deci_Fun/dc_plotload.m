function [Subjects,info] =  dc_plotload(Deci,info)

for  subject_list = 1:length(Deci.SubjectList)
    
    display(['Loading Plottor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
    
    for Conditions = 1:length(Deci.Plot.CondTitle)
        for Channel = 1:length(info.Chois)
            
            load([Deci.Folder.Analysis filesep info.extension filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.CondTitle{Conditions} filesep info.Chois{Channel} '.mat'],info.variable);
            
            evalc(['var =' info.variable ';']);
            
            info.toi = round(var.time,4) >= info.Tois(1) & round(var.time,4) <= info.Tois(2);
            
            if isfield(var,'freq')
                foi = var.freq >= round(info.Fois(1),4) & var.freq <= round(info.Fois(2),4);
                
                try
                Foi(subject_list,:) = var.freq(foi);
                catch
                  error('mismatch in frequencies, previous subjects have \n%s \n while current subject %s has \n %s \n try reanalyzing?', regexprep(num2str(round(Foi(subject_list-1,:),2)),'\s+',' '),Deci.SubjectList{subject_list},regexprep(num2str(round(var.freq(foi),2)),'\s+',' '))
                end
                
                var.freq = var.freq(foi);
                var.(info.parameter) = var.(info.parameter)(:,foi,:);
                
                Chans{Channel} = var;
                
            end
        end
        
        if isfield(var,'trllength')
            info.trllen(subject_list,Conditions) = var.trllength;
        else
            info.trllen(subject_list,Conditions) = nan;
        end
        
        if isfield(var,'lockers')
            info.lockers(subject_list,Conditions,:) = var.lockers;
        else
            info.lockers(subject_list,Conditions,:) = nan;
        end
        
        acfg.parameter = info.parameter;
        acfg.appenddim = 'chan';
        
        if isfield(var,'freq')
        Subjects{subject_list,Conditions} = rmfield(ft_appendfreq(acfg,Chans{:}),'cfg');
        else
        Subjects{subject_list,Conditions} = rmfield(ft_appendtimelock(acfg,Chans{:}),'cfg');    
        end
        %Subjects{subject_list,Conditions}.dimord = 'chan_freq_time';
    end
    clear Chans;
end

end