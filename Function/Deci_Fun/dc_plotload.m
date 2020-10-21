function [Subjects,info] =  dc_plotload(Deci,info)
%JC 6/3/20: replaced all instances of 'var' with 'variable' since 'var' is
%reserved

for  subject_list = 1:length(Deci.SubjectList)
    
    display(['Loading Plottor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
    
    for Conditions = 1:length(Deci.Plot.CondTitle)
        for Channel = 1:length(info.Chois)
            
            eval([info.variable '= [];']);
            
            if exist([Deci.Folder.Analysis filesep info.extension filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.CondTitle{Conditions} filesep info.Chois{Channel} '.mat']) == 2
                load([Deci.Folder.Analysis filesep info.extension filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.CondTitle{Conditions} filesep info.Chois{Channel} '.mat'],info.variable);
            
            else
                display(['could not find' [Deci.Folder.Analysis filesep info.extension filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.CondTitle{Conditions} filesep info.Chois{Channel} '.mat']])
                continue
            end
            
            evalc(['variable =' info.variable ';']);
            if isfield(variable,'time')
                info.toi = round(variable.time,4) >= info.Tois(1) & round(variable.time,4) <= info.Tois(2);
            end
            
            if isfield(variable,'freq')
                foi = variable.freq >= round(info.Fois(1),4) & variable.freq <= round(info.Fois(2),4);
                
                try
                    Foi(subject_list,:) = variable.freq(foi);
                catch
                    if subject_list > 1
                        error('mismatch in frequencies, previous subjects have \n%s \n while current subject %s has \n %s \n try reanalyzing?', regexprep(num2str(round(Foi(subject_list-1,:),2)),'\s+',' '),Deci.SubjectList{subject_list},regexprep(num2str(round(variable.freq(foi),2)),'\s+',' '))
                    end
                end
                
                if ~strcmp(info.variable,'MI_full')
                variable.freq = variable.freq(foi);
                variable.(info.parameter) = variable.(info.parameter)(:,foi,:);
                end
                
                Chans{Channel} = variable;
                
            end
        end
        
        if isempty(eval(info.variable))
            continue;
        end
            
        if isfield(variable,'trllength')
            info.trllen(subject_list,Conditions) = variable.trllength;
        else
            info.trllen(subject_list,Conditions) = nan;
        end
        
        if isfield(variable,'lockers')
            LockNum = Deci.Analysis.Locks(ismember(Deci.Analysis.LocksTitle,Deci.Plot.Lock));
            info.lockers(subject_list,Conditions,:) = variable.lockers - variable.lockers(LockNum);
        else
            info.lockers(subject_list,Conditions,:) = nan;
        end
        
        acfg.parameter = info.parameter;
        acfg.appenddim = 'chan';
        
        %this was breaking - JC 4/27/2020
        if isfield(variable,'freq')
            Subjects{subject_list,Conditions} = rmfield(ft_appendfreq(acfg,Chans{:}),'cfg');
        else
            Subjects{subject_list,Conditions} = rmfield(ft_appendtimelock(acfg,Chans{:}),'cfg');
        end
        %Subjects{subject_list,Conditions}.dimord = 'chan_freq_time';
    end
    clear Chans;
end

end