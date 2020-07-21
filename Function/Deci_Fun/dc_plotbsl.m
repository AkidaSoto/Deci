function Subjects = dc_plotbsl(Deci,Subjects,info)


for Conditions = 1:size(Subjects,2)
    for subject_list = 1:size(Subjects,1)
        
        if ~strcmpi(Deci.Plot.BslRef,Deci.Plot.Lock) || ~isempty(Deci.Plot.LockCond)
            
            if ~isempty(Deci.Plot.LockCond)
                BslCond =    Deci.Plot.CondTitle{Deci.Plot.LockCond(Conditions)};
            else
                BslCond =    Deci.Plot.CondTitle{Conditions};
            end
            
            for Channel = 1:length(info.Chois)
                
                load([Deci.Folder.Analysis filesep info.extension filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.BslRef filesep BslCond filesep info.Chois{Channel} '.mat'],info.variable);
                evalc(['vari =' info.variable]);
                
                if isfield(vari,'freq')
                    foi = vari.freq >= round(info.Fois(1),4) & vari.freq <= round(info.Fois(2),4);
                    
                    vari.freq = vari.freq(foi);
                    vari.(info.parameter) = vari.(info.parameter)(:,foi,:);
                    
                    Chans{Channel} = vari;
                    
                end
            end
            
            acfg.parameter = info.parameter;
            acfg.appenddim = 'chan';
            
            if isfield(vari,'freq')
                Bsl{subject_list,Conditions} = rmfield(ft_appendfreq(acfg,Chans{:}),'cfg');
            else
                Bsl{subject_list,Conditions} = rmfield(ft_appendtimelock(acfg,Chans{:}),'cfg');
            end
            
        else
            Bsl{subject_list,Conditions} =Subjects{subject_list,Conditions};
        end
        
        toi = Bsl{subject_list,Conditions}.time >= round(Deci.Plot.Bsl(1),4) & Bsl{subject_list,Conditions}.time <= round(Deci.Plot.Bsl(2),4);
        
        if info.isfreq
            Bsl{subject_list,Conditions}.(info.parameter) = nanmean(Bsl{subject_list,Conditions}.(info.parameter)(:,:,toi),3);
            Bsl{subject_list,Conditions}.(info.parameter) = repmat(Bsl{subject_list,Conditions}.(info.parameter),[1 1 size(Subjects{subject_list,Conditions}.(info.parameter),3)]);
        else
            Bsl{subject_list,Conditions}.(info.parameter) = nanmean(Bsl{subject_list,Conditions}.(info.parameter)(:,toi),2);
            Bsl{subject_list,Conditions}.(info.parameter) = repmat(Bsl{subject_list,Conditions}.(info.parameter),[1 size(Subjects{subject_list,Conditions}.(info.parameter),3)]);
        end
        

        switch Deci.Plot.BslType
            case 'none'
            case 'absolute'
                Subjects{subject_list,Conditions}.(info.parameter) =  Subjects{subject_list,Conditions}.(info.parameter) - Bsl{subject_list,Conditions}.(info.parameter);
            case 'relative'
                Subjects{subject_list,Conditions}.(info.parameter)=  Subjects{subject_list,Conditions}.(info.parameter) ./ Bsl{subject_list,Conditions}.(info.parameter);
            case 'relchange'
                Subjects{subject_list,Conditions}.(info.parameter) = ( Subjects{subject_list,Conditions}.(info.parameter) - Bsl{subject_list,Conditions}.(info.parameter)) ./ Bsl{subject_list,Conditions}.(info.parameter);
            case 'db'
                Subjects{subject_list,Conditions}.(info.parameter) = 10*log10( Subjects{subject_list,Conditions}.(info.parameter) ./ Bsl{subject_list,Conditions}.(info.parameter));
        end
        
        Subjects{subject_list,Conditions}.time = Subjects{subject_list,Conditions}.time(info.toi);
        
        if info.isfreq
            Subjects{subject_list,Conditions}.(info.parameter) = Subjects{subject_list,Conditions}.(info.parameter)(:,:,info.toi);
        else
            Subjects{subject_list,Conditions}.(info.parameter) = Subjects{subject_list,Conditions}.(info.parameter)(:,info.toi); 
        end

    end
end

end