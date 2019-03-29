function Plottor_ERP(Deci)

%% File Checks

if isequal(Deci.Plot.ERP.Channel,'Reinhart-All')
    Deci.Plot.ERP.Channel = [{'AF3'  } {'AF4'  } {'AF7'  } ...
        {'AF8'  } {'AFz'  } {'C1'   } {'C2'   } {'C3'   } {'C4'   } {'C5'   } ...
        {'C6'   } {'CP1'  } {'CP2'  } {'CP3'  } {'CP4'  } {'CP5'  } {'CP6'  } ...
        {'CPz'  } {'Cz'   } {'F1'   } {'F2'   } {'F3'   } {'F4'   } {'F5'   } ...
        {'F6'   } {'F7'   } {'F8'   } {'FC1'  } {'FC2'  } {'FC3'  } {'FC4'  } ...
        {'FC5'  } {'FC6'  } {'FCz'  } {'FT7'  } {'FT8'  } {'Fz'   } {'O1'   } ...
        {'O2'   } {'Oz'   } {'P1'   } {'P2'   } {'P3'   } {'P4'   } {'P5'   } ...
        {'P6'   } {'P7'   } {'P8'   } {'PO3'  } {'PO4'  } {'PO7'  } {'PO8'  } ...
        {'POz'  } {'Pz'   } {'T7'   } {'T8'   } {'TP10' } {'TP7'  } {'TP8'  } ...
        {'TP9'  } ] ;
end

end

if ~isempty(Deci.Plot.PRP)
    Deci.Plot.GA = 0;
    
    if isempty(Deci.Plot.PRP.label)
        Deci.Plot.Freq.BslType =  0;
        warning('Cannot Find label for PRP plot, setting as 0');
    end
    
    if isempty(Deci.Plot.PRP.ScatterToi)
        Deci.Plot.Freq.ScatterToi =  [-inf inf];
        warning('Cannot Find ScatterToi for PRPplot, setting as  [-inf inf]');
    end
end

%% Load

for  Condition =  1:length(Deci.Plot.Conditions)
    for  subject_list = 1:length(Deci.SubjectList)
        
        time = [];
        load([Deci.Folder.Analysis filesep 'Volt_ERP' filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Conditions{Condition}],'time');
        
        toi = round(time.time,4) >= Deci.Plot.ERP.Toi(1) & round(time.time,4) <= Deci.Plot.ERP.Toi(2);
        time.time =  time.time(toi);
        time.avg  =time.avg(:,toi);
        
        if ~ischar(Deci.Plot.ERP.Channel)
            if ~any(ismember(Deci.Plot.ERP.Channel,time.label))
                error('Channel Parameter has additional channels not found in Data')
            end
            
            tcfg.channel = Deci.Plot.ERP.Channel;
            time = ft_selectdata(tcfg,time);
        end
        
        Subjects{subject_list} = time;
        clear time;
    end
    
    if Deci.Plot.GA
        facfg.parameter =  'avg';
        TimeData{Condition,1} = rmfield(ft_timelockgrandaverage(facfg,Subjects{:}),'cfg');
        
    else
        TimeData(Condition,:) = Subjects(:);
    end
    clear Subjects;
    
end


if Deci.Plot.Math.Type ~= 0
    
    if isequal(Deci.Plot.Math.Form,'Hemifield')
        
        for subj = 1:size(TimeData,2)
            HemiData = TimeData(:,subj);
            Hemis = length(HemiData);
            
            for Hem = 1:Hemis
                HemiData{Hem} = hemifieldflip(HemiData{Hem});
            end
            
            Hemidata = [TimeData(:,subj);HemiData];
            
            for Hem = 1:Hemis
                scfg.operation = ['x' num2str(Hem) '-x' num2str(Hem+Hemis)];
                
                scfg.parameter = 'avg';
                
                MathData{Hem,subj} = ft_math(scfg,Hemidata{:});
                
                scfg.channel =  MathData{Hem,subj}.label(cell2mat(cellfun(@(c) any(rem(c(isstrprop(c,'digit')),2)),  MathData{Hem,subj}.label,'un',0)));
                MathData{Hem,subj} = ft_selectdata(scfg, MathData{Hem,subj});
            end
            
        end
        
        TimeData = MathData;
        
    else
        
        for cond = 1:length(Deci.Plot.Math.Form)
            for subj = 1:size(TimeData,2)
                scfg.parameter = 'avg';
                scfg.operation = Deci.Plot.Math.Form{cond};
                MathData{subj} = ft_math(scfg,TimeData{:,subj});
            end
            
            TimeData(length(Deci.Plot.Conditions)+cond,:) = MathData;
        end
        
        if Deci.Plot.Math.Type == 1
            TimeData = TimeData(length(Deci.Plot.Conditions)+1:end,:);
        end
    end
    
    
    
end


end