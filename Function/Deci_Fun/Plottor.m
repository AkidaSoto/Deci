function Plottor(Deci)


%% File Checks

for subject_list = 1:length(Deci.SubjectList)
    
    if Deci.Run.Freq 
        if ~isempty(Deci.Plot.Freq) 
            if ~isdir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}])
                error(['Freq Analysis not found for '  Deci.SubjectList{subject_list}])
            end
        end
    end
    
    if Deci.Run.ERP
        if ~isempty(Deci.Plot.ERP) 
            if ~isdir([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list}])
                error(['ERP Analysis not found for '  Deci.SubjectList{subject_list}])
            end
        end
    end
    
    if Deci.Run.Behavior
        if ~isempty(Deci.Plot.Behv)
            
            switch Deci.Plot.Behv.Source
                case 'Definition'
                    if exist([Deci.Folder.Version filesep 'Definition' filesep Deci.SubjectList{subject_list} '.mat'],'file') ~= 2
                        error(['Definition (Behv) Analysis not found for '  Deci.SubjectList{subject_list}])
                    end
                case 'Freq'
                    if ~isempty(Deci.Plot.Freq) ||  ~isempty(Deci.Plot.PRP) ||  ~isempty(Deci.Plot.CFC)
                        if ~isdir([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep Deci.SubjectList{subject_list}])
                            error(['Freq Analysis not found for '  Deci.SubjectList{subject_list}])
                        end
                    end
            end
            
        end
    end
end

%% Split

if Deci.Run.Freq 
    Plottor_Freq(Deci);
end

if Deci.Run.ERP 
    Plottor_ERP(Deci);
end

if Deci.Run.Behavior 
    Plottor_Behv(Deci);
end

if Deci.Run.Extra && ~isempty(Deci.Plot.Extra)
    for funs = find(Deci.Plot.Extra.List)
        feval(Deci.Plot.Extra.Functions{funs},Deci,Deci.Plot.Extra.Params{funs}{:});
    end
end


end

