function Plottor_Behv(Deci)

for subject_list = 1:length(Deci.SubjectList)
    
    data = [];
    
    switch Deci.Plot.Behv.Source
        case 'Freq'
            chan1 = CleanDir([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep Deci.SubjectList{subject_list}]);
            load([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep Deci.SubjectList{subject_list} filesep chan1{1}]);
        case 'Definition'
            load([Deci.Folder.Version filesep 'Definition' filesep Deci.SubjectList{subject_list}]);
            data = cfg;
            data.trialinfo = data.trl(:,end-length(Deci.DT.Locks)+1:end);
    end
    
    
    if ~isempty(Deci.Plot.Behv.RT)
        
        if length(Deci.DT.Locks) <2 || length(Deci.Plot.Behv.RT.Locks) ~= 2
            error('DT.Locks seems to not have enough to calculate Locks(2)-Locks(1)')
        end
        
        
        if isempty(Deci.Plot.Behv.RT.Block)
            
            
            for draw = 1:length(Deci.Plot.Behv.RT)
                
                draws = Deci.Plot.Conditions(Deci.Plot.Behv.RT.Draw{draw});
                
                maxt = max(sum(ismember(data.event,[draws{:}]),2));
                trl = sum(ismember(data.event,[draws{:}]),2) == maxt;
                
                RT(draw,subject_list) =  mean([-data.trialinfo(trl,Deci.Plot.Behv.RT.Locks(2)) - -data.trialinfo(trl,Deci.Plot.Behv.RT.Locks(1))]);
                
            end
            
        end
        disp(RT)
    end
    
    
    if ~isempty(Deci.Plot.Behv.Acc)
        Exist(Deci.Plot.Behv.Acc,'Total');
        Exist(Deci.Plot.Behv.Acc,'Subtotal');
        
                
        if isempty(Deci.Plot.Behv.RT.Block)
            
            for draw = 1:length(Deci.Plot.Behv.Acc.Total)
                
                draws = Deci.Plot.Conditions(Deci.Plot.Behv.Acc.Total{draw});
                
                maxt = max(sum(ismember(data.event,[draws{:}]),2));
                trl = sum(ismember(data.event,[draws{:}]),2) == maxt;
                Total = data.event(trl,:);
                
                subdraws = Deci.Plot.Conditions(Deci.Plot.Behv.Acc.Subtotal{draw});
                
                maxt = max(sum(ismember(data.event,[subdraws{:}]),2));
                subtrl = sum(ismember(Total,[subdraws{:}]),2) == maxt;
                
                Acc(draw,subject_list) = length(find(subtrl))/length(find(trl));
                
            end
            
        else
            
        end

        disp(Acc);
    end
    
    
end

end