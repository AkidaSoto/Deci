function Subjects = dc_math(Deci,Subjects,info)   

display(' ')
    display(['Doing ' num2str(length(Deci.Plot.Math)) ' Maths'] )

    for conds = 1:length(Deci.Plot.Math)
        for subj = 1:size(Subjects,1)
            scfg.parameter = info.parameter;
            scfg.operation = Deci.Plot.Math{conds};
            evalc('MathData{subj} = ft_math(scfg,Subjects{subj,:})');
        end 
        
        Subjects(:,size(Subjects,2)+1) = MathData;
    end
    
end