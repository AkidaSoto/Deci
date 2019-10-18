function PCPreProcessor(Deci)

if Deci.PCom
    
    for subject_list = 1:length(Deci.SubjectList)
        PCom(subject_list)= parfeval(@PreProcessor,0,Deci,subject_list);
    end
else
    for subject_list = 1:length(Deci.SubjectList)
        PreProcessor(Deci,subject_list);
    end
end

end