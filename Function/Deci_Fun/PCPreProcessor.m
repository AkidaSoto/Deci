function PCPreProcessor(Deci)

if Deci.PCom
    global PCom

    if isempty(PCom)
       PCom = parfeval(@numel,0,1);
    end
    
    for subject_list = 1:length(Deci.SubjectList)
        PCom(end+1)= parfeval(@PreProcessor,0,Deci,subject_list);
    end

    delete(timerfindall)
    dc_PComTimer(Deci);
else
    for subject_list = 1:length(Deci.SubjectList)
        PreProcessor(Deci,subject_list);
    end
end

end