function PCAnalyzor(Deci)

    if Deci.PCom
        global PCom
        
        if isempty(PCom)
            PCom = parfeval(@numel,0,1);
        end
        
        for subject_list = 1:length(Deci.SubjectList)
           PCom(end+1)= parfeval(gcp,@Analyzor,0,Deci,subject_list);
        end
        
        delete(timerfindall)
        dc_PComTimer(Deci);
    else
      TimerAn = clock;
       for subject_list = 1:length(Deci.SubjectList)
            TimerSub = clock;
            Analyzor(Deci,subject_list);
            disp(['Analyzed for ' Deci.SubjectList{subject_list} ' Time: ' num2str(etime(clock,TimerSub))])
       end
       disp(['Full Analysis Time:' num2str(etime(clock,TimerAn))])
    end
end


