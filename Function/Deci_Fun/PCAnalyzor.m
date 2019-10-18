function PCAnalyzor(Deci)
    if Deci.PCom
        for subject_list = 1:length(Deci.SubjectList)
           PCom(subject_list)= parfeval(gcp,@Analyzor,0,Deci,subject_list);
        end
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


