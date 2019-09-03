function Analyzor(Deci)

    if Deci.PCom
        %#function parallel.internal.queue.evaluateRequest
        for subject_list = 1:length(Deci.SubjectList)
           PCCom(subject_list)= parfeval(gcp,@PCAnalyzor4,0,Deci,subject_list);
        end
    
        wait(PCCom);
        disp(['Mean Time is ' char(mean(minus([PCCom.FinishDateTime],[PCCom.StartDateTime])))]);
        disp(['Total Time is ' char(max(minus([PCCom.FinishDateTime],[PCCom.StartDateTime])))]);
        
    else
      TimerAn = clock;
       for subject_list = 1:length(Deci.SubjectList)
            TimerSub = clock;
            PCAnalyzor4(Deci,subject_list);
            disp(['Analyzed for ' Deci.SubjectList{subject_list} ' Time: ' num2str(etime(clock,TimerSub))])
       end
       disp(['Full Analysis Time:' num2str(etime(clock,TimerAn))])
    end

end


