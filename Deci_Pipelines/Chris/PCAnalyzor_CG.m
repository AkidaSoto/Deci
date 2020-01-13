function PCAnalyzor_CG(Deci)

global PCom

if Deci.PCom
    
if isempty(PCom)
   PCom = parfeval(gcp,@numel,0,1); 
end

% if Deci.PCom
    for subject_list = 1:length(Deci.SubjectList)
        PCom(end+1)= parfeval(gcp,@Analyzor_CG,0,Deci,subject_list);
    end
else
    TimerAn = clock;
    for subject_list = 1:length(Deci.SubjectList)
        TimerSub = clock;
        Analyzor_CG(Deci,subject_list);
        disp(['Analyzed for ' Deci.SubjectList{subject_list} ' Time: ' num2str(etime(clock,TimerSub))])
    end
    disp(['Full Analysis Time:' num2str(etime(clock,TimerAn))])
end
end


