   
function AnalyzePart2(Deci,subject_list,freq,Lock,Chan)

        
        mkdir([Deci.Folder.Analysis filesep 'Fourier_Full' filesep  Deci.SubjectList{subject_list} filesep num2str(Lock)]);
        save([Deci.Folder.Analysis filesep 'Fourier_Full' filesep Deci.SubjectList{subject_list} filesep num2str(Lock) filesep Chan],'freq','-v7.3');
        
        
end