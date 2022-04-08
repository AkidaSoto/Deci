function PCArtifactor(Deci)
global PCom

if isempty(PCom)
    PCom = parfeval(@numel,0,1);
end

for subject_list = 1:length(Deci.SubjectList)
    PCom(end+1)= parfeval(@Artifactor,0,Deci,subject_list);
end

delete(timerfindall)
dc_PComTimer(Deci);

end