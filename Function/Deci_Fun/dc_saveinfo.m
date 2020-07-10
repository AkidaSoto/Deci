function dc_saveinfo(Folder)

Files = CleanDir(Folder);

for subj = 1:length(Files)
   
    load([Folder filesep Files{subj}]);
    
    info = rmfield(data,'trial');
    
    save([Folder filesep Files{subj}],'data','info');
    
end


end