function TxT(Deci,info,freq,params)

k = 0;


if ismember(info.Channels{info.ChanNum},params.Channels)
    
    
    for fh = 1:length(params.freq)
        switch params.freq{fh}
            case 'theta'
                HF(:,fh) = [4 8];
            case 'beta'
                HF(:,fh) = [12.5 30];
            case 'alpha'
                HF(:,fh) =[8 12.5];
            case 'lowgamma'
                HF(:,fh) =[30 55];
            case 'highgamma'
                HF(:,fh) = [55 80];
        end
    end
    
    for 1:length(params.freq)
        
    end
        
        save([Deci.Folder.Analysis filesep 'Txt_' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep Chan{info.ChanNum}],'freq','-v7.3');
        
    end
end

end