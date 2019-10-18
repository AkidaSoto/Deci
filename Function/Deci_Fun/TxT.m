function TxT(Deci,info,freq,params)

k = 0;


%freq.powspctrm      = permute(abs(mean(freq.fourierspctrm./abs(freq.fourierspctrm),1)),[2 3 4 1]);


%This (below) colapses across trials and gives a (channel x Freq x time)
%matrix (64,40,1001)
%freq.powspctrm = permute(mean(abs(freq.fourierspctrm).^2 ,1),[2 3 4 1]);


save([Deci.Folder.Analysis filesep 'Txt_' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep Chan{info.ChanNum}],'freq','-v7.3');
              
end