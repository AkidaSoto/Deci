function AnalyzePart1(Deci,subject_list,data,Lock)    

    condinfo = data.condinfo;
    cfg.offset = data.condinfo{1}(:,Lock);
    cfg.toilim = Deci.Analysis.Freq.Toilim;
    
    data = ft_datashift2(cfg,data);
    
    if Deci.Analysis.Laplace
        [elec.label, elec.elecpos] = elec_1020select(data.label);
        ecfg.elec = elec;
        data = ft_scalpcurrentdensity(ecfg, data);
    end
    
    data.condinfo = condinfo;
    
    mkdir([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list} filesep num2str(Lock)]);
    save([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list} filesep num2str(Lock) ],'data');
    
    if ~isfield(Deci.Analysis.Freq,'Toi')
        Deci.Analysis.Toi = [-inf inf];
        warning('Parameter for Toi not found, presuming [-inf inf]')
    end
    
    fcfg = Deci.Analysis.Freq;
    fcfg.output='fourier';
    fcfg.pad = 'maxperlen';
    fcfg.scc = 0;
    fcfg.keeptapers = 'no';
    fcfg.keeptrials = 'yes';
    fcfg.toi = Deci.Analysis.Freq.Toi(1):round(diff([data.time{1}(1) data.time{1}(2)]),5):Deci.Analysis.Freq.Toi(2);
    
    Fourier = rmfield(ft_freqanalysis(fcfg, data),'cfg');
    Fourier.condinfo = condinfo;
end
    