function PCAnalyzor3(Deci,subject_list)

if ~isfield(Deci.Analysis,'Channels')
    Deci.Analysis.Channels = 'all';
    warning('Parameter for Channels not found, presuming all')
end


if ~isfield(Deci.Analysis,'Laplace')
    Deci.Analysis.Laplace = 0;
    warning('Parameter for Laplace not found, presuming not wanted')
end


if  isempty(Deci.Analysis.Freq) &&  ~Deci.Analysis.ERP
    error('No analysis step was called for.')
end

if ~isfield(Deci.Analysis.Freq,'Toi')
    Deci.Analysis.Toi = [-inf inf];
    warning('Parameter for Toi not found, presuming [-inf inf]')
end

Deci.Analysis = Exist(Deci.Analysis,'DownSample',[]);

data = [];
load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data');
condinfo = data.condinfo;

if Deci.Analysis.Laplace
    
    if ~exist([Deci.Folder.Raw  filesep Deci.SubjectList{subject_list} '.bvct'])
        error([Deci.SubjectList{subject_list} ' does not have bvct file']);
    end
    
    [elec.label, elec.elecpos] = CapTrakMake([Deci.Folder.Raw  filesep Deci.SubjectList{subject_list} '.bvct']);
    ecfg.elec = elec;
    data = ft_scalpcurrentdensity(ecfg, data);
end

if ~isempty(Deci.Analysis.DownSample)
    data = ft_resampledata(struct('resamplefs',Deci.Analysis.DownSample,'detrend','no'),data);
end


dataplaceholder =data;

for Lock = 1:length(Deci.Analysis.Locks)
    
    cfg.offset = condinfo{1}(:,Deci.Analysis.Locks(Lock));
    cfg.toilim = Deci.Analysis.Freq.Toilim;
    
    data = ft_datashift2(cfg,dataplaceholder);
    
    data.condinfo = condinfo;
    
    mkdir([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list} filesep num2str(Deci.Analysis.Locks(Lock))]);
    save([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list} filesep num2str(Deci.Analysis.Locks(Lock)) ],'data');
    
    fcfg = Deci.Analysis.Freq;
    fcfg.output='fourier';
    fcfg.pad = 'maxperlen';
    fcfg.scc = 0;
    fcfg.keeptapers = 'no';
    fcfg.keeptrials = 'yes';
    fcfg.toi = Deci.Analysis.Freq.Toi(1):round(diff([data.time{1}(1) data.time{1}(2)]),5):Deci.Analysis.Freq.Toi(2);
    
    Fourier = rmfield(ft_freqanalysis(fcfg, data),'cfg');
    Fourier.condinfo = condinfo;
    Chan = Fourier.label;
    
    if Deci.Analysis.Var
        a = figure;
        for Channel = 1:length(Fourier.label)
            
            plot(reshape(squeeze(10*log10(mean(mean(abs(Fourier.fourierspctrm(:,Channel,:,:)).^2,3),2)))',[1 size(Fourier.fourierspctrm,1)*size(Fourier.fourierspctrm,4)]));
            hold on
        end
        a.Visible = 'on';
        title([Deci.SubjectList{subject_list} ' Chan']);
        
        waitfor(a)
        
        b = figure;
        for Channel = 1:length(Fourier.freq)
            
            plot(reshape(squeeze(10*log10(mean(mean(abs(Fourier.fourierspctrm(:,:,Channel,:)).^2,3),2)))',[1 size(Fourier.fourierspctrm,1)*size(Fourier.fourierspctrm,4)]));
            hold on
        end
        b.Visible = 'on';
        title([Deci.SubjectList{subject_list} ' Freq'])
        waitfor(b)
        
    end
    
    
    if Deci.Analysis.Clean
        fakedata = [];
        cfg =[];
        fakedata.trial = squeeze(mat2cell(10*log10(squeeze(permute(squeeze(mean(abs(Fourier.fourierspctrm).^2,3)),[2 3 1]))),size(Fourier.fourierspctrm,2),size(Fourier.fourierspctrm,4),ones(1,size(Fourier.fourierspctrm,1))));
        fakedata.time = repmat({Fourier.time},[size(Fourier.fourierspctrm,1) 1 ]);
        fakedata.label = Fourier.label;
        cfg.demean        = 'yes';
        cfg.baselinewindow = [0 0];
        [fakedata] = ft_preprocessing(cfg, fakedata);
        
        cfg =[];
        cfg.method = 'summary';
        cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
        %freqdata = ft_rejectvisual(cfg,fakedata);
        
        cfg =[];
        cfg.artfctdef.zvalue.channel = Fourier.label;
        cfg.artfctdef.zvalue.cutoff = 18.5;
        %cfg.artfctdef.zvalue.interactive = 'yes';
        %cfg.artfctdef.zvalue.cumulative = 'no';
        [cfg, ecfg.artfctdef.muscle.artifact] = ft_artifact_abszvalue(cfg,fakedata);
        
        cfg =[];
        fakedata.trial = squeeze(mat2cell(10*log10(squeeze(permute(squeeze(mean(abs(Fourier.fourierspctrm).^2,2)),[2 3 1]))),size(Fourier.fourierspctrm,3),size(Fourier.fourierspctrm,4),ones(1,size(Fourier.fourierspctrm,1))));
        fakedata.time = repmat({Fourier.time},[size(Fourier.fourierspctrm,1) 1 ]);
        fakedata.label = strsplit(num2str(Fourier.freq),' ');
        cfg.demean        = 'yes';
        cfg.baselinewindow = [0 0];
        [fakedata] = ft_preprocessing(cfg, fakedata);
        %             cfg =[];
        %             cfg.method = 'summary';
        %             cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
        %             chandata = ft_rejectvisual(cfg,fakedata);
        
        cfg =[];
        cfg.artfctdef.zvalue.channel = strsplit(num2str(Fourier.freq),' ');
        cfg.artfctdef.zvalue.cutoff = 18.5;
        %cfg.artfctdef.zvalue.interactive = 'yes';
        %cfg.artfctdef.zvalue.cumulative = 'no';
        [cfg, ecfg.artfctdef.eog.artifact] = ft_artifact_abszvalue(cfg,fakedata);
        
        
        fakedata = ft_rejectartifact(ecfg, fakedata);
        
    else
        fakedata.saminfo = logical(ones(1,size(freq.condinfo{2},1)));
    end
    
    scfg.trials = logical(fakedata.saminfo);
    
    Fourier = ft_selectdata(scfg,Fourier);
    
    %Fourier.condinfo = Fourier.condinfo{2}(logical(fakedata.saminfo),:);
    
    tic;
    for i = 1:length(Chan)
        
        dcfg = [];
        dcfg.channel = Chan(i);
        freq = ft_selectdata(dcfg,Fourier);
        freq.condinfo = Fourier.condinfo;
        
        %         mkdir([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}  filesep num2str(Deci.Analysis.Locks(Lock))]);
        %         mkdir([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep  Deci.SubjectList{subject_list}  filesep num2str(Deci.Analysis.Locks(Lock))]);
        %
        %
        %         freq = freqplaceholder;
        %         freq.dimord = 'chan_freq_time';
        %         freq.powspctrm      = permute(abs(mean(freq.fourierspctrm./abs(freq.fourierspctrm),1)),[2 3 4 1]);         % divide by amplitude
        %         freq  = rmfield(freq,'fourierspctrm');
        %         save([Deci.Folder.Analysis filesep 'Freq_ITPC' filesep Deci.SubjectList{subject_list}  filesep num2str(Deci.Analysis.Locks(Lock)) filesep Chan{i}],'freq','-v7.3');
        %
        %         freq = freqplaceholder;
        %         freq.powspctrm = permute(mean(abs(freq.fourierspctrm).^2 ,1),[2 3 4 1]);
        %         freq.dimord = 'chan_freq_time';
        %         freq  = rmfield(freq,'fourierspctrm');
        %         save([Deci.Folder.Analysis filesep 'Freq_TotalPower' filesep Deci.SubjectList{subject_list}  filesep num2str(Deci.Analysis.Locks(Lock)) filesep Chan{i}],'freq','-v7.3');
        %
        %         if Deci.Analysis.Freq.kptrls
        %freq = freqplaceholder;
        %freq.dimord = 'rpt_chan_freq_time';
        %label.dimord = 'rpt_chan_freq_time';
        mkdir([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep  Deci.SubjectList{subject_list}  filesep num2str(Deci.Analysis.Locks(Lock))]);
        save([Deci.Folder.Analysis filesep 'Four_TotalPower' filesep Deci.SubjectList{subject_list}  filesep num2str(Deci.Analysis.Locks(Lock)) filesep Chan{i}],'freq','-v7.3');
        %         end
        
        
    end
    toc;
    clear Fourier
end



end


