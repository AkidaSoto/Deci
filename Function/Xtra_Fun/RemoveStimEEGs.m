function RemoveStimEEGs(Dir,datatype,Start,End)

AllFiles = CleanDir(Dir);

AllFiles = unique(cellfun(@(d) d(1),cellfun(@(c) strsplit(c,'.'),AllFiles,'un',0)));

mkdir(strcat(Dir,'_new'));
for Each = 1:length(AllFiles)
    
    Deci.DT.Starts     = {Start};                                                                       % Cell Array of Markers for Start codes.
    Deci.DT.Ends       = {End};
    Deci.DT.Type = 'EEG_Polymerase';
    Deci.DT.Toi        = [-2 5];
    Deci.DT.Markers    = {[Start]};
    Deci.DT.Block = [];
    Deci.DT.Locks = Start;
    
    cfg = [];
   
    cfg.dataset = [Dir filesep AllFiles{Each} '.' datatype];

    cfg.DT = Deci.DT;
    cfg.trialfun = cfg.DT.Type;
    cfg.file_ending = datatype;
    Trialcfg = ft_definetrial(cfg);

    Trialcfg.detrend = 'yes';
    data = ft_preprocessing(Trialcfg);
    
    Thresh =    cellfun(@(c) mean(range(c,2)),data.trial)  < .2e4;

    figure;
    plot(cellfun(@(c) mean(range(c,2)),data.trial));
    hold on
    plot(cellfun(@(c) max(var(c,[],2)),data.trial))

    TrialThresh = find(Thresh,1,'first');
    TimeThresh = data.sampleinfo(TrialThresh+10,1);
    
    
    hdr   = ft_read_header([Dir filesep AllFiles{Each} '.' datatype]);
    event = ft_read_event([Dir filesep AllFiles{Each} '.' datatype]);
    dat   = ft_read_data([Dir filesep AllFiles{Each} '.' datatype]);
    
    event = StandardizeEventMarkers(event);

    FirstTrial = find([event.sample] > TimeThresh,1,'first');
    FirstMrk = FirstTrial + find(ismember({event(FirstTrial:end).value},{num2str(Start)}),1,'first') - 1;

    event = event(FirstMrk:end);
    dat = dat(:,event(1).sample-5000:end);
    event= arrayfun(@(c)  setfield(c,'sample',c.sample - [event(1).sample - 5000]),event);
    ft_write_data([Dir '_new' filesep  AllFiles{Each}],dat,'header',hdr,'dataformat','brainvision_eeg','event',event);
    
    %%plot
    figure;
    cfg = [];
   
    cfg.dataset = [Dir '_new' filesep AllFiles{Each} '.' datatype];

    cfg.DT = Deci.DT;
    cfg.trialfun = cfg.DT.Type;
    cfg.file_ending = datatype;
    Trialcfg = ft_definetrial(cfg);
    

    Trialcfg.detrend = 'yes';
    data = ft_preprocessing(Trialcfg);
    plot(cellfun(@(c) median(range(c,2)),data.trial));

end

end