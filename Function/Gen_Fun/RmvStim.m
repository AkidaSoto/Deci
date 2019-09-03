function RmvStim(folder,chans)

        if ~isempty(dir([folder filesep '*.dat']))
            exp_files = dir([folder filesep '*.dat']);
        elseif  ~isempty(dir([folder filesep '*.bdf']))
            exp_files = dir([folder filesep '*.bdf']);
        elseif ~isempty(dir([folder filesep '*.eeg']))
            exp_files = dir([folder filesep '*.eeg']);
        elseif ~isempty(dir([folder filesep '*.cnt']))
            exp_files = dir([folder filesep '*.cnt']);
        end
        
mkdir([folder '_nonstim'])

for k = 1:length(exp_files)
    
    hdr   = ft_read_header([folder filesep exp_files(k).name]);
    event = ft_read_event([folder filesep exp_files(k).name]); %Mk<Marker number>=<Type>,<Description>,<Position in data points>,
    dat   = ft_read_data([folder filesep exp_files(k).name]);
    cfg.dataset = [folder filesep exp_files(k).name];
    
    happy = 0;
    while happy == 0
    cfg.channel = chans;
    rmvd = ft_databrowser(cfg);
    max(rmvd.artfctdef.visual.artifact) < [event.sample];
    
    happy = input('happy? (0/1)');
    end
    
    
    event = event(max(rmvd.artfctdef.visual.artifact) < [event.sample]);
    
    string = strsplit(exp_files(k).name,'.');
    ft_write_data([folder '_nonstim' filesep  string{1}],dat,'header',hdr,'dataformat','brainvision_eeg','event',event)

    
end