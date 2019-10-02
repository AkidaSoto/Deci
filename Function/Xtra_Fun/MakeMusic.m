function MakeMusic(File,bpm,Channel,Type)



%% Load Data Segment Data
cfg = [];
cfg.dataset = File;

% event = ft_read_event(cfg.dataset);
% hdr   = ft_read_header(cfg.dataset);
% data  = ft_read_data(cfg.dataset);

%data = ft_definetrial(cfg);

cfg.detrend = 'yes';
cfg.hpfilter = 'yes';
cfg.hpfreq = 2;
data = ft_preprocessing(cfg);

tcfg.latency = [0.001 [16*3]+10];
tcfg.channel = Channel;
data = ft_selectdata(tcfg,data);
data.time{1} = data.time{1} - 5;


SampleRate = 44100;
%%
fcfg.method        = 'wavelet';                                                  % Currently only uses 'wavelet' and 'hilbert'
fcfg.foi           = [6 12 35];                           % Frequency of Interest
fcfg.width         = 7 ;                                                         % Width
fcfg.gwidth        = 3;
fcfg.pad       = 'nextpow2';

if isempty(bpm)
    fcfg.foi           = [9:14];
    fcfg.toi           = [-.5:.05:16*3];                                                  % Time Range to save
    freq = ft_freqanalysis(fcfg,data);
    fcfg.baseline = [-.5 0];
    fcfg.baselinetype = 'db';
    freq = ft_freqbaseline(fcfg,freq);
    
    [~,I] = max(mean(freq.powspctrm(ismember(freq.label,Channel(1)),:,:),3));
    
    bpm = freq.freq(I)*10;
end

bpm = bpm*2;
bpsec = bpm/60;
secpb = 1/bpsec;

fcfg.foi           = [6 12 35];
% Gwidth
fcfg.toi           = [-.5:secpb*1:16*3];                                                  % Time Range to save
freq = ft_freqanalysis(fcfg,data);
fcfg.baseline = [-.5 0];
fcfg.baselinetype = 'db';
freq = ft_freqbaseline(fcfg,freq);


fcfg.latency = [0.01 16*3];
freq = ft_selectdata(fcfg,freq);

% figure; imagesc(squeeze(freq.powspctrm));

Scale = {'D' 'E' 'F#' 'G' 'A' 'B' 'C#'};
Sequence = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1]';

noteseq = cell([length(freq.label) size(squeeze(freq.powspctrm),3)]);

for instru = 1:length(Channel)
    
    
    
    notes = nan(size(squeeze(freq.powspctrm(ismember(freq.label,Channel(instru)),:,:))));
    squeezed = squeeze(freq.powspctrm(ismember(freq.label,Channel(instru)),:,:));
    
    notes(squeezed == max(squeezed,[],1)) = 1;
    notes(squeezed == min(squeezed,[],1)) = 0;
    notes(isnan(notes)) = squeezed(isnan(notes)) > mean(squeezed,1)';
    notes(:,std(squeezed,[],1) < 3) = repmat(1,size(notes(:,std(squeezed,[],1) < 3)));
    
    
    N = cell2mat(cellfun(@(c) ismember(notes',Sequence(:,ismember(Scale,c))','rows'),Scale,'un',0));
    
    
    for n = 1:length(noteseq)
        noteseq{instru,n} = Scale{find(N(n,:))};
    end
end
%%

for instru = 1:length(Channel)
    
    tcfg.latency = [0.001 [16*3]];
    tcfg.channel = Channel(instru);
    datahold = ft_selectdata(tcfg,data);
    
    datahold.trial{1} = datahold.trial{1} - min(datahold.trial{1});
    
    datahold.time{1} = 1:prod([floor(datahold.fsample/bpsec) length(datahold.trial{1})/[datahold.fsample/bpsec]]);
    datahold.trial{1} = datahold.trial{1}(1:prod([floor(datahold.fsample/bpsec) length(datahold.trial{1})/[datahold.fsample/bpsec]]));
    
    
    datahold.time{1} = mean(reshape(datahold.time{1},[floor(datahold.fsample/bpsec) ceil(length(datahold.trial{1})/[datahold.fsample/bpsec])]),1);
    datahold.trial{1} = mean(reshape(datahold.trial{1},[floor(datahold.fsample/bpsec) ceil(length(datahold.trial{1})/[datahold.fsample/bpsec])]),1);
    datahold.fsample = bpsec;
    datahold = rmfield(datahold,'hdr');
    
    tcfg.keeptrials = 'yes';
    tcfg.vartrllength = 2;
    ERP = ft_timelockanalysis(tcfg,datahold);
    
    AverageEnergy = squeeze(mean(ERP.avg,2));
    
    [~,NotePosition{instru}] = findpeaks(ERP.avg); %+ .5*std(ERP.avg);
end
%figure;plot(ERP.avg)
%%
FinalPiano = cell(size(noteseq(1,:)));
FinalPiano(NotePosition{1}) = noteseq(1,NotePosition{1});

FinalGuitar = cell(size(noteseq(3,:)));
FinalGuitar(NotePosition{3}) = noteseq(3,NotePosition{3});


Scale = {'D' 'E' 'F#' 'G' 'A' 'B' 'C#'};
Frequencies = [296.66 329.63 369.99 392 440 493.88 554.37];

path = 'C:\Users\User\Desktop\Notes\';
PianoNotes = [{audioread([path 'piano-d.wav'])} {audioread([path 'piano-e.wav'])} ...
    {audioread([path 'piano-fsharp.wav'])} {audioread([path 'piano-g.wav'])} ...
    {audioread([path 'piano-a.wav'])} {audioread([path 'piano-b.wav'])} ...
    {audioread([path 'piano-csharp.wav'])} ];
PianoNotes = cellfun(@(c) c(1:67936,:)*.30,PianoNotes,'UniformOutput',false);

GuitarChords = [{audioread([path 'DMajor.wav'])} {audioread([path 'EMajor.wav'])} ...
    {audioread([path 'FSharp.wav'])} {audioread([path 'GMajor.wav'])} ...
    {audioread([path 'AMajor.wav'])} {audioread([path 'BMajor.wav'])} ...
    {audioread([path 'CSharp.wav'])}];
GuitarChords = cellfun(@(c) c(1:44032,:)*3,GuitarChords,'UniformOutput',false);



cadence = round([SampleRate/bpsec]);
Song = zeros([[16*3*SampleRate]+67936-[cadence] 2]);



for Hz = 1:length(FinalPiano)
    
    
    if Type(1)
        if ~isempty(FinalPiano{Hz})
            %         x1 = -.25*sin(2*pi*Frequencies(ismember(Scale,FinalNotes{Hz}))*[0:1/20000:.5]);
            %
            %         x2 = .25*sin(pi*Frequencies(ismember(Scale,FinalNotes{Hz}))*2*[0:1/20000:.5]);
            %         x3 = [sqrt(3)/2]*cos(pi*Frequencies(ismember(Scale,FinalNotes{Hz}))*4*[0:1/20000:.5]);
            %
            %         Song = [Song  x1+x2+x3];
            
            Note = Song([[Hz-1]*cadence ]+1:[[Hz-1]*cadence]+67936,1:2)+PianoNotes{ismember(Scale,FinalPiano{Hz})};
            
            if FinalPiano{Hz} == noteseq{2,Hz}
                
                chord1 = rem(find(ismember(Scale,FinalPiano{Hz})) + 2,7);
                
                if chord1 == 0
                    chord1 = 7;
                end
                
                chord2 = rem(find(ismember(Scale,FinalPiano{Hz})) + 4,7);
                
                if chord2 == 0
                    chord2 = 7;
                end
                
                Note = Note + PianoNotes{chord1} + PianoNotes{chord2};
            end
            
            
            Song([[Hz-1]*cadence ]+1:[[Hz-1]*cadence ]+67936,1:2) = Note;
        else
            %Song([[Hz-1]*48000]+1:[Hz]*48000,1:2) = 0;
        end
    end
    
    if Type(2)
        if ~isempty(FinalGuitar{Hz})
            Note = Song([[Hz-1]*cadence ]+1:[[Hz-1]*cadence]+44032,1:2)+GuitarChords{ismember(Scale,FinalGuitar{Hz})};
            Song([[Hz-1]*cadence ]+1:[[Hz-1]*cadence ]+44032,1:2) = Note;
        end
    end
end

sound(Song,SampleRate)


FinalPiano = array2table(FinalPiano');


end