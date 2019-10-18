function MakeMusic(File,bpm,Channel,Type,Image)



%% Load Data Segment Data
cfg = [];
cfg.dataset = File;

% event = ft_read_event(cfg.dataset);
% hdr   = ft_read_header(cfg.dataset);
% data  = ft_read_data(cfg.dataset);

%data = ft_definetrial(cfg);

cfg.trialdef.eventtype = 'Stim';
cfg.trialdef.eventvalue = 'S 10';
cfg.trialdef.poststim = 54;
cfg.trialdef.prestim = 0;


cfg = ft_definetrial(cfg);

cfg.detrend = 'yes';
cfg.hpfilter = 'yes';
cfg.hpfreq = 2;
data = ft_preprocessing(cfg);

ecfg.trials = logical(zeros(1,5));

ecfg.trials(Image) = 1;
ecfg.latency = [0.001 [16*3]+6];
ecfg.channel = Channel;
data = ft_selectdata(ecfg,data);
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
    
    Topo{instru} = notes;
    
end
%%
% 
% SquareData.powspctrm = double(permute(cat(3,Topo{:}),[4 3 1 2]));
% SquareData.freq = [6 12 25];
% SquareData.time = 1:size(Topo{1},2);
% SquareData.label = Channel;
% SquareData.dimord = 'rpt_chan_freq_time';
% SquareData.trialinfo = 1;
% cfg        = [];
% cfg.layout = 'C:\Users\User\Documents\GitHub\Deci\Function\Gen_Fun\easycap_rob_eye.mat';
% cfg.channel = Channel(1)';
% cfg.trials = 1;
% 
% 
% Framer = 0:1:size(Topo{1},2);
% for a = 1:length(Framer)
%     topofig = figure;
%     
%     ft_singleplotTFR(cfg,SquareData);
%     hold on
%     xticks(linspace(1,size(Topo{1},2),size(Topo{1},2)))
%     xticklabels(noteseq(1,:))
%     yticks([6 18 25]);
%     yticklabels({'Theta' 'Alpha' 'Beta'});
%     
%     timer = plot([Framer(a) Framer(a)],topofig.Children(3).YLim,'Color','r','LineWidth',4);
%     drawnow
%     F(a) = getframe; 
%     close all
% end


 fig = figure;
for instru = 1:length(Channel)
    
    tcfg.latency = [0.001 [16*3]];
    tcfg.channel = Channel(instru);
    datahold = ft_selectdata(tcfg,data);
    
    datahold.time{1} = 1:prod([floor(datahold.fsample/bpsec) length(datahold.trial{1})/[datahold.fsample/bpsec]]);
    datahold.trial{1} = datahold.trial{1}(1:prod([floor(datahold.fsample/bpsec) length(datahold.trial{1})/[datahold.fsample/bpsec]]));
    
    
    datahold.time{1} = mean(reshape(datahold.time{1},[floor(datahold.fsample/bpsec) ceil(length(datahold.trial{1})/floor(datahold.fsample/bpsec))]),1);
    datahold.trial{1} = mean(reshape(datahold.trial{1},[floor(datahold.fsample/bpsec) ceil(length(datahold.trial{1})/floor(datahold.fsample/bpsec))]),1);
    datahold.fsample = bpsec;
    datahold = rmfield(datahold,'hdr');
    datahold.trial{1} = datahold.trial{1} - min(datahold.trial{1});
        
    tcfg.keeptrials = 'yes';
    tcfg.vartrllength = 2;
    ERP = ft_timelockanalysis(tcfg,datahold);
    
    AverageEnergy = squeeze(mean(ERP.avg,2));
    
    [~,NotePosition{instru}] = findpeaks(datahold.trial{1}); %+ .5*std(ERP.avg);
    
    findpeaks(datahold.trial{1});
    hold on;
    
end

%figure;plot(ERP.avg)
%%
FinalMelody = cell(size(noteseq(1,:)));
FinalMelody(NotePosition{1}) = noteseq(1,NotePosition{1});

FinalBass = cell(size(noteseq(3,:)));
FinalBass(NotePosition{3}) = noteseq(3,NotePosition{3});


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
%GuitarChords = cellfun(@(c) c(1:44032,:)*3,GuitarChords,'UniformOutput',false);

path = 'C:\Users\User\Desktop\Notes\Vibrophone\';

VibroPhone = [{audioread([path 'D.wav'])} {audioread([path 'E.wav'])} ...
    {audioread([path 'FSharp.wav'])} {audioread([path 'G.wav'])} ...
    {audioread([path 'A.wav'])} {audioread([path 'B.wav'])} ...
    {audioread([path 'CSharp.wav'])}];
VibroPhone = cellfun(@(c) c*1,VibroPhone,'UniformOutput',false);

path = 'C:\Users\User\Desktop\Notes\TenorSax\';
TenorSax = [{audioread([path 'D.wav'])} {audioread([path 'E.wav'])} ...
    {audioread([path 'FSharp.wav'])} {audioread([path 'G.wav'])} ...
    {audioread([path 'A.wav'])} {audioread([path 'B.wav'])} ...
    {audioread([path 'CSharp.wav'])}];
TenorSax = cellfun(@(c) c*.75,TenorSax,'UniformOutput',false);


path = 'C:\Users\User\Desktop\Notes\Atmosphere\';
Atmosphere = [{audioread([path 'D.wav'])} {audioread([path 'E.wav'])} ...
    {audioread([path 'FSharp.wav'])} {audioread([path 'G.wav'])} ...
    {audioread([path 'A.wav'])} {audioread([path 'B.wav'])} ...
    {audioread([path 'CSharp.wav'])}];
Atmosphere = cellfun(@(c) c*2,Atmosphere,'UniformOutput',false);

path = 'C:\Users\User\Desktop\Notes\Goblin\';
Goblin = [{audioread([path 'D.wav'])} {audioread([path 'E.wav'])} ...
    {audioread([path 'FSharp.wav'])} {audioread([path 'G.wav'])} ...
    {audioread([path 'A.wav'])} {audioread([path 'B.wav'])} ...
    {audioread([path 'CSharp.wav'])}];
Goblin = cellfun(@(c) c*1,Goblin,'UniformOutput',false);


path = 'C:\Users\User\Desktop\Notes\FenderJazz\';
FenderJazz = [{audioread([path 'D.wav'])} {audioread([path 'E.wav'])} ...
    {audioread([path 'FSharp.wav'])} {audioread([path 'G.wav'])} ...
    {audioread([path 'A.wav'])} {audioread([path 'B.wav'])} ...
    {audioread([path 'CSharp.wav'])}];
FenderJazz = cellfun(@(c) c*2,FenderJazz,'UniformOutput',false);


path = 'C:\Users\User\Desktop\Notes\Cello\';
Cello = [{audioread([path 'D.wav'])} {audioread([path 'E.wav'])} ...
    {audioread([path 'FSharp.wav'])} {audioread([path 'G.wav'])} ...
    {audioread([path 'A.wav'])} {audioread([path 'B.wav'])} ...
    {audioread([path 'CSharp.wav'])}];
Cello = cellfun(@(c) c*.5,Cello,'UniformOutput',false);

path = 'C:\Users\User\Desktop\Notes\Cello2\';
Cello2 = [{audioread([path 'D.wav'])} {audioread([path 'E.wav'])} ...
    {audioread([path 'FSharp.wav'])} {audioread([path 'G.wav'])} ...
    {audioread([path 'A.wav'])} {audioread([path 'B.wav'])} ...
    {audioread([path 'CSharp.wav'])}];
Cello2 = cellfun(@(c) c*.8,Cello2,'UniformOutput',false);

MelodyInstru = PianoNotes;
BassInstru = TenorSax;

cadence = round([SampleRate/bpsec]);
Song = zeros([[16*3*SampleRate]+218400-[cadence] 2]);



for Hz = 1:length(FinalMelody)
    
    
    if Type(1)
        if ~isempty(FinalMelody{Hz})
            %         x1 = -.25*sin(2*pi*Frequencies(ismember(Scale,FinalNotes{Hz}))*[0:1/20000:.5]);
            %
            %         x2 = .25*sin(pi*Frequencies(ismember(Scale,FinalNotes{Hz}))*2*[0:1/20000:.5]);
            %         x3 = [sqrt(3)/2]*cos(pi*Frequencies(ismember(Scale,FinalNotes{Hz}))*4*[0:1/20000:.5]);
            %
            %         Song = [Song  x1+x2+x3];
            
            Note = Song([[Hz-1]*cadence ]+1:[[Hz-1]*cadence]+size(MelodyInstru{ismember(Scale,FinalMelody{Hz})},1),1:2)+MelodyInstru{ismember(Scale,FinalMelody{Hz})};
            
            Song([[Hz-1]*cadence ]+1:[[Hz-1]*cadence ]+size(MelodyInstru{ismember(Scale,FinalMelody{Hz})},1),1:2) = Note;
            
            
            if FinalMelody{Hz} == noteseq{2,Hz}
                
                chord1 = rem(find(ismember(Scale,FinalMelody{Hz})) + 2,7);
                if chord1 == 0
                    chord1 = 7;
                end
                Song([[Hz-1]*cadence ]+1:[[Hz-1]*cadence]+size(MelodyInstru{chord1},1),1:2) = Song([[Hz-1]*cadence ]+1:[[Hz-1]*cadence]+size(MelodyInstru{chord1},1),1:2) + MelodyInstru{chord1};
                
                chord2 = rem(find(ismember(Scale,FinalMelody{Hz})) + 4,7);
                if chord2 == 0
                    chord2 = 7;
                end
                Song([[Hz-1]*cadence ]+1:[[Hz-1]*cadence]+size(MelodyInstru{chord2},1),1:2) = Song([[Hz-1]*cadence ]+1:[[Hz-1]*cadence]+size(MelodyInstru{chord2},1),1:2) + MelodyInstru{chord2};
            end

        else
            %Song([[Hz-1]*48000]+1:[Hz]*48000,1:2) = 0;
        end
    end
    
    if Type(2)
        if rem(Hz,ceil(length(FinalMelody)/[size(Song,1)/SampleRate])*ceil(max(cellfun(@(C) size(C,1),BassInstru))/SampleRate)) == 1
            
            if ~isempty(find(cellfun(@(c) ~isempty(c),FinalBass(Hz:min([Hz+3 length(FinalBass)]))),1,'first'))
            LongNote = BassInstru{find(cellfun(@(c) ~isempty(c),FinalBass(Hz:min([Hz+3 length(FinalBass)]))),1,'first')};
            
            
            Note = Song([[Hz-1]*cadence ]+1:[[Hz-1]*cadence]+size(LongNote,1),1:2)+LongNote;
            Song([[Hz-1]*cadence ]+1:[[Hz-1]*cadence ]+size(LongNote,1),1:2) = Note;
            
%             if FinalMelody{Hz} == noteseq{2,Hz}
%             
%             end
            
            
            end
        end
    end
end

%sound(Song,SampleRate)

% xticks(linspace(1,size(Topo{1},2),49))
% xticklabels(0:48)
% 
% legend({'Time' 'Melody' ' ' 'Chord' ' ' 'Bass' ' '})
% timer = plot([1 1],fig.Children(2).YLim,'Color',[.5 .5 .5],'LineWidth',4);
% uistack(timer,'bottom') 
% legend({'Time' 'Melody' ' ' 'Chord' ' ' 'Bass' ' '})
% tic
% a = toc;
% 
% while a < 16*3
%     timer.XData = [[a/48]*size(Topo{1},2) [a/48]*size(Topo{1},2)];
%     pause(.001);
%     a = toc;
% end



audiowrite(['C:\Users\User\Desktop\Music\' num2str(Image) '.wav'],Song,SampleRate)



FinalMelody = array2table(FinalMelody');


end