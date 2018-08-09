% ----------Welcome to Deci-------
% a fluid pipeline to analyze EEG Data with ease.

% Deci uses Fieldtrip (Robert Oostenveld) as the backbone of its analysis with
% a few codes modified so it may not be compatible with other versions of Fieldtrip.

% Written and Compiled by @John Nguyen of Robert Reinhart Lab on 7/17/2018, contact: Akida@bu.edu.


% --------Disclaimer-----------
%Deci is dwarf version of a large code, Decima, compiled recently so it
% may come with some bugs and error, please contact me anytime.
% Additionally, Deci may only work with newer version of MATLAB

%---------Manual-------

%To Start, Open MainMenu and go the MainMenu Section

% In order to have a good experience, I suggest turn on 'Sections' in Code Folding in Preferences (under Home Tab)
% You'll then be able press [Ctrl and =] to Collapse all Sections to a easier-to-read format

%Deci has 6 Stages, each segmented by it's own section with it's own list of parameters
%MainMenu - Information used to find Raw Files director and Output Directory

%To start, You'll need a Raw Data Folder with your Data with .Folder.Raw
%Types of Raw formats compatible are [.dat, .eeg, .bdf, .cnt]

%Specify your output folder as Folder.Version

%When MainMenu runs, a snapshot of parameter will be saved as "Parameter.txt" in your .Version Directory
%Deci will ask you which steps you want to start on:



%Step 1:DefineTrial
%Deci comes with a rudimentary auto trial definiton, however this comes with a few assumptions.
% You want every trial to be same beg and end, and therefore the same length
% There are unique markers signifying the start and end of your trials
% Any one of your Start Code will appear first in the Marker file before a End Code appears
% There are Equal Amounts of Start-Code Markers and End-Marker Code found in your data.
% There are Markers between your start/end codes specify the type of trial that has occured
% There is a specific marker you will timelock to and set your time axis as 0

%If you find that the above assumptions conflict with your needs,
%you'll have to refer to Fieldtrip's ft_definetrial tutorial to
%create your cfg data

% .Markers,.Conditions,.Locks works respectively, ie: the first
% element of the .Conditions correspond to the first element of
% .Markers

% If you find the need that there are two types of trials that
% you want to be considered as the same Condition, use the
% .Condition field
%Step 2: PreProcessing
% Segments data specified by DefineTrial
%Currently does:
%Offline Rereference of Data, ex: TP10 was online mastoid reference so TP10:TP9 is the rereference of all channels over both mastoid averaged
%Ocular reference, saving new channel as the name of first channel, ie TVEOG:BVEOG averages to become TVEOG
%High-bandpass
%Demean/baseline correction

%Step 3: Artifacts
%Individual Component Analysis
% ICA will directly modify PreProc Data so if you run ICA twice, it might break,
% because it cannot decompose data again into another 20 Components.
% If you wish to revisit ICA, Redo PreProc Data
%Trial Rejection
% Both Eye and Muscle Artifacts. Check code for additional parameters
% Trial Artifacts are not removed from Preproc Data!!!
% Instead, bad trials are tagged and you can specify in Analysis whether or not you want them removed.
%Step 4: Analysis
%Currently Deci is able to compute:
%Frequency data: TotalPower and ITPC [only through wavelet decomposition]
%ERP Data
%During Analysis, Deci will collapse trial information so you'll need to
%specify whether or not you want bad trial removed prior to running.

%Step 5:
%Currently Deci is able to compute:
%Frequency data
%Power-Related Potential, an ReinhartLab analysis that
%looks at correlation between ERP and Freq data.
%% MainMenu
Deci                    = [];
Deci.File = mfilename;
Deci.Folder.Raw         = ['C:\Users\User\Desktop\John\RawData\TimeEst EEG-PS_all_new'];        %Folder with Raw Data files
Deci.SubjectList        = ['all'];                                      %Cell Array of Subjects Files or 'all' to specify all in Deci.Folder.Raw or 'gui' to choose.

Deci.Folder.Version     = ['C:\Users\User\Desktop\John\ProcessedData\TimeEst'];     % This is your Output Directory, it does not have to exist yet

Deci.Layout.Noeye       = 'easycap_rob_noeye.mat';                      % Layout without Ocular Channels, Default is ReinhartLab 1020 system
Deci.Layout.eye         = 'easycap_rob_binocular.mat';                        % Layout with Ocular Channels

Deci.Step = input('Start from which step? DefineTrial = 1, PreProcess = 2, Artifacts = 3, Analysis = 4, Plots = 5 >>');

Deci.Proceed            = 0;                                             %0 or 1, Do you want Deci to automatically proceed to next steps? Suggest 0 for first-time users.
Deci.Debug = 1;
Deci = Checkor(Deci);

%% 1. Define Trial
if Deci.Step <= 1
    
    %Experimental Trials are specfic to the experiment.
    %Collected Trials are specfic to the Trial Definition and may exceed the time from Experimental Trials
    
    Deci.DT.Type = 'timeestfun2';                       % Define trial via 'Manual' with definitions below
    % or file name with output trl [Store them in Xtra_Fun]
    
    %Manual Definition
    Deci.DT.Starts     = {10 11};       % Markers Represent Start of Experimental* Trial
    Deci.DT.Ends       = {200 201 202}; % Markers Represent End of   Experimental* Trial
    
    Deci.DT.Markers    = {[10] [11]};   % 1xConds Cell Array where specified Marker(s) within Experimental* Trial indicates Collection* for a specifc Cond
    Deci.DT.Conditions = [1 2];         % 1xConds Vector where each Cond is a numerical value to show if Conds should be Collapsed
    Deci.DT.Locks      = [10 11];       % 1xConds Vector for each Collected*
    
    Deci.DT.Toi        = [-2 2];        % [Beg End], Time of Interest for Collected* Trial respective to it's Lock Marker
    
    % if you get a "cannot read data after the end of the file" error, this
    % is because the time length request after the Lock is further than the
    % file's end. Suggestion is tp remove the last Trial and make sure to
    % leave data recording to a few seconds after last trial finishes.
    
    DefineTrialor(Deci);
    if ~Deci.Proceed; return; end
end
%% 2. PreProcessing Steps
% PreProcess
if Deci.Step <= 2
    %Leave any option empty ,[], if not desired
    
    Deci.PP.ScalingFactor       = [];                       % If the data needs to scale
    
    %ReinhartLab Settings, Template Format
    Deci.PP.Imp                = 'TP10:TP9';                   % 'Impicit Reference?[Implicit:RefChannels without Implicit]'
    Deci.PP.Ocu                = 'TVEOG:BVEOG,LHEOG:RHEOG';    % 'Ocular Reref?[Ref-Channels without Ref]'
    
    %Rereference channel to average of Implicit Electrode and Offline Reference
    %Deci.PP.Imp                 = 'RM:LM';                      % 'Implicit Reference?[Implicit:RefChannels without Implicit]
    %Ocular channel reference for average of Horizontal/Vertial average. Takes the naming of the Ref channel.
    %Deci.PP.Ocu                 = [];                           % 'Ocular Reref?[Ref:Channels without Ref]'
    
    Deci.PP.HBP                 = [];                           % 'High-bandpass filter, in Hz.'
    Deci.PP.Demean              = [-.5 -.2];                      % 'Baseline Window for Correction'
    Deci.PP.Repair              = [];                           % 'View and Repair Channels'
    Deci.PP.Repair.Type         = 'Auto';                        %'Manual' or 'Auto'
    Deci.PP.Repair.Auto         = {{} {} {'Pz'} {} {'FC4'} {} {} {} {} {} {} {} {} {'FCz'} {} {} {} {}};      %If Auto, Nx1 Subject cell array of Mx1 cells of Channels to reject
    
    
    PreProcessor(Deci);
    if ~Deci.Proceed; return; end
end
%% 3. ICA and Artifacts
if Deci.Step <= 3
    
    % Individual Component Analysis, 20 Components.
    Deci.Art.ICA                     = 0; % Whether or Not to do ICA
    
     % Trial Artifacts
    Deci.Art.TR.Toi              = [-.5 -.2];          % Time Range of Interest to look for Artifacts
    Deci.Art.TR.rToi           = [-.5 1];
    
    
    Deci.Art.TR.Eye                  = [];               % Comment the rest of .Eye to not do Eye Artifact Trial Rejection
    Deci.Art.TR.Eye.Interactive      = 0;                % Enable Interactive Trial Artifact Rejection
    Deci.Art.TR.Eye.Chans            = {'BVEOG','RHEOG','AF7','AF8'};    % 1xChan Cell Array of Ocular Eye Channels

    Deci.Art.TR.Muscle               = [];               % Comment the rest of .Muscle to not do Muscle Artifact Trial Rejection
    Deci.Art.TR.Muscle.Interactive   = 0;                % Enable Interactive Trial Artifact Rejection


    Deci.Art.Manual = [];
    Deci.Art.Manual.Toi              = [-.5 -.2]; 
    Deci.Art.Manual.rToi              =[-.5 1];

    Artifactor(Deci);
    if ~Deci.Proceed; return; end
end
%% 4. Analysis
if Deci.Step <= 4
    
     %Deci.Folder.Version      = [Deci.Folder.Version '_nolaplace'];
     %Deci.Folder.Analysis     = [Deci.Folder.Version filesep 'Analysis'];

    Deci.Analysis.ArtifactReject     = 1;                    % Specify whether or not to Apply Trial Rejection
    Deci.Analysis.Channels           = 'all';                % Cell array of channels or 'all'
    Deci.Analysis.Laplace            = 1;                   % Leave .Laplace as 0 to not do a Laplacian Transformation based on file realistic_1005.txt format
    
    Deci.Analysis.Redefine = [];                            % Parameters if a Redefinition was Created
    Deci.Analysis.Redefine.Bsl = [-.5 -.2];
    Deci.Analysis.Redefine.ERPToi = -.5;
    
    % Fourier analysis to collect Trial-collapsed TotalPower and ITPC
    % This .Freq Structure follows exactly ft_freqanalysis but may need to
    % modify due to the Freq.Toi parameter.
    Deci.Analysis.Freq               = [];           % Comment the rest of .Freq to not do Fourier Transformation
    Deci.Analysis.Freq.method        = 'wavelet';    % Currently only uses Wavelet Decomp
    Deci.Analysis.Freq.foi           = exp(linspace(log(2),log(40),40));       % Frequency of Interest
    Deci.Analysis.Freq.width         = 7 ;           % Width
    Deci.Analysis.Freq.gwidth        = 4;            % Gwidth
    Deci.Analysis.Freq.Toi           = [-.5 1];       % Time Range
    
    Deci.Analysis.Freq.Redefine = 1;   % 1 for Both BSL and Freq, or 2 just BSL
    
    %Event Related Potential Analysis
    %Note that entire Time range will be analyzed.
    Deci.Analysis.ERP                = 1;
    
    Deci.Analysis.EvokedPower        = 0; %Not Available
    
    Analyzor(Deci);
    if ~Deci.Proceed; return; end
end
%% 5. Plotting
if Deci.Step <= 5
    
    Deci.Folder.Version      = [Deci.Folder.Version '_nolaplace'];
    Deci.Folder.Analysis     = [Deci.Folder.Version filesep 'Analysis'];
    Deci.Folder.Plot         = [Deci.Folder.Version filesep 'Plot'];
    
    Deci.Plot.GA = 1;                       % 0 for Individual Plots, 1 for GrandAverage
    Deci.Plot.Scale = 'log';
    
    Deci.Plot.Save = [];
    Deci.Plot.Save.Format = 'jpg';
    
    Deci.Plot.Math = [];
    Deci.Plot.Math.Form = {'[x1+x3]','[x2+x4]','x6-x5'};                    % Not yet Available
    Deci.Plot.Math.Type = 1;                    % 0 = no maths, 1 = only math plots, 2 = both math and single plots
    
    Deci.Plot.Toi = [-.2 .5];                % Time Range of Interest
    Deci.Plot.Toi = [];
    
    Deci.Plot.Freq = [];
    Deci.Plot.Freq.Plot = 1;                % 1 or 0 to do Freq Plots
    Deci.Plot.Freq.Type = 'TotalPower';     % Use 'TotalPower','ITPC' or 'EvokedPower'
    Deci.Plot.Freq.Foi = [4 8];             % Leave empty to reach all
    Deci.Plot.Freq.Foi = [3 inf];
    
    Deci.Plot.Freq.Roi = 'maxabs';          % use 'maxmin','maxabs or [min max] for settin
    %Deci.Plot.Freq.Roi = [-1 3];
    Deci.Plot.Freq.BslType = 'db';    % Use 'absolute','relative','relchange','db','normchange'
    Deci.Plot.Freq.Bsl = [-.2 0];           % Baseline Correction Time
    Deci.Plot.Freq.Channel = 'all';              % Cell array of channels or 'all'
    %Deci.Plot.Freq.Channel = {'FCz'};
    
    Deci.Plot.Freq.Topo = 1;
    Deci.Plot.Freq.Square = 1;
     
     Deci.Plot.Freq.Wires = [];
%     Deci.Plot.Freq.Wires.Collapse = 'Time';        %Which Dim to Collapse 'Time' or 'Freq'
%     Deci.Plot.Freq.Wires.avgoverfreq = 'no';
%     Deci.Plot.Freq.Wires.avgovertime= 'yes';
%     Deci.Plot.Freq.Wires.frequency = [5 8];
%     Deci.Plot.Freq.Wires.latency = [.1  .5];
%     Deci.Plot.Freq.Wires.nanmean = 'no';
%     
    
    Deci.Plot.ERP = []; %Not Available to plot
    %     Deci.Plot.ERP.Channel = 'all';              % Cell array of channels or 'all'
    
    
    Deci.Plot.PRP = [];                     %Requires .Freq/.ERP to be filled out
    %Power related Potential
    %Plot 1 = Collapsed-Frequency/Channel Freq.Type Data
    %plot 2 = Collapsed-Channel           ERP Data
    %plot 3 = Condition-Specific Timepoint by Timepoint Freq-ERP coorelation across subjects
    %plot 4 = Collapsed-Time/Freq/Channel Condition-Specific Subject by Subject plotted along Freq-YAxis, ERP-XAxis
    %Changes .GA into 0
    %     Deci.Plot.PRP.label = 1;                %Whether or Not to add subject info on plot 4
    %     Deci.Plot.PRP.ScatterToi = [.3 .7];     %Time range for plot 4
    
    Plottor(Deci);
    if ~Deci.Proceed; return; end
end