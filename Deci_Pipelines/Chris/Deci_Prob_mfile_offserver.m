%  Welcome to Deci
% A fluid pipeline to analyze EEG Data with ease
% Written and Compiled by @John Nguyen of Robert Reinhart Lab on 10/15/2018, contact: Akida@bu.edu
% 
% 0. MainMenu

addpath(genpath('/projectnb/crc-nak/chartove/Julia/Deci'))
addpath(genpath('/projectnb/crc-nak/chartove/Julia/OurFieldTrip'))

Deci = [];
Deci.Folder.Raw         = ['G:\Julia\New_Datasets\3_9_2020'];              % Raw Data Folder Directory
Deci.SubjectList        = 'gui';                                                % Cell Array of strings, 'all' or 'gui'
Deci.Step               = 5;                                            % Which Step to implement 1-TD, 2-PP, 3-Art, 4-Analysis, 5-Plot
Deci.Proceed            = false;                                                                  % Continue to next step when current step is done?
Deci.PCom               = true;                                                       % Activates Parallel Computing for PP and Analysis only
Deci.Folder.Version     = ['G:\Julia\ProcessedData'];        % Output Folder Directory

%% 1. Trial Definitions

Deci.DT.Type = 'Prob_Exp2';
Deci.DT.Starts     = {5};
Deci.DT.Ends       = {6};         %Cell Array of Markers for ETs End.
Deci.DT.Markers    = {[14 15 16] [41 42] [51 52] [31 34] [37 38] [22 23]};
Deci.DT.Locks      = [20];
Deci.DT.Toi        = [-4 5];

%MrkMake('C:\Users\User\Desktop\Prob_Stim\RawData',{{'111' '113' '115'} {'112' '114' '116'}},{'110' '117'},[-1/1000 1/1000])

Deci.DT.Block.Markers = {[100]};  
%% 2. PreProcessing Steps
Deci.PP.filter.bpfreq = [.5 100];
Deci.PP.filter.bpfilter = 'yes';
Deci.PP.filter.detrend = 'yes';

Deci.ICA.do = true;
%% 3. Artifact Rejection
Deci.Art.do = true;
Deci.Art.crittoilim = [-.5 1.5];   
Deci.Art.AddComponents =  true;

%Deci.Art.More.channel = {'FCz'};

%% 4. Analysis
Deci.Analysis.Laplace       = false;    
Deci.Analysis.LaplaceFile = ['standard_1020.elc'];
% Uses Captrak .bvct file with same name in RawData Folder
Deci.Analysis.Locks         = [1];                                                          % Which Lock to Analyze        
Deci.Analysis.LocksTitle    = {'Fdb Onset'};                      % Folder Title to save each Lock as
Deci.Analysis.DownSample    = 500;                                              % Downsample Pre-Analysis

% Rsp_Opt = 300; % Use Prob_Exp for trial defs
% Rsp_Wst = 301;
% 
% Fdb_Correct = 22;
% Fdb_Incorrect = 23;

Deci.Analysis.Conditions    = {[22] [23]};
Deci.Analysis.CondTitle     = {'Cor' 'Inc'}; 

Deci.Analysis.Channels = {'FCz'};
Deci.Analysis.Toi           = [-1 2];                                                  % Time Range to save
Deci.Analysis.Toilim        = [-2 3];


Deci.Analysis.Connectivity.Sets = {{{'FCz'} {'FCz'} 'theta' 'beta' 'mi'}};
Deci.Analysis.Connectivity.DownSample = 100;
Deci.Analysis.Connectivity.Zscore.do = false;
Deci.Analysis.Connectivity.Zscore.Runs = 500;
Deci.Analysis.Connectivity.list = [true];
Deci.Analysis.Connectivity.Functions = {'dc_connectivity'};
Deci.Analysis.Connectivity.Params = {{Deci.Analysis.Connectivity}};

%% 5. Plotting

%%  6. Run

Deci.Run.Behavior = false;
Deci.Run.Freq =false;
Deci.Run.ERP =false;
Deci.Run.Extra = true;

if Deci.Step == 4
        Deci_Wavelet;
        
        Deci.Analysis.Freq.method        = 'wavelet';                                                  % Currently only uses 'wavelet' and 'hilbert'
        Deci.Analysis.Freq.foi           = exp(linspace(log(2),log(40),40));                           % Frequency of Interest
        Deci.Analysis.Freq.width         = exp(linspace(log(3),log(13),40));                                                         % Width
        Deci.Analysis.Freq.gwidth        = 4;                                                          % Gwidth

        Deci.Analysis.ERP.do  = false;
        Deci.Analysis.Freq.do  = true;
        Deci.Analysis.Extra.do  = false;
        Deci.Analysis.Connectivity.do = true;
        
        Deci_Backend(Deci);
        Deci = rmfield(Deci,'Analysis');
else
    
    Deci.Plot.BslRef ='Fdb Onset';
    Deci.Plot.Lock = 'Fdb Onset';
    
    Deci.Plot.GrandAverage = true;
    
    Deci.Analysis.CondTitle     = {'Cor' 'Inc'}; 
    
    Deci.Plot.IndexTitle = {'Cor' 'Inc'};
    Deci.Plot.Math         = {'x2-x1'};           % Condition Math done after Bsl, Condition Indexes are appended on.
    Deci.Plot.Draw         = {[1:2] [3]};                   % Cell array of Condition Index for each figure
    Deci.Plot.Figures      = [true true];                             % Which figure to plot currently
    Deci.Plot.Title        = {'All Trials' 'Subtraction'};          % Title for each figure
    Deci.Plot.Subtitle     = {{'Cor' 'Inc'} {'Inc - Cor'}};     % Cell array of strings of subtitles for each Condition
    
    
    Deci.Plot.Topo.do    =false;
    Deci.Plot.Topo.Foi     = [3 8];                   % Frequency of Interest
    Deci.Plot.Topo.Toi     = [];                   % Time of Interest
    Deci.Plot.Topo.Channel = ['Reinhart-All'];                   % Channel of Interest
    
    Deci.Plot.Square.do  =true;
    Deci.Plot.Square.Foi     = [3 inf];                 % Frequency of Interest
    Deci.Plot.Square.Toi     = [];                   % Time of Interest
    Deci.Plot.Square.Channel = [{'Cz'}];                % Channel of Interest
    
    Deci.Plot.Wire.do    =false;
    Deci.Plot.Wire.Foi     = [3 8];           % Frequency of Interest
    Deci.Plot.Wire.Toi     = [];                   % Time of Interest
    Deci.Plot.Wire.Channel = [{'Cz'}];              % Channel of Interest
    
    Deci.Plot.Bar.do    =false;
    Deci.Plot.Bar.Foi     =  [3 8];                 % Frequency of Interest
    Deci.Plot.Bar.Toi     = [.5 1];                   % Time of Interest
    Deci.Plot.Bar.Channel =  [{'Cz'}];
    
    Deci.Plot.Freq.Type    = 'TotalPower'; 
    
    if  strcmpi(Deci.Plot.BslRef, 'Response Onset')
        Deci.Plot.Bsl     = [-.5 -.3];
    else
        Deci.Plot.Bsl     = [-.5 0];
    end
    
    if  Deci.Plot.Square.do  
    Omnibus = true;
    
    Deci.Plot.Stat.do = true;
    if Omnibus == true
        Deci.Plot.Stat.Comp = 'Bsl';
        Deci.Plot.Stat.alpha = .01;
        Deci.Plot.Stat.correctm = 'fdr';
        Deci.Plot.BslType = 'db'
        Deci.Plot.Stat.FPlots =  false;
    else
        Deci.Plot.Stat.do = true;
        Deci.Plot.Stat.Type = 'Randomize Permutation';
        Deci.Plot.Stat.alpha = .05;
        Deci.Plot.Stat.correctm = 'fdr';
    end
    Deci.Plot.Stat.FPlots =  false;
    end
    
    

    
    if Deci.Plot.Topo.do
    Perma = true;
    
    if Perma == true
        Deci.Plot.Stat.Type = 'Randomize Permutation';
        Deci.Plot.Stat.alpha = .05;
        Deci.Plot.Stat.FPlots =  false;
    end
    end
    
    if Deci.Plot.Wire.do
        Deci.Plot.Stat.do = true;
        Deci.Plot.Stat.Type = 'Anova/T-test';
        Deci.Plot.Stat.alpha = .05;
        Deci.Plot.Stat.FPlots =  false;
        Deci.Plot.Stat.correctm = 'no';
        Deci.Plot.BslType = 'relchange';
    end
    
    if Deci.Plot.Bar.do
        Deci.Plot.Stat.do = true;
        Deci.Plot.Stat.Type = 'Anova/T-test';
        Deci.Plot.Stat.alpha = .05;
        Deci.Plot.Stat.FPlots =  false;
        Deci.Plot.Stat.correctm = 'fdr';
        Deci.Plot.BslType = 'relchange';
    end
    
    
Deci.Plot.Extra.Conn.List = Deci.Analysis.Connectivity.Sets([true]);
Deci.Plot.Extra.Conn.FL_FH.do =  true;
Deci.Plot.Extra.Conn.FL_FH.toi = [];

Deci.Plot.Extra.Conn.FL_time.do =  false;
Deci.Plot.Extra.Conn.FH_time.do =  false;

Deci.Plot.Extra.Conn.CL_time.do =  false;
Deci.Plot.Extra.Conn.CH_time.do =  false;

Deci.Plot.Extra.List = [true];
Deci.Plot.Extra.Functions = { 'Plottor_Conn'};
Deci.Plot.Extra.Params = {{Deci.Plot.Extra.Conn}};
    
Deci_Backend(Deci);
end


