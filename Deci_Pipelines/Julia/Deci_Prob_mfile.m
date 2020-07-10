%  Welcome to Deci
% A fluid pipeline to analyze EEG Data with ease
% Written and Compiled by @John Nguyen of Robert Reinhart Lab on 10/15/2018, contact: Akida@bu.edu
% 
% 0. MainMenu

addpath(genpath('/projectnb/crc-nak/chartove/Julia/Deci'))
addpath(genpath('/projectnb/crc-nak/chartove/Julia/OurFieldTrip'))

Deci = [];
Deci.Folder.Raw         = ['/projectnb/crc-nak/chartove/Julia/data/data_3_9_20'];              % Raw Data Folder Directory
Deci.SubjectList        = 'gui';                                                % Cell Array of strings, 'all' or 'gui'
Deci.Step               = 2;                                            % Which Step to implement 1-TD, 2-PP, 3-Art, 4-Analysis, 5-Plot
Deci.Proceed            = false;                                                                  % Continue to next step when current step is done?
Deci.PCom               = false;                                                       % Activates Parallel Computing for PP and Analysis only
Deci.Folder.Version     = ['/projectnb/crc-nak/chartove/Julia/data/ProcessedData'];        % Output Folder Directory

%% 1. Trial Definitions

Deci.DT.Type = 'Prob_Exp';
Deci.DT.Starts     = {5};
Deci.DT.Ends       = {6};         %Cell Array of Markers for ETs End.
Deci.DT.Markers    = {[14 15 16] [41 42] [51 52] [31 34] [37 38] [22 23]};
Deci.DT.Locks      = [12 30 20];
Deci.DT.Toi        = [-2 3];

%MrkMake('C:\Users\User\Desktop\Prob_Stim\RawData',{{'111' '113' '115'} {'112' '114' '116'}},{'110' '117'},[-1/1000 1/1000])

Deci.DT.Block.Markers = {[100]};  
%%

% 2. PreProcessing Steps
Deci.PP.filter.bpfreq = [1 50];
Deci.PP.filter.bpfilter = 'yes';
Deci.PP.filter.detrend = 'yes';

Deci.ICA.do = true;

% 3. Artifact Rejection
Deci.Art.do = true;
Deci.Art.crittoilim = [-.5 1.5];   
Deci.Art.AddComponents =  true;

Deci.Art.More.channel = {'FCz'};

% 4. Analysis


Deci.Analysis.Laplace       =false;                                                               % Uses Captrak .bvct file with same name in RawData Folder
Deci.Analysis.Locks         = [1 2 3];                                                          % Which Lock to Analyze        
Deci.Analysis.LocksTitle    = {'Stim Onset','Response Onset','Fdb Onset'};                      % Folder Title to save each Lock as
Deci.Analysis.DownSample    = 500;                                              % Downsample Pre-Analysis


Deci.Analysis.Conditions    = {[14 300] [14 301] [15 300] [15 301] [16 300] [16 301]};
Deci.Analysis.CondTitle     = {'AB_Cor' 'AB_Inc' 'CD_Cor' 'CD_Inc' 'EF_Cor' 'EF_Inc'}; 

Deci.Analysis.Channels = {'FCz'};
Deci.Analysis.Toi           = [-1 2];                                                  % Time Range to save
Deci.Analysis.Toilim        = [-2 3];


Deci.Analysis.Connectivity.Sets = {{{'FCz'} {'FCz'} 'theta' 'theta' 'mi'}};

Deci.Analysis.Connectivity.toi = [-.5 1.5]; 
Deci.Analysis.Connectivity.DownSample = 100;

Deci.Analysis.Connectivity.list = [true];
Deci.Analysis.Connectivity.Functions = {'dc_connectivity'};
Deci.Analysis.Connectivity.Params = {{Deci.Analysis.Connectivity}};


% Time-Frequency Analysis

% 5. Plotting

%  6. Run

Deci.Run.Behavior = false;
Deci.Run.Freq =false;
Deci.Run.ERP =false;
Deci.Run.Extra = true;

if Deci.Step == 4
        Deci_Wavelet;

        Deci.Analysis.ERP.do  = false;
        Deci.Analysis.Freq.do  = false;
        Deci.Analysis.Extra.do  = false;
        Deci.Analysis.Connectivity.do = true;
        
        Deci_Backend(Deci);
        Deci = rmfield(Deci,'Analysis');
else
    
    Deci.Plot.BslRef = 'Stim Onset';
    Deci.Plot.Lock = 'Stim Onset';
    
    Deci.Plot.GrandAverage = true;
    
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
    

Deci_Backend(Deci);
end


