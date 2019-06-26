%% *Welcome to Deci*
% *A fluid pipeline to analyze EEG Data with ease*
% 
% Written and Compiled by @John Nguyen of Robert Reinhart Lab on 10/15/2018, 
% contact: Akida@bu.edu
%% 0. MainMenu

Deci = [];
Deci.Folder.Raw         = ['/Users/researcher/Documents/SRTData/SRT2Data/013BU'];    
Deci.SubjectList = 'gui';  
Deci.Step               = 2;
Deci.Proceed            = false;  
Deci.PCom               = false;                           % Activates Parallel Computing (if availible)
Deci.Folder.Version     = ['/Users/researcher/Documents/SRTData/SRT2Data/013BU'];    
%% 1. Trial Definitions

    Deci.DT.Type = 'Manual';             
    Deci.DT.Starts     = {11 12};         %Cell Array of Markers for ETs Start.
    Deci.DT.Ends       = {30};         %Cell Array of Markers for ETs End.
    Deci.DT.Markers    = {[11 12] [21 22 23 24] [26 27 28 29] [42 43]};   
    Deci.DT.Locks      = [52 53];
    Deci.DT.Toi        = [-2 3]; 
    
    Deci.DT.Block.Start = {10};
    Deci.DT.Block.End = {13};
    Deci.DT.Block.Markers = [];
%% 2. PreProcessing Steps

Deci.PP.Demean              = [];        % [Beg End], Secs Relative to .DT.Locks, Baseline Window for Correction.
disp('Preproc done')
%% 3. Artifact Rejection

Deci.Art.crittoilim = [-.6 .8]; 
Deci.Art.bpf = [1 100];
%% 4. Analysis
% General Set-Up

Deci.Analysis.Laplace            =true;                   % Leave .Laplace as 0 to not do a Laplacian Transformation based on file realistic_1005.txt format
Deci.Analysis.Locks = [1 2]; 
Deci.Analysis.DownSample = 500;
% Time-Frequency Analysis

Deci.Analysis.Freq.do = true;
Deci.Analysis.Freq.method        = 'wavelet';    % Currently only uses Wavelet Decomp
Deci.Analysis.Freq.foi           = exp(linspace(log(2),log(100),40));       % Frequency of Interest
Deci.Analysis.Freq.width         = 7;           % Width
Deci.Analysis.Freq.gwidth        = 3;            % Gwidth
Deci.Analysis.Freq.Toi           = [-.6 .6];       % Time Range
Deci.Analysis.Freq.Toilim        = [-2 3]; 


Deci.Analysis.Var =false;
Deci.Analysis.Clean = true;
% ERP

Deci.Analysis.ERP.do  = false;
Deci.Analysis.ERP.Toi = [-.5 1.5];
%% 5. Plotting

    Deci.Plot.Lock = '2';
    Deci.Plot.Conditions = {[11 42] [12 42] [11 43] [12 43]};
    Deci.Plot.Var = false;
    Deci.Plot.GA = false;
    Deci.Plot.Scale = 'linear';
% Math

%Deci.Plot.Math.Form = {'[x4+x2]/2','[x7+x5]/2','[x2+x5]/2','[x4+x7]/2' '[x17-x18]' '[x19-x20]'...
 %                      '[x2-x4]','[x5-x7]','[x4-x7]','[x2-x5]' '[[x2+x3+x4+x13+x14+x15]/6]' '[[x5+x6+x7+x10+x11+x12]/6]' '[x27-x28]' ...
  %                     '[x4+x2+x12+x10]/2','[x7+x5+x15+x13]/2','[x2+x5+x10+x13]/2','[x4+x7+x12+x15]/2' '[x2+x3+x4]/3' '[x13+x14+x15]/3' '[x5+x6+x7]/3' '[x10+x11+x12]/3'  ...
   %                     '[x2+x10]' '[x3+x11]/2' '[x4+x12]/2' '[x5+x13]/2' '[x6+x14]/2' '[x7+x15]/2'};   
                    
                        
                      Deci.Plot.Math.Form = {'[x1 - x2]'};    
%Deci.Plot.Math.Type = 0;

% Conditions

Deci.Plot.Draw = {[1 2] [5]};
Deci.Plot.Title = {'Seq and Rand' 'Seq - Rand'};
Deci.Plot.Subtitle = {{'Seq' 'Rand'} {'Seq - Rand'}};
Deci.Plot.Figures = [true true ];
% Time-Frequency

%Deci.Plot.Freq.Foi = [];
Deci.Plot.Freq.Foi = [2 inf];
Deci.Plot.Freq.Toi = [];
Deci.Plot.Freq.Channel = {'FCz'};
Deci.Plot.Freq.Type = 'TotalPower';     % Use 'TotalPower','ITPC' or 'EvokedPower'

Deci.Plot.Freq.Bsl = [-.6 -.4];           % Baseline Correction Time
Deci.Plot.Freq.Roi = 'maxabs';          % use 'maxmin','maxabs or [min max] to set.
Deci.Plot.Freq.BslType ='db';    % Use 'absolute','relative','relchange','db','normchange'

Deci.Plot.Freq.Topo =true ;
Deci.Plot.Freq.Square =true ;

Deci.Plot.Freq.Wires.avg =  'freq';
Deci.Plot.Freq.Wires.errorbars = true ;

Deci.Plot.CFC.chanlow =  {'FCz'};
Deci.Plot.CFC.chanhigh = {'CP3'};
Deci.Plot.CFC.latencyhigh = [-.2 1]; %no implementation yet to check for edge artifacts so be careful
Deci.Plot.CFC.latencylow = [-.2 1];
Deci.Plot.CFC.freqhigh = 'theta';
Deci.Plot.CFC.freqlow = 'beta';

Deci.Plot.CFC.timebin = 6;
Deci.Plot.CFC.methods = {'mi' 'plv'};  % can only currently handle '' or 'mi'(pac)
Deci.Plot.CFC.errorbars = true ;
Deci.Plot.CFC.Topo = 0;
Deci.Plot.CFC.Square = 0;
Deci.Plot.CFC.Hist = 0;
Deci.Plot.CFC.Wire = 1;
Deci.Plot.CFC.Roi = 'maxmin';          % Power Range of Interest

%% 
% *ERP*

Deci.Plot.ERP.Toi = [0 .8];
Deci.Plot.ERP.Channel = {'FCz'};

Deci.Plot.ERP.Bsl = [-.5 -.2];           % Baseline Correction Time

Deci.Plot.ERP.Topo =0;
Deci.Plot.ERP.Wires =1;
Deci.Plot.ERP.errorbars =  true;
% Behavioral

Deci.Plot.Behv = [];
Deci.Plot.Behv.Source = 'Definition';

% Deci.Plot.Behv.Acc.Total = {[1 3] [1:16]};
% Deci.Plot.Behv.Acc.Subtotal = {[1:8] [9:16]};
% Deci.Plot.Behv.Acc.Title = {'All Trials Percent' 'All Trials Percent'};
% Deci.Plot.Behv.Acc.Subtitle = {'Opt Percent' 'Wst Percent'};
% 
% Deci.Plot.Behv.Acc.Block = [1:6];
% Deci.Plot.Behv.Acc.Collapse.Trial =true;
% Deci.Plot.Behv.Acc.Collapse.Block = true;
% Deci.Plot.Behv.Acc.Collapse.Subject = true;
% Deci.Plot.Behv.Acc.Collapse.Uneven = 'maxlength:nans';
% Deci.Plot.Behv.Acc.Collapse.Movmean =  false;

Deci.Plot.Behv.RT.Draw = {[1] [2]};
Deci.Plot.Behv.RT.Title = {'Seq' 'Rand'};
Deci.Plot.Behv.RT.Subtitle = {'Seq' 'Rand'};
Deci.Plot.Behv.RT.Locks = [1 2];

Deci.Plot.Behv.RT.Block = [];
Deci.Plot.Behv.RT.Collapse.Trial =true;
Deci.Plot.Behv.RT.Collapse.Block = false;
Deci.Plot.Behv.RT.Collapse.Subject = false;
Deci.Plot.Behv.RT.Collapse.Uneven = 'maxlength:nans';
%% 
% * Extra*

% Deci.Plot.Extra.List = true;
% Deci.Plot.Extra.Functions = {@QL};
% Deci.Plot.Extra.Params = {{}};

%% 
% Options

Deci.Run.Behavior = false;
Deci.Run.Freq =true;
Deci.Run.ERP =false;
Deci.Run.CFC = false;
Deci.Run.Extra = false;
%% * 6. Run*
% 

Deci_Backend(Deci);