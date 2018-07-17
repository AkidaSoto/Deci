%Ctrl and = to Collapse all to a legible format
%turn on Code Folding in Preferences (under Home)
%% MainMenu
Deci = [];
Deci.Folder.Raw          = ['C:\Users\Researcher\Documents\Raw'];
Deci.SubjectList = ['all']; %Cell Array or 'all' to specify all in Deci.Folder.Raw

Deci.Folder.Version = ['C:\Users\Researcher\Documents\test2\']; %Change this to start over in a few folder

Deci.Layout.Noeye = 'easycap_rob_noeye.mat';
Deci.Layout.eye = 'easycap_rob_eye.mat';

%Steps Are DefineTrial = 1, PreProcess = 2, Artifacts = 3, Analysis = 4,
%Plots = 5
Deci.Step = input('Start from which step? DefineTrial = 1, PreProcess = 2, Artifacts = 3, Analysis = 4, Plots = 5 >>');

Deci = Checkor(Deci);
%% 1. Define Trial
if Deci.Step <= 1
    
    Deci.DT.Toi = [-2 2]; 
    Deci.DT.Starts     = {10 11};
    Deci.DT.Ends       = {200 201 202};
    
    Deci.DT.Markers    = {[10] [11]};
    Deci.DT.Conditions = [1 2];
    Deci.DT.Locks      = [10 11];
    
    % if you get a "cannot read data after the end of the file" error, this
    % is because the time length request after the Lock is further than the
    % file's end. Suggestion is remove the last Trial and make sure to
    % leave data recording to a few seconds after last trial finishes.
    
    DefineTrialor(Deci);
end
%% 2. PreProcessing Steps

% PreProcess
if Deci.Step <= 2
    %Leave any option empty ,[], if not desired
    Deci.PP.ScalingFactor = .00819; %If the data needs to scale
    
    %ReinhartLab Settings
    %Deci.PP.Imp = 'TP10:TP9'; %'Impicit Reference?[Implicit:RefChannels-Implicit]'
    %Deci.PP.Ocu  = 'TVEOG:BVEOG,LHEOG:RHEOG'; %'Ocular Reref?[Ref:Channels-Ref]'
    
    Deci.PP.Imp = 'RM:LM'; %'Impicit Reference?[Implicit:RefChannels-Implicit]'
    Deci.PP.Ocu  = []; %'Ocular Reref?[Ref:Channels-Ref]'
        
    Deci.PP.hbp = [];  %'High-bandpass filter'
    Deci.PP.Demean = [-.2 0];  %'Baseline Window?'
    PreProcessor(Deci);
end
%% 3. ICA and Artifacts
if Deci.Step <= 3
   Deci.Art.ICA = 0;
   Deci.Art.TR.Eye  = [];
   Deci.Art.TR.Muscle = [];
   
   Deci.Art.TR.Eye.Interactive      = 0;
   Deci.Art.TR.Eye.Chans            = {'VEM','HEM'};
   Deci.Art.TR.Eye.Toi              = [-.2 1];
   
   Deci.Art.TR.Muscle.Interactive      = 0;
   Deci.Art.TR.Muscle.Toi              = [-.2 1];
   
   Artifactor(Deci);
end
%% 4. Analysis
if Deci.Step <= 4
   
   Deci.Analysis.ArtifactReject     = 1;
   Deci.Analysis.Toi                = [-2 2];
   Deci.Analysis.Channels           = 'all';
   Deci.Analysis.Laplace            = [];
%  Deci.Analysis.Laplace.EyeChans   = {'RHEOG' 'BVEOG'};
   
   Deci.Analysis.Freq = []; 
   Deci.Analysis.Freq.method        = 'wavelet'; %Currently only uses Wavelet Decomp
   Deci.Analysis.Freq.foi           = 1:1:60;
   Deci.Analysis.Freq.width         = 7 ;
   Deci.Analysis.Freq.gwidth        = 4;
   
   Deci.Analysis.ERP                = 1;
   Deci.Analysis.EvokedPower        = 0; %Not Available
   
   Analyzor(Deci);
end
%% 5. Plotting
if Deci.Step <= 5
    
Deci.Plot.GA = 0; % 0 for Individual Plots, 1 for GrandAverage
Deci.Plot.Math = []; % Not yet Available

Deci.Plot.Channel = 'all';
Deci.Plot.Toi = [-.2 1];

Deci.Plot.Freq = [];
Deci.Plot.Freq.Plot = 1;
Deci.Plot.Freq.Type = 'TotalPower'; % Use 'TotalPower','ITPC' or 'EvokedPower'
Deci.Plot.Freq.Foi = [4 8]; %Leave empty to reach all
Deci.Plot.Freq.Roi = 'maxmin';% use 'maxmin','maxabs or [min max] for settin
Deci.Plot.Freq.BslType = 'relative'; % Use 'absolute','relative','relchange','db','normchange'
Deci.Plot.Freq.Bsl = [-.2 0]; %Baseline Correction Time  

Deci.Plot.ERP = 0; %Not Available

Deci.Plot.PRP = []; %Requires .Freq to be filled out
                    %Also Requires Deci.Plot.GA == 0; for 3rd plot
% Deci.Plot.GA = 0;                                      
% Deci.Plot.PRP.label = 0;
% Deci.Plot.PRP.ScatterToi = [.3 .7];

Plottor(Deci);
end