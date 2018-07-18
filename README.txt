% ----------Welcome to Deci------- 
% a fluid pipeline to analyze EEG Data with ease.
 
% Make sure OurFieldTrip/Deci are in your matlab paths.
% Deci uses Fieldtrip (Robert Oostenveld) as the backbone of its analysis with
% a few codes modified so it may not be compatible with other versions of Fieldtrip.
 
 
% Written and Compiled by @John Nguyen of Robert Reinhart Lab on 7/17/2018, contact: Akida@bu.edu.
 
% In order to have a good experience, I suggest turn on 'Sections' in Code Folding in Preferences (under Home Tab)
% You'll then be able press [Ctrl and =] to Collapse all Sections to a easier-to-read format
 
% --------Disclaimer----------- 
%Deci is dwarf version of a large code, Decima, compiled recently so it
% may come with some bugs and error, please contact me anytime.
% Additionally, Deci may only work with newer version of MATLAB
 
%---------Manual-------
 
%Deci has 6 Stages
    %MainMenu - Information used to find Raw Files director and Output Directory
        %Types of Raw formats compatible are [.dat, .eeg, .bdf, .cnt]
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
