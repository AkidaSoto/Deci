
function Plottor_Freq(Deci)

%% Deci Checks
Dims = {'Topo' 'Square' 'Wire' 'Bar'};
[Deci, info] = dc_plotcheck(Deci,Dims);
info.isfreq = true;

%% Load
disp('----------------------');
display(' ')

display(['Plotting ' Deci.Plot.Freq.Type]);

info.extension = ['Freq_' Deci.Plot.Freq.Type];
info.parameter = 'powspctrm';
info.variable = 'freq';

[Subjects,info] =  dc_plotload(Deci,info);

%% Baseline Correction
display(' ');
display(['Using Lock: ' Deci.Plot.Lock]); 
display(['Using Ref: ' Deci.Plot.BslRef ' at times ' strrep(regexprep(num2str(Deci.Plot.Bsl),' +',' '),' ','-')]);

Subjects = dc_plotbsl(Deci,Subjects,info);

%% Math
if ~isempty(Deci.Plot.Math)
    Subjects = dc_math(Deci,Subjects,info);
end

%% Data Management
if size(Subjects,1) == 1
    Deci.Plot.GrandAverage = false;
end

if Deci.Plot.GrandAverage
    if any(~isnan(info.trllen))
        info.trlstd = nanstd(info.trllen,[],1);
        info.trllen = nanmean(info.trllen,1);
    end
    
    if any(~isnan(info.lockers))
        info.lockersstd = nanstd(info.lockers,[],1);
        info.lockers = nanmean(info.lockers,1);
    end
end

%% Plot Types
if Deci.Plot.Topo.do
    dc_plottopo(Deci,Subjects,info);
end

if Deci.Plot.Square.do
    dc_plotsquare(Deci,Subjects,info);
end

if Deci.Plot.Wire.do
    dc_plotwire(Deci,Subjects,info)
end

if Deci.Plot.Bar.do
    dc_plotbar(Deci,Subjects,info)
end
end
