function dc_Correlation(Deci,params)


%% Deci Checks
Dims = {'Topo' 'Square' 'Wire' 'Bar'};
[Deci, info] = dc_plotcheck(Deci,Dims);
info.isfreq = true;
info.isconn = false;

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
    [Deci, Subjects] = dc_math(Deci,Subjects,info);
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

for draw = 1:length(Deci.Plot.Draw)
    
    for subdraw = 1:length(Deci.Plot.Draw{draw})
       
        
        
        subdata = arrayfun(@(c) c.powspctrm,[Subjects{:,Deci.Plot.Draw{draw}(subdraw)}],'un',0);
        subdata = cat(4,subdata{:});
        
        corrparam = params.Values{draw}(subdraw);
        Corrdata = cell2mat(arrayfun(@(c) corrcoef(c,corrparam,'un',0),subdata));
        
        
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