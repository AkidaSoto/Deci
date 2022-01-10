function Plottor_RawFreq(Deci,Params)

%% Load
cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';

disp('----------------------');
display(' ')

%% Deci Checks
Dims = {'Topo' 'Wire' 'Bar'};
Tois = [];
for Dim = 1:length(Dims)
    if Deci.Plot.(Dims{Dim}).do
        if isequal(Deci.Plot.(Dims{Dim}).Channel,'Reinhart-All')
            Deci.Plot.(Dims{Dim}).Channel = dc_getchans('noeyes');
        end
        
        if isequal(Deci.Plot.(Dims{Dim}).Channel,'Reinhart-All_eyes')
            Deci.Plot.(Dims{Dim}).Channel = dc_getchans('all');
        end
        
        if isequal(Deci.Plot.(Dims{Dim}).Channel,'PD-All')
            Deci.Plot.(Dims{Dim}).Channel = dc_getchans('PD-All');
        end
        
        Tois{Dim} = Deci.Plot.(Dims{Dim}).Toi;
        Chois{Dim} = Deci.Plot.(Dims{Dim}).Channel;
    end
end

if ~any(cellfun(@(c) Deci.Plot.(c).do,Dims))
    error('No Figure Type has been checked')
end

Tois = sort([Tois{:}]);

Tois = [Tois(1) Tois(end)];
Chois = unique([Chois{:}]);

if length(Deci.Plot.Topo.Channel) == 1
    Deci.Plot.Topo.do = 0;
end

%% Load
disp('----------------------');
display(' ')

display(['Plotting ERP']);

for  subject_list = 1:length(Deci.SubjectList)
    
    display(['Loading Plottor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
    
    for Conditions = 1:length(Deci.Plot.CondTitle)
        for Channel = 1:length(Chois)
            
            load([Deci.Folder.Analysis filesep 'Extra/RawFreq' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.CondTitle{Conditions} filesep Chois{Channel} '.mat'],'rawfreq');

            Chans{subject_list,Conditions,Channel} = rawfreq;
            Chans{subject_list,Conditions,Channel}.label = Chois(Channel);
           
            clear rawfreq;
        end 
    end

end

for Channel = 1:size(Chans,3)
    for Conditions = 1:length(Deci.Plot.CondTitle)
        a = figure;
        a.Visible = 'on';
        for  subject_list = 1:length(Deci.SubjectList)
          for block = 1:length(Chans{subject_list,Conditions,Channel}.trialpow)
            scatter(Chans{subject_list,Conditions,Channel}.trialnums{block},Chans{subject_list,Conditions,Channel}.trialpow{block})
            hold on
          end
        end
        
    end
end

end
