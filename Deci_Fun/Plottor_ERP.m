function Plottor_ERP(Deci)

cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';


%% Deci Checks

if isequal(Deci.Plot.ERP.Channel,'Reinhart-All')
    Deci.Plot.ERP.Channel = [{'AF3'  } {'AF4'  } {'AF7'  } ...
        {'AF8'  } {'AFz'  } {'C1'   } {'C2'   } {'C3'   } {'C4'   } {'C5'   } ...
        {'C6'   } {'CP1'  } {'CP2'  } {'CP3'  } {'CP4'  } {'CP5'  } {'CP6'  } ...
        {'CPz'  } {'Cz'   } {'F1'   } {'F2'   } {'F3'   } {'F4'   } {'F5'   } ...
        {'F6'   } {'F7'   } {'F8'   } {'FC1'  } {'FC2'  } {'FC3'  } {'FC4'  } ...
        {'FC5'  } {'FC6'  } {'FCz'  } {'FT7'  } {'FT8'  } {'Fz'   } {'O1'   } ...
        {'O2'   } {'Oz'   } {'P1'   } {'P2'   } {'P3'   } {'P4'   } {'P5'   } ...
        {'P6'   } {'P7'   } {'P8'   } {'PO3'  } {'PO4'  } {'PO7'  } {'PO8'  } ...
        {'POz'  } {'Pz'   } {'T7'   } {'T8'   } {'TP10' } {'TP7'  } {'TP8'  } ...
        {'TP9'  } ] ;
end

if length(Deci.Plot.ERP.Channel) == 1
    Deci.Plot.ERP.Topo = 0;
end

%% Load

if ~isempty(Deci.Plot.PRP)
    Deci.Plot.GA = 0;
    
    if isempty(Deci.Plot.PRP.label)
        Deci.Plot.Freq.BslType =  0;
        warning('Cannot Find label for PRP plot, setting as 0');
    end
    
    if isempty(Deci.Plot.PRP.ScatterToi)
        Deci.Plot.Freq.ScatterToi =  [-inf inf];
        warning('Cannot Find ScatterToi for PRPplot, setting as  [-inf inf]');
    end
end

%% Load
for  subject_list = 1:length(Deci.SubjectList)
    
    tic;
    
    ERP.Channels       = CleanDir([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Lock]);
    ERP.Channels       = cellfun(@(c) c(1:end-4),ERP.Channels,'un',0);
    
    if ~ischar(Deci.Plot.Freq.Channel)
        if ~any(ismember(Deci.Plot.ERP.Channel,ERP.Channels))
            error('Freq Channel Parameter has additional channels not found in Data')
        end
        
        ERP.Channels = Deci.Plot.ERP.Channel(ismember(Deci.Plot.ERP.Channel,ERP.Channels));
    end
    
    for Channel = 1:length(ERP.Channels)
        
        display(['Channel ' ERP.Channels{Channel} ' of ' num2str(length(ERP.Channels))])
        
        erp = [];
        load([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list} filesep Deci.Plot.Lock filesep ERP.Channels{Channel}],'erp');

        Chans{Channel} = erp;
        %Chans{Channel}.time =  Chans{Channel}.time(toi);
        Chans{Channel}.label = ERP.Channels(Channel);
        
    end
    
    
    toc;
    
    if length(Chans) > 1
    acfg.appenddim = 'chan';
    Chans = rmfield(ft_appenddata(acfg,Chans{:}),'cfg');
    
    else
        
    Chans = Chans{1};
    
    end
    
    bpcfg.bpfilter = 'yes';
    bpcfg.bpfreq = [.1 15];
    bpcfg.bpfilttype = 'fir';
    Chans = ft_preprocessing(bpcfg,Chans);
    
    for Conditions = 1:length(Deci.Plot.Conditions)
        
        maxt = max(sum(ismember(erp.condinfo{2},Deci.Plot.Conditions{Conditions}),2));
        
        trl = sum(ismember(erp.condinfo{2},Deci.Plot.Conditions{Conditions}),2) == maxt;
        
        tcfg = [];
        tcfg.trials = trl;
        [timelock] = ft_timelockanalysis(tcfg, Chans);
        
        toi = round(timelock.time,4) >= Deci.Plot.ERP.Toi(1) & round(timelock.time,4) <= Deci.Plot.ERP.Toi(2);
        
        tlcfg.latency = [timelock.time(find(toi,1,'first')) timelock.time(find(toi,1,'last'))];
        
        timelock= ft_selectdata(tlcfg,timelock);
        
        TrialCount{subject_list,Conditions} = length(find(trl));
        
        tcfg.baseline = Deci.Plot.ERP.Bsl;
        
        Subjects{subject_list,Conditions} = ft_timelockbaseline(tcfg,timelock);

        
    end
    clear Chans;
    
end




if Deci.Plot.Math.Type ~= 0
    
    
    for cond = 1:length(Deci.Plot.Math.Form)
        for subj = 1:size(Subjects,1)
            scfg.parameter = 'avg';
            scfg.operation = Deci.Plot.Math.Form{cond};
            MathData{subj} = ft_math(scfg,Subjects{subj,:});
            
        end
        
        Subjects(:,length(Deci.Plot.Conditions)+cond) = MathData;
        TrialCount(:,length(Deci.Plot.Conditions)+cond) = num2cell(nan(size(Subjects,1),1));
        
        
    end
    

end

mkdir([Deci.Folder.Plot filesep 'ERP'])
save([Deci.Folder.Plot filesep 'ERP' filesep 'Subjects'],'Subjects')

if Deci.Plot.GA
    
    for conds = 1:size(Subjects,2)
        facfg.parameter =  'avg';
        facfg.type = 'mean';
        
        ERPData{conds} = rmfield(ft_timelockgrandaverage(facfg,Subjects{:,conds}),'cfg');
        
        TotalCount{conds} = mean([TrialCount{:,conds}]);
        
        facfg.type = 'sem';
        ERPStd{conds} = rmfield(ft_timelockgrandaverage(facfg,Subjects{:,conds}),'cfg');
        

    end
    
    save([Deci.Folder.Plot filesep 'GrandAverageSubFreq'],'ERPData','ERPStd');
    

else
    ERPData = Subjects;
    TotalCount = TrialCount;
    Deci.Plot.Freq.Wires.errorbars = 0;
end
clear Subjects;

%% Plot


if Deci.Plot.GA
    Deci.SubjectList = {'Group 1'};
end

for cond = 1:length(Deci.Plot.Draw)
    for subj = 1:size(ERPData,1)
        
        if Deci.Plot.ERP.Topo
            topo(subj)  = figure;
        end
        
        if Deci.Plot.ERP.Wires
            wire(subj)  = figure;
        end
        
        for subcond = 1:length(Deci.Plot.Draw{cond})
            
            if Deci.Plot.ERP.Topo
                if length(ERP.Channels) ~= 1
                    
                    
                    set(0, 'CurrentFigure', topo(subj))
                    topo(subj).Visible = 'on';
                    cirky(subj,subcond)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                    ft_topoplotER(cfg, ERPData{subj,Deci.Plot.Draw{cond}(subcond)});
                    
                    title([Deci.SubjectList{subj} ' '  Deci.Plot.Subtitle{cond}{subcond} ' trial count: ' num2str(TotalCount{subj,Deci.Plot.Draw{cond}(subcond)})]);
                    colorbar('vert');
                    map = colormap('jet'); %'hot' 'gray'
                    colormap(map);
                    
                else
                    close topo
                end
            end
            
            if  Deci.Plot.ERP.Wires
                
                set(0, 'CurrentFigure', wire(subj) )
                wire(subj).Visible = 'on';
                
                if ~Deci.Plot.ERP.errorbars
                    
                    h = plot(ERPData{subj,Deci.Plot.Draw{cond}(subcond)}.time,squeeze(mean(ERPData{subj,Deci.Plot.Draw{cond}(subcond)}.avg,1)));
                else
                    top = squeeze(mean(ERPData{subj,Deci.Plot.Draw{cond}(subcond)}.avg,1)) + squeeze(mean(squeeze(ERPStd{subj,Deci.Plot.Draw{cond}(subcond)}.avg),1))';
                    bot = squeeze(mean(ERPData{subj,Deci.Plot.Draw{cond}(subcond)}.avg,1)) - squeeze(mean(squeeze(ERPStd{subj,Deci.Plot.Draw{cond}(subcond)}.avg),1))';
                    
                    pgon = polyshape([ERPData{subj,Deci.Plot.Draw{cond}(subcond)}.time fliplr(ERPData{subj,Deci.Plot.Draw{cond}(subcond)}.time)],[top fliplr(bot)]);
                    b = plot(pgon,'HandleVisibility','off');
                    hold on
                    b.EdgeAlpha = 0;
                    b.FaceAlpha = .15;
                    h = plot(ERPData{subj,Deci.Plot.Draw{cond}(subcond)}.time,squeeze(mean(ERPData{subj,Deci.Plot.Draw{cond}(subcond)}.avg,1)));
                    h.Color = b.FaceColor;
                    h.LineWidth = 1;
                    
                    
                end
                ylabel('Amplitude')
                hold on
                
                legend(h.Parent,Deci.Plot.Subtitle{cond})
                title(h.Parent,[Deci.SubjectList{subj} ' ' Deci.Plot.Title{cond} ' Wire'])
                
                if ~isempty(Deci.Folder.Plot)
                    saveas(wire(subj),[Deci.Folder.Plot filesep Deci.Plot.Title{cond} filesep Deci.SubjectList{subj} '_wire'],Deci.Plot.Save.Format);
                end
                
            end
            
            
        end
        
    end
    
    for subj = 1:size(ERPData,1)
       
        
        if length(ERP.Channels) ~= 1
            
            if Deci.Plot.ERP.Topo
                set(0, 'CurrentFigure', topo(subj) )
                topo(subj).Visible = 'on';
                suptitle(Deci.Plot.Title{cond});
                
                for r = 1:length(cirky(:))
                    
                        cirky(r).CLim = [-1*max(abs([cirky.CLim])) max(abs([cirky.CLim]))];

                    
                end
                
                
                
                if ~isempty(Deci.Folder.Plot)
                    mkdir([Deci.Folder.Plot filesep Deci.Plot.Title{cond} filesep Deci.SubjectList{subj} '_topo']);
                    saveas(topo(subj),[Deci.Folder.Plot filesep Deci.Plot.Title{cond} filesep Deci.SubjectList{subj} '_topo'],Deci.Plot.Save.Format);
                end
                
            end
            
            
        end
        
        
        
    end
    
end


end