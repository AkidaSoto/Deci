function Plottor_ERP(Deci)


cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';

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
            
            load([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.CondTitle{Conditions} filesep Chois{Channel} '.mat'],'raw');
            
            
            
            Chans{Channel} = raw;
            Chans{Channel}.label = Chois(Channel);
            
            if isfield(Deci.Plot.ERP,'PP')
                ppcfg = Deci.Plot.ERP.PP;
                
                Chans{Channel} = ft_preprocessing(ppcfg,  Chans{Channel});
                
            end
            
            clear raw;
        end
        
        if isfield(Chans{Channel},'trllength')
            trllen(subject_list,Conditions) = Chans{Channel}.trllength;
        else
            trllen(subject_list,Conditions) = nan;
        end
        
        if isfield(Chans{Channel},'lockers')
            LockNum = Deci.Analysis.Locks(ismember(Deci.Analysis.LocksTitle,Deci.Plot.Lock));
            lockers(subject_list,Conditions,:) = Chans{Channel}.lockers - Chans{Channel}.lockers(LockNum);
        else
            lockers(subject_list,Conditions,:) = nan;
        end
        
        
        acfg.parameter = 'avg';
        acfg.appenddim = 'chan';
        Subjects{subject_list,Conditions} = rmfield(ft_appendtimelock(acfg,Chans{:}),'cfg');
        Subjects{subject_list,Conditions}.dimord = 'chan_time';
        
        
    end
    clear Chans;
end

%% Bsl correct


for Conditions = 1:size(Subjects,2)
    for subject_list = 1:size(Subjects,1)
        
        if ~strcmpi(Deci.Plot.BslRef,Deci.Plot.Lock) || ~isempty(Deci.Plot.LockCond)
            
            if ~isempty(Deci.Plot.LockCond)
                BslCond =    Deci.Plot.CondTitle{Deci.Plot.LockCond(Conditions)};
            else
                BslCond =    Deci.Plot.CondTitle{Conditions};
            end
            
            for Channel = 1:length(Chois)
                
                load([Deci.Folder.Analysis filesep 'Volt_ERP' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.BslRef filesep BslCond filesep Chois{Channel} '.mat'],'raw');
                
                BslChan{Channel} = raw;
                BslChan{Channel}.label = Chois(Channel);
            end
            
            %acfg.parameter = 'avg';
            acfg.appenddim = 'chan';
            Bsl{subject_list,Conditions} = rmfield(ft_appendtimelock(acfg,BslChan{:}),'cfg');
            Bsl{subject_list,Conditions}.dimord = 'chan_time';
        else
            Bsl{subject_list,Conditions} =Subjects{subject_list,Conditions};
        end
        
        toi = round(Bsl{subject_list,Conditions}.time,4) >= Tois(1) & round(Bsl{subject_list,Conditions}.time,4) <= Tois(2);
        toi2 = Bsl{subject_list,Conditions}.time >= round(Deci.Plot.Bsl(1),4) & Bsl{subject_list,Conditions}.time <= round(Deci.Plot.Bsl(2),4);
        
        Bsl{subject_list,Conditions}.trial = nanmean(Bsl{subject_list,Conditions}.trial(:,:,toi2),3);
        Bsl{subject_list,Conditions}.trial = repmat(Bsl{subject_list,Conditions}.trial,[1 1 size(Subjects{subject_list,Conditions}.trial,3)]);
        
        Subjects{subject_list,Conditions}.trial =  Subjects{subject_list,Conditions}.trial - Bsl{subject_list,Conditions}.trial;
        Subjects{subject_list,Conditions}.time = Subjects{subject_list,Conditions}.time(toi);
        Subjects{subject_list,Conditions}.trial = Subjects{subject_list,Conditions}.trial(:,:,toi);
        
        Subjects{subject_list,Conditions}.avg = permute(mean(Subjects{subject_list,Conditions}.trial,1),[2 3 1]);
    end
end


%% Math
if ~isempty(Deci.Plot.Math)
    
    display(' ')
    display(['Doing ' num2str(length(Deci.Plot.Math)) ' Maths'] )
    
    for conds = 1:length(Deci.Plot.Math)
        for subj = 1:size(Subjects,1)
            scfg.parameter = 'avg';
            scfg.operation = Deci.Plot.Math{conds};
            evalc('MathData{subj} = ft_math(scfg,Subjects{subj,:})');
        end
        
        Subjects(:,size(Subjects,2)+1) = MathData;
    end
end

Subs = Deci.SubjectList;

if Deci.Plot.GrandAverage
    Deci.SubjectList = {'Group Average'};
end

%% Hemifield
if Deci.Plot.Hemiflip.do
    display(' ')
    display(['Applying Hemifield Flipping'] )
    
    Deci.Plot.Hemiflip = Exist(Deci.Plot.Hemiflip,'Type','Subtraction');
    
    if Deci.Plot.Topo.do
        error('cannot plot hemifield with topo')
    end
    
    
    for conds = 1:size(Subjects,2)
        for subj = 1:size(Subjects,1)
            
            hcfg.parameter = 'avg';
            hcfg.operation = 'x2 - x1';
            
            ContraCfg.channel = dc_getchans('even');
            ContraData{subj,conds} = ft_selectdata(ContraCfg,Subjects{subj,conds});
            ContraData{subj,conds} = hemifieldflip(ContraData{subj,conds});
            
            IpsiCfg.channel = dc_getchans('odd');
            IpsiData{subj,conds} = ft_selectdata(IpsiCfg,Subjects{subj,conds});
            
            if strcmpi(Deci.Plot.Hemiflip.Type,'Subtraction')
                Subjects{subj,conds} = ft_math(hcfg,IpsiData{subj,conds},ContraData{subj,conds});
            end
        end
        
    end
    
    if strcmpi(Deci.Plot.Hemiflip.Type,'Both')
        Subjects = cat(2,IpsiData,ContraData);
        
        %Deci.SubjectList = cat(2,cellfun(@(c) [c ' Ipsilateral'],Deci.SubjectList,'un',0),cellfun(@(c) [c ' Contralateral'],Deci.SubjectList,'un',0));
        
        drawlength =  length(Deci.Plot.Draw) + length(Deci.Plot.Math);
        
        Deci.Plot.Draw =  cellfun(@(c) [c arrayfun(@(d) d+drawlength,c,'un',1)],Deci.Plot.Draw,'un',0);
        Deci.Plot.Subtitle =  cellfun(@(c) [cellfun(@(d) [d ' Ipsilateral'],c,'un',0) cellfun(@(d) [d ' Contralateral'],c,'un',0)],Deci.Plot.Subtitle,'un',0);
        
    else
        Deci.Plot.Title = cellfun(@(c) [c ' Contra - Ipsilateral'],Deci.Plot.Title,'un',0);
    end
end
%% Data Management
if size(Subjects,1) == 1
    Deci.Plot.GrandAverage = false;
end

if Deci.Plot.GrandAverage
    if any(~isnan(trllen))
        trlstd = nanstd(trllen,[],1);
        trllen = nanmean(trllen,1);
    end
    
    if any(~isnan(lockers))
        lockersstd = nanstd(lockers,[],1);
        lockers = nanmean(lockers,1);
    end
end

for conds = 1:size(Subjects,2)
    
    if Deci.Plot.GrandAverage
        facfg.parameter =  'avg';
        facfg.type = 'mean';
        facfg.keepindividual = 'yes';
        evalc('ErpData{1,conds} = ft_timelockgrandaverage(facfg,Subjects{:,conds});');
        ErpData{1,conds}.avg = ErpData{1,conds}.individual;
        ErpData{1,conds} = rmfield(ErpData{1,conds},'cfg');
        
    else
        Deci.Plot.Stat.do = false;
        ErpData(:,conds) = Subjects(:,conds);
    end
    
    for subj = 1:size(ErpData(:,conds),1)
        if Deci.Plot.Topo.do
            tcfg = [];
            tcfg.nanmean = Deci.Plot.nanmean;
            
            tcfg.latency = Deci.Plot.Topo.Toi;
            tcfg.channel = Deci.Plot.Topo.Channel;
            
            topodata{subj,conds} = ft_selectdata(tcfg,ErpData{subj,conds});
            
            tcfg.avgovertime = 'yes';
            topotdata{subj,conds} = ft_selectdata(tcfg,topodata{subj,conds});
        end
        
        if Deci.Plot.Wire.do
            tcfg = [];
            tcfg.nanmean = Deci.Plot.nanmean;
            
            tcfg.latency = Deci.Plot.Wire.Toi;
            tcfg.channel = Deci.Plot.Wire.Channel;
            
            wiredata{subj,conds} = ft_selectdata(tcfg,ErpData{subj,conds});
            
            tcfg.avgoverchan = 'yes';
            wiretdata{subj,conds} = ft_selectdata(tcfg,wiredata{subj,conds});
            
        end
        
        if Deci.Plot.Bar.do
            tcfg = [];
            tcfg.nanmean = Deci.Plot.nanmean;
            
            tcfg.latency = Deci.Plot.Bar.Toi;
            tcfg.channel = Deci.Plot.Bar.Channel;
            
            bardata{subj,conds} = ft_selectdata(tcfg,ErpData{subj,conds});
            
            
            tcfg.avgoverchan = 'yes';
            tcfg.avgovertime = 'yes';
            
            bartdata{subj,conds} = ft_selectdata(tcfg,bardata{subj,conds});
            
        end
    end
end

if Deci.Plot.Bar.do

    bardatas = num2cell(cell2mat(cellfun(@(c) mean(c.avg,3),bardata,'UniformOutput',false)));
    save([Deci.Folder.Plot filesep  Deci.Plot.Title{1} ' Bar Data.mat' ],'bardatas');
end

if Deci.Plot.Stat.do
    
    display(' ')
    display('Calculating Statistics')
    
    neighbours       = load('easycap_neighbours','neighbours');
    Deci.Plot.Stat.neighbours = neighbours.neighbours;
    Deci.Plot.Stat.ivar = 1;
    
    for conds = 1:length(Deci.Plot.Draw)
        design = [];
        
        switch Deci.Plot.Stat.Comp
            case 'Conditions'
                if length(Deci.Plot.Draw{conds}) ~= 1
                    Deci.Plot.Stat.uvar = 2;
                    
                    for subcond = 1:length(Deci.Plot.Draw{conds})
                        for subj = 1:size(Subjects,1)
                            design(1,subj+size(Subjects,1)*[subcond-1]) =  subcond;
                            design(2,subj+size(Subjects,1)*[subcond-1]) = subj;
                        end
                    end
                    
                    Deci.Plot.Stat.design = design;
                    
                    if length(Deci.Plot.Draw{conds}) > 2
                        Deci.Plot.Stat.tail = 1;
                        Deci.Plot.Stat.statistic = 'depsamplesFmultivariate';
                        Deci.Plot.Stat.clustertail      = 1;
                    else
                        Deci.Plot.Stat.statistic = 'depsamplesT';
                        Deci.Plot.Stat.tail = 0;
                        Deci.Plot.Stat.clustertail      = 0;
                    end
                    
                    if Deci.Plot.Topo.do
                        [topostat{conds}] = ft_timelockstatistics(Deci.Plot.Stat, topotdata{Deci.Plot.Draw{conds}});
                    end
                    
                    if Deci.Plot.Wire.do
                        [wirestat{conds}] = ft_timelockstatistics(Deci.Plot.Stat, wiretdata{:,Deci.Plot.Draw{conds}});
                    end
                    
                    if Deci.Plot.Bar.do
                        [barstat{conds}] = ft_timelockstatistics(Deci.Plot.Stat, bartdata{Deci.Plot.Draw{conds}});
                    end
                else
                    
                    if Deci.Plot.Topo.do
                        topostat{conds}.mask = permute(ones(size(topotdata{Deci.Plot.Draw{conds}}.avg(1,:,:,:))),[2 3 4 1]);
                    end
                    
                    if Deci.Plot.Wire.do
                        [wirestat{conds}.mask] = permute(ones(size(wiretdata{:,Deci.Plot.Draw{conds}}.avg(1,:,:,:))),[2 3 4 1])*0;
                    end
                    
                    if Deci.Plot.Bar.do
                        [barstat{conds}.mask] = permute(ones(size(bartdata{Deci.Plot.Draw{conds}}.avg(1,:,:,:))),[2 3 4 1]);
                    end
                    
                end
            case 'Bsl'
                Deci.Plot.Stat.tail = 0;
                Deci.Plot.Stat.statistic = 'indepsamplesT';
                
                Deci.Plot.Stat.design = ones(size(Subjects,1));
                
                if Deci.Plot.Topo.do
                    alltopotdata = ft_timelockgrandaverage([],topotdata{Deci.Plot.Draw{conds}});
                    alltopotdata.dimord = 'rpt_chan_time';
                    [topostat{conds}] = ft_timelockstatistics(Deci.Plot.Stat, alltopotdata);
                end
                
                
                if Deci.Plot.Wire.do
                    allwiretdata = ft_timelockgrandaverage([],wiretdata{Deci.Plot.Draw{conds}});
                    allwiretdata.dimord = 'rpt_chan_time';
                    [wirestat{conds}] = ft_timelockstatistics(Deci.Plot.Stat, allwiretdata);
                end
                
                if Deci.Plot.Bar.do
                    allbartdata = ft_timelockgrandaverage([],bartdata{Deci.Plot.Draw{conds}});
                    allbartdata.dimord = 'rpt_chan_time';
                    [barstat{conds}] = ft_timelockstatistics(Deci.Plot.Stat, allbartdata);
                end
            case 'Corr'
                
                for subcond = 1:length(Deci.Plot.Draw{conds})
                    for subj = 1:size(Subjects,1)
                        design(1,subj+size(Subjects,1)*[subcond-1]) =  subcond;
                    end
                end
                
                Deci.Plot.Stat.design = design;
                Deci.Plot.Stat.statistic = 'ft_statfun_correlationT';
                
                
                if Deci.Plot.Topo.do
                    [topostat{conds}] = ft_timelockstatistics(Deci.Plot.Stat, topotdata{Deci.Plot.Draw{conds}});
                end
                
                if Deci.Plot.Wire.do
                    [wirestat{conds}] = ft_timelockstatistics(Deci.Plot.Stat, wiretdata{:,Deci.Plot.Draw{conds}});
                end
                
                if Deci.Plot.Bar.do
                    [barstat{conds}] = ft_timelockstatistics(Deci.Plot.Stat, bartdata{Deci.Plot.Draw{conds}});
                end
        end
        
    end
end

%% Plot

for cond = 1:length(Deci.Plot.Draw)
    if Deci.Plot.Stat.do
        
        if Deci.Plot.Topo.do
            tacfg = cfg;
            tacfg.parameter = 'stat';
            topostat{cond}.mask = double(topostat{cond}.mask);
            topostat{cond}.mask(topostat{cond}.mask == 0) = .2;
            tacfg.clim = 'maxmin';
            tacfg.maskparameter ='mask';
            tacfg.colormap = Deci.Plot.ColorMap;
            
            if Deci.Plot.Stat.FPlots
                topot(cond)  = figure;
                topot(cond).Visible = 'on';
                
                ft_topoplotTFR(tcfg, topostat{cond});
                title([Deci.Plot.Stat.Type ' ' Deci.Plot.Title{cond} ' (alpha = ' num2str(Deci.Plot.Stat.alpha) ')']);
                dc_pmask(topot(cond))
            end
        end
        
        
    end
    
    clear cirky
    for subj = 1:size(ErpData,1)
        if Deci.Plot.Topo.do
            topo(subj)  = figure;
            
            if Deci.Plot.Stat.do
                dc_pmask(topo)
            end
        end
        
        if Deci.Plot.Wire.do
            wire(subj)  = figure;
        end
        
        for subcond = 1:length(Deci.Plot.Draw{cond})
            if Deci.Plot.Topo.do
                set(0, 'CurrentFigure', topo(subj))
                topo(subj).Visible = 'on';
                cirky(subj,subcond)    =  subplot(length(Deci.Plot.Draw{cond}),1,subcond);
                
                pcfg = cfg;
                if Deci.Plot.Stat.do
                    pcfg.clim = 'maxmin';
                    pcfg.maskparameter ='mask';
                    topodata{subj,Deci.Plot.Draw{cond}(subcond)}.mask = repmat(topostat{cond}.mask',[size(topodata{subj,Deci.Plot.Draw{cond}(subcond)}.avg,1) 1 length(topodata{subj,Deci.Plot.Draw{cond}(subcond)}.freq) length(topodata{subj,Deci.Plot.Draw{cond}(subcond)}.time)]);
                end
                
                pcfg.imagetype = 'contourf';
                pcfg.comment = 'no';
                pcfg.style = 'fill';
                pcfg.markersymbol = '.';
                pcfg.colormap = Deci.Plot.ColorMap;
                pcfg.colorbar = 'yes';
                pcfg.contournum = 15;
                ft_topoplotER(pcfg, topodata{subj,Deci.Plot.Draw{cond}(subcond)});
                
                if Deci.Plot.Draw{cond}(subcond) <= size(trllen,2)
                    title([Deci.SubjectList{subj} ' '  Deci.Plot.Subtitle{cond}{subcond} ' (' num2str(trllen(subj,Deci.Plot.Draw{cond}(subcond))) ')']);
                    
                else
                    title([Deci.SubjectList{subj} ' '  Deci.Plot.Subtitle{cond}{subcond}]);
                end
                
                colorbar('vert');
                
            end
            
            if Deci.Plot.Wire.do
                set(0, 'CurrentFigure', wire(subj) )
                wire(subj).Visible = 'on';
                
                
                if Deci.Plot.GrandAverage
                    top = squeeze(nanmean(nanmean(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.avg,1),2)) + squeeze(nanstd(nanmean(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.avg,2),[],1))/sqrt(size(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.avg,1));
                    bot = squeeze(nanmean(nanmean(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.avg,1),2)) - squeeze(nanstd(nanmean(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.avg,2),[],1))/sqrt(size(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.avg,1));
                    
                    pgon = polyshape([wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.time fliplr(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.time)],[top' fliplr(bot')],'Simplify', false);
                    b(subcond) = plot(pgon,'HandleVisibility','off');
                    hold on
                    b(subcond).EdgeAlpha = 0;
                    b(subcond).FaceAlpha = .15;
                end
                
                if Deci.Plot.Stat.do
                    %pcfg.maskparameter ='mask'
                    
                    sigs  = double(wirestat{cond}.mask);
                    sigs(sigs == 0) = nan;
                    
                    colors = {'r' 'g' 'b'};
                    
                    mod = min(cell2mat(arrayfun(@(c) min(mean(c.avg,1),[],'all'),[wiredata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)));
                    modfunc = @minus;
                    if max(cell2mat(arrayfun(@(c) max(abs(mean(c.avg,1)),[],'all'),[wiredata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false))) < 0
                        mod = max(cell2mat(arrayfun(@(c) max(mean(c.avg,1),[],'all'),[wiredata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)));
                        modfunc = @plus;
                    end
                    diffy = diff(ylim);
                    hold on
                    for z = 1:size(sigs,1)
                        plot(wiredata{subj,Deci.Plot.Draw{cond}(subcond)}.time,sigs(z,:)*modfunc(mod,diffy*[.01*z]),'HandleVisibility','off','LineWidth',5,'Color',colors{z});
                    end
                end
            end
            
            
        end
        
        
        
        if Deci.Plot.Wire.do
            set(0, 'CurrentFigure', wire(subj) )
            wire(subj).Visible = 'on';
            
            pcfg = cfg;
            pcfg.clim = 'maxmin';
            
            
            %pcfg.ylim = ylim;
            pcfg.graphcolor = lines;
            pcfg.linewidth = 1;
            ft_singleplotER(pcfg,wiredata{subj,Deci.Plot.Draw{cond}});
            set(gca, 'YDir', 'reverse');
            
            if Deci.Plot.GrandAverage
                arrayfun(@(c) uistack(c,'top'),b);
                clear b
            end
            
            axis tight
            hold on
            plot([wiredata{cond}.time(1), wiredata{cond}.time(end)], [0 0], 'k--') % hor. line
            plot([0 0], ylim, 'k--') % vert. l
            
            if Deci.Plot.Stat.do
                boxes = wire(subj).Children(2).Children.findobj('Type','Patch');
                for bb = 1:length(boxes)
                    if ~isempty(boxes)
                        boxes(bb).YData(boxes(bb).YData == min(boxes(bb).YData)) = min(ylim);
                        boxes(bb).YData(boxes(bb).YData == max(boxes(bb).YData)) = max(ylim);
                        
                        boxes(bb).FaceAlpha = .35;
                        uistack(boxes(bb),'bottom')
                        boxes(bb).HandleVisibility = 'off';
                    end
                end
            end
            
            if max(Deci.Plot.Draw{cond}) <= size(trllen,2)
                legend(arrayfun(@(a,b) [a{1} ' (' num2str(b) ')'] ,Deci.Plot.Subtitle{cond},trllen(subj,Deci.Plot.Draw{cond}),'UniformOutput',false));
            else
                legend([ Deci.Plot.Subtitle{cond}]);
            end
            
            title([Deci.SubjectList{subj} ' ' Deci.Plot.Title{cond}], 'Interpreter', 'none');
            set(legend, 'Interpreter', 'none')
            xlim([wiredata{cond}.time(1) wiredata{cond}.time(end)])
            xlabel('Time');
            
            
            for subcond = 1:length(Deci.Plot.Draw{cond})
                if Deci.Plot.Draw{cond}(subcond) <= size(lockers,2)
                    xlims = xlim;
                    ylims = ylim;
                    
                    for locks = 1:length([lockers(subj,Deci.Plot.Draw{cond}(subcond),:)])
                        hold on
                        
                        locktime = [lockers(subj,Deci.Plot.Draw{cond}(subcond),locks)/1000];
                        
                        
                        if locktime > xlims(1) && locktime < xlims(2)
                            plotlock = plot([locktime locktime], ylims,'LineWidth',2,'Color','k','LineStyle','--','HandleVisibility','off');
                            plotlock.Color(4) = .2;
                            
                            if Deci.Plot.GrandAverage
                                lockstd = [lockersstd(subj,Deci.Plot.Draw{cond}(subcond),locks)/1000];
                                
                                
                                if [locktime - lockstd] < xlims(1)
                                    lockpstd(1) = xlims(1);
                                else
                                    lockpstd(1) = [locktime - lockstd];
                                end
                                
                                if [locktime + lockstd] > xlims(2)
                                    lockpstd(2) = xlims(2);
                                else
                                    lockpstd(2) = [locktime + lockstd];
                                end
                                
                                lockpgon = polyshape([lockpstd fliplr(lockpstd)],sort([ylims ylims]),'Simplify', false);
                                lockb = plot(lockpgon,'HandleVisibility','off');
                                hold on
                                lockb.EdgeAlpha = 0;
                                lockb.FaceAlpha = .15;
                                lockb.FaceColor = plotlock.Color;
                                
                                arrayfun(@(c) uistack(c,'top'),lockb);
                                clear lockb
                                
                            end
                        end
                        ylim(ylims)
                    end
                end
                
                if Deci.Plot.Wire.do
                    
                    lims = minmax(cell2mat(arrayfun(@(c) c.Children.findobj('Type','Axes').YLim,wire,'UniformOutput',false)));
                    for r = 1:length(wire(:))
                        set(0, 'CurrentFigure', wire(r))
                        wire(r).Children.findobj('Type','Axes').YLim = lims;
                    end
                end
                
            end
            
        end
        
        if Deci.Plot.Bar.do
            barfig(subj)  = figure;
            set(0, 'CurrentFigure', barfig(subj) )
            barfig(subj).Visible = 'on';
            
            CleanBars(mean(cell2mat(arrayfun(@(c) nanmean(nanmean(nanmean(c.avg,2),3),4),[bardata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)),1), ...
                nanstd(cell2mat(arrayfun(@(c) nanmean(nanmean(nanmean(c.avg,2),3),4),[bardata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)),[],1) ...
                /sqrt(size(cell2mat(arrayfun(@(c) nanmean(nanmean(nanmean(c.avg,2),3),4),[bardata{subj,Deci.Plot.Draw{cond}}],'UniformOutput',false)),1)));
            
            l = legend(Deci.Plot.Subtitle{cond});
            title([Deci.SubjectList{subj} ' ' Deci.Plot.Title{cond} ' Wire'])
            set(l, 'Interpreter', 'none')
            if Deci.Plot.Stat.do
                hold on
                
                if barstat{cond}.prob < Deci.Plot.Stat.alpha
                    plot([.75 1.25],[max(ylim) max(ylim)]*.90, '-k', 'LineWidth',2);
                    plot([1],[max(ylim)]*.95, '*k');
                end
            end
        end
    end
end

for subj = 1:size(ErpData,1)
    if length(Deci.Plot.Topo.Channel) ~= 1
        
        if Deci.Plot.Topo.do
            set(0, 'CurrentFigure', topo(subj) )
            suptitle(Deci.Plot.Title{cond});
            for r = 1:length(cirky(:))
                if length(Deci.Plot.Roi) == 2 && isnumeric(Deci.Plot.Roi)
                    cirky(r).CLim = Deci.Plot.Freq.Roi;
                elseif strcmp(Deci.Plot.Roi,'maxmin')
                    if ~ isempty(cirky(r).Children.UserData)
                        cirky(r).CLim = [min(arrayfun(@(c) min(c.Children.UserData(:)),cirky(:))) max(arrayfun(@(c) max(c.Children.UserData(:)),cirky(:)))];
                    else
                        cirky(r).CLim= [min(arrayfun(@(c) min(c.Children.ZData(:)),cirky(:))) max(arrayfun(@(c) max(c.Children.ZData(:)),cirky(:)))];
                    end
                    %cirky(r).Children.LevelList = [linspace(cirky(r).CLim(1),cirky(r).CLim(2),12)];
                elseif strcmp(Deci.Plot.Roi,'maxabs')
                    if ~isempty(cirky(r).Children.findobj('Type','Contour').UserData)
                        cirky(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.findobj('Type','Contour').UserData(:))),cirky(:))) max(arrayfun(@(c) max(abs(c.Children.findobj('Type','Contour').UserData(:))),cirky(:)))];
                    else
                        cirky(r).CLim = [-1*max(arrayfun(@(c) max(abs(c.Children.findobj('Type','Contour').ZData(:))),cirky(:))) max(arrayfun(@(c) max(abs(c.Children.findobj('Type','Contour').ZData(:))),cirky(:)))];
                    end
                    %cirky(r).Children.findobj('Type','Contour').LevelList = [linspace(cirky(r).CLim(1),cirky(r).CLim(2),12)];
                end
            end
            
            if ~isempty(Deci.Folder.Plot)
                %mkdir([Deci.Folder.Plot filesep Deci.Plot.Title{cond}]);
                %saveas(topo(subj),[Deci.Folder.Plot filesep Deci.Plot.Title{cond} filesep Deci.SubjectList{subj} '_topo'],Deci.Plot.Save.Format);
            end
            
        end
        
        
        
    end
    
end
end

