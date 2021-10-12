function Plottor_PhaseDiffs(Deci,Params)

%% Load
cfg        = [];
cfg.layout = Deci.Layout.eye;
cfg.channel = 'all';
cfg.interactive = 'yes';

disp('----------------------');
display(' ')

display(['Plotting PhaseDiffs']);
warning off
for ConnList = 1:length(Params.List)
    
    
    Current = Params.List{ConnList};
    chanl = Current{1};
    chanh = Current{2};
    freqlow = Current{3};
    freqhigh = Current{4};
    
    %% do plots for 1 conoi at a time

        clear Subjects
        if ischar(chanl)
            chanl = {chanl};
        end
        
        if ischar(chanh)
            chanh = {chanh};
        end
        
        if isequal(chanl,{'Reinhart-All'})
            chanl = dc_getchans('noeyes');
        end
        
        if isequal(chanh,{'Reinhart-All'})
            chanh = dc_getchans('noeyes');
        end
        
        if ischar(freqlow)
            freqlow = {freqlow};
        end
        
        if ischar(freqhigh)
            freqhigh = {freqhigh};
        end
 
        if isequal(chanl,chanh)
            chancmb = [chanl; chanh]';
        else
            chancmb = [chanl chanh];
            
            if size(chancmb,2) ~= 1
                chancmb = chancmb(combvec(1:length(chanl),[1:length(chanh)]+length(chanl)))';
            end
        end

       chancmb = chancmb(cellfun(@(a,b) ~isequal(a,b),chancmb(:,1),chancmb(:,2)),:);
  

        if isequal(freqlow,freqhigh)
            freqcmb = [freqlow; freqhigh]';
        else
            freqcmb = [freqlow freqhigh];
            freqcmb = freqcmb(combvec(1:length(freqlow),[1:length(freqhigh)]+length(freqlow)))';
        end
        
        if isequal(size(freqcmb), [2 1])
            freqcmb = freqcmb';
        end
        
        
        %% load 1 full conoi
        clear sub_cond Foi tempdata
        

        for Conditions = 1:length(Deci.Plot.CondTitle)
            
            for  subject_list = 1:length(Deci.SubjectList)
                
                display(['Loading Plottor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
                
                for choicmb = 1:size(chancmb,1)
                    
                    for foicmb = 1:size(freqcmb,1)
                        
                        connfile = strjoin([chancmb(choicmb,1) chancmb(choicmb,2) freqcmb(foicmb,1) freqcmb(foicmb,2) 'phasediff'],'_');
                        
                        
                        load([Deci.Folder.Analysis filesep 'Extra' filesep 'Conn' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock filesep Deci.Plot.CondTitle{Conditions} filesep connfile  '.mat'],'phase_diff');
                        
                        sub_freq{subject_list,Conditions,choicmb,foicmb} = phase_diff;


                         clear phase_diff
                    end
                end
            end
            
            if Deci.Plot.Extra.Phase_diff.toivbsl.do
                
                bsl = cellfun(@(c) c.toi(:,2),sub_freq(:,Conditions,:,:),'UniformOutput',false);
                shift = cellfun(@mean,bsl,'un',0);
                bsl = cellfun(@(a,b) a-b,bsl,shift,'un',0);
                bsl = cat(1,bsl{:});
                bsl = bsl(:);
                
                toi = cellfun(@(c) c.toi(:,1),sub_freq(:,Conditions,:,:),'UniformOutput',false);
                toi = cellfun(@(a,b) a-b,toi,shift,'un',0);
                toi = cat(1,toi{:});
                toi = toi(:);
                
                [wwtest] = circ_wwtest(toi,bsl);
                a = figure;
                a.Visible = 'on';
                polarhistogram(toi);
                hold on
                polarhistogram(bsl);
                
                title([connfile ' toi v bsl ' Deci.Plot.CondTitle{Conditions} ' p value ' num2str(wwtest)],'Interpreter','none')
                legend({'toi','bsl'})
            end
            
        end
        
%% plot

for pair = 1:length(Deci.Plot.Extra.Phase_diff.pairings)
    if length(Deci.Plot.Extra.Phase_diff.pairings{pair}) ~= 2
       warning('PhaseDiffs can only be done with pairings') 
       continue
    end
    
    if iscell(Deci.Plot.Extra.Phase_diff.pairings{pair}(1))
        toi1 = cellfun(@(c) c.toi(:,1),sub_freq(:,Deci.Plot.Extra.Phase_diff.pairings{pair}{1},:,:),'UniformOutput',false);
        shift = cellfun(@mean,toi1,'un',0);
        toi1 = cellfun(@(a,b) a-b,toi1,shift,'un',0);
        toi1 = cat(1,toi1{:});
        toi1 = toi1(:);
    end
    
    if iscell(Deci.Plot.Extra.Phase_diff.pairings{pair}(2))
        toi2 = cellfun(@(c) c.toi(:,1),sub_freq(:,Deci.Plot.Extra.Phase_diff.pairings{pair}{2},:,:),'UniformOutput',false);
        toi2 = cellfun(@(a,b) a-b,toi2,shift,'un',0);
        toi2 = cat(1,toi2{:});
        toi2 = toi2(:);
    end
    
    [wwtest] = circ_wwtest(toi1,toi2);
    a = figure;
    a.Visible = 'on';
    polarhistogram(toi1);
    hold on
    polarhistogram(toi2);
    

    title([Deci.Plot.Extra.Phase_diff.title{pair} ' p value ' num2str(wwtest)],'Interpreter','none')
    legend(Deci.Plot.Extra.Phase_diff.subtitle{pair})
    
end

clear sub_freq
end