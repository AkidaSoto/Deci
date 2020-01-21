%Prob_Learning_Behavioral_Analysis
%
% by CTGill
% last edited 1/17/20
%
% Prob_Learning_Behavioral_Analysis.m assess the percentage of optimal
% responses made by each subject per block per condition and determines whether 
% each subject has achieved threshold performance in each stimulus
% condition. Threshold performance is established for each stimulus condition 
% independently by averaging performance across the final 3 blocks.
% Threshold performance level is set as 60% optimal responses. 
%
% Bar plots are generated and saved for each subject showing performance
% across all six blocks. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_figs = false;
save_figs = false;
save_sbjPerformance_data = true;

pt = 'C:\Users\CTGill\Documents\Probablistic_Learning\ProcessedData\Definition';
dirlist = dir('C:\Users\CTGill\Documents\Probablistic_Learning\ProcessedData\Definition\*');
Subject_List = {dirlist(~[dirlist.isdir]).name}';


save_fig_pt = 'C:\Users\CTGill\Documents\Probablistic_Learning\ProcessedData\Plot\Behv_performance\percent_corr_bar_graphs';
save_data_pt = 'C:\Users\CTGill\Documents\Probablistic_Learning\ProcessedData\Behavioral_Performance';

for sbj = 1 :size(Subject_List,1)
    load([pt filesep Subject_List{sbj}])
    
    cond_info = cfg.event(1:360,1:6);
    cond_info(:,7) = cfg.event(1:360,15);
    
    Blk_ends = [60,120,180,240,300,360];
    
    for i = 1 : 6
        if i == 1
            Blk{i} = cond_info(1:Blk_ends(i),:);
        else
            Blk{i} = cond_info(Blk_ends(i-1)+1:Blk_ends(i),:);
        end
    end
    
    ab_opt_total = 0;
    cd_opt_total = 0;
    ef_opt_total = 0;
    
    %these are for assessing whether subjects achieved threshold
    %performance in each of the three stimulus condititions. 
    %threshold performance is based on average performance over the last 3
    %blocks
    ab_thresh_count = 0;
    cd_thresh_count = 0;
    ef_thresh_count = 0;
    
    ab_thresh_corr = 0;
    cd_thresh_corr = 0;
    ef_thresh_corr = 0;
    
    for blk = 1:6
        ab_count = 0;
        cd_count = 0;
        ef_count = 0;
        ab_opt = 0;
        cd_opt = 0;
        ef_opt = 0;
        
        for trl = 1:60
            if Blk{blk}(trl,1) == 14
                ab_count = ab_count + 1;
                if Blk{blk}(trl,7) == 300
                    ab_opt = ab_opt + 1;
                end
            elseif Blk{blk}(trl,1) == 15
                cd_count = cd_count + 1;
                if Blk{blk}(trl,7) == 300
                    cd_opt = cd_opt + 1;
                end
            else
                ef_count = ef_count + 1;
                if Blk{blk}(trl,7) == 300
                    ef_opt = ef_opt + 1;
                end
            end
        end
        Stats{blk}.AB_percent_corr = ab_opt/ab_count;
        Stats{blk}.CD_percent_corr = cd_opt/cd_count;
        Stats{blk}.EF_percent_corr = ef_opt/ef_count;
        
        ab_opt_total = ab_opt_total + ab_opt;
        cd_opt_total = cd_opt_total + cd_opt;
        ef_opt_total = ef_opt_total + ef_opt;
        
        if blk > 3
            ab_thresh_count = ab_thresh_count + ab_count;
            cd_thresh_count = cd_thresh_count + cd_count;
            ef_thresh_count = ef_thresh_count + ef_count;
            
            ab_thresh_corr = ab_thresh_corr + ab_opt;
            cd_thresh_corr = cd_thresh_corr + cd_opt;
            ef_thresh_corr = ef_thresh_corr + ef_opt;
        end
        
    end
    Sbj_data{sbj}=Stats;
    
    ab_avg_corr(sbj) = ab_opt_total/120;
    cd_avg_corr(sbj) = cd_opt_total/120;
    ef_avg_corr(sbj) = ef_opt_total/120;
    
    performance(sbj,1) = ab_thresh_corr/ab_thresh_count;
    performance(sbj,2) = cd_thresh_corr/cd_thresh_count;
    performance(sbj,3) = ef_thresh_corr/ef_thresh_count;
end

%determine whether threshold was met for each sbj- 0 = not met 1 = met
for sbj = 1:size(performance,1)
    if performance(sbj,1) >= .6
        thresh_performance(sbj,1) = 1;
    else
        thresh_performance(sbj,1) = 0;
    end
    
    if performance(sbj,2) >= .6
        thresh_performance(sbj,2) = 1;
    else
        thresh_performance(sbj,2) = 0;
    end
    
    if performance(sbj,3) >= .6
        thresh_performance(sbj,3) = 1;
    else
        thresh_performance(sbj,3) = 0;
    end

    thresh_performance(sbj,4) = thresh_performance(sbj,1)+thresh_performance(sbj,2)+thresh_performance(sbj,3); 
end

if save_sbjPerformance_data
    for sbj = 1:size(performance,1)
        Thresh_performance.AB = thresh_performance(sbj,1);
        Thresh_performance.CD = thresh_performance(sbj,2);
        Thresh_performance.EF = thresh_performance(sbj,3);
        
        Thresh_performance.Group = thresh_performance(sbj,4);

        
        save([save_data_pt filesep 'Behv_Thresh_Performance_' Subject_List{sbj}],'Thresh_performance'); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting

if plot_figs == 1
    
    for i = 1:size(Subject_List,1)
        X=[Sbj_data{1,i}{1,1}.AB_percent_corr,Sbj_data{1,i}{1,1}.CD_percent_corr,Sbj_data{1,i}{1,1}.EF_percent_corr;...
            Sbj_data{1,i}{1,2}.AB_percent_corr,Sbj_data{1,i}{1,2}.CD_percent_corr,Sbj_data{1,i}{1,2}.EF_percent_corr;...
            Sbj_data{1,i}{1,3}.AB_percent_corr,Sbj_data{1,i}{1,3}.CD_percent_corr,Sbj_data{1,i}{1,3}.EF_percent_corr;...
            Sbj_data{1,i}{1,4}.AB_percent_corr,Sbj_data{1,i}{1,4}.CD_percent_corr,Sbj_data{1,i}{1,4}.EF_percent_corr;...
            Sbj_data{1,i}{1,5}.AB_percent_corr,Sbj_data{1,i}{1,5}.CD_percent_corr,Sbj_data{1,i}{1,5}.EF_percent_corr;...
            Sbj_data{1,i}{1,6}.AB_percent_corr,Sbj_data{1,i}{1,6}.CD_percent_corr,Sbj_data{1,i}{1,6}.EF_percent_corr;...
            performance(i,1),performance(i,2),performance(i,3)];
        
        fig = figure;
        bar(X)
        ylabel('Percent Correct');
        xlabel('Block numer');
        ylim([0 1])
        title([Subject_List{i}(size(Subject_List{i},2)-8:end-4) ' - Behv Performance']);
        xlim=get(gca,'xlim');
        xticks([1 2 3 4 5 6 7])
        xticklabels({'1','2','3','4','5','6','Thresh'})
        hold on;
        plot(xlim,[.6 .6], 'r--','LineWidth',2)
        a = annotation('textbox', [0.9, 0.625, 0, 0], 'string', 'Thresh');
        a.Color = 'red';
        
        prompt = 'Press Enter to continue ';
        x = input(prompt);
        
        
        if save_figs == true
            savefig(fig,[save_fig_pt filesep sprintf('%s_behv_performance.fig',Subject_List{i})]);
        end
        close all
    end
end

