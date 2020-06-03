%Prob_Learning_Behavioral_Analysis
%
% by CTGill
% last edited 3/10/20
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
%% settings
clear;
clc;

plot_figs = false;
save_figs = false;
save_sbjPerformance_data = true;
display_group_stats = true;
plot_RTs = true;
save_Sbj_RTs = false;
save_group_RTs = false;

pt = 'C:\Users\ctgil\Documents\Probablistic_Learning\ProcessedData\Definition';
dirlist = dir('C:\Users\ctgil\Documents\Probablistic_Learning\RawData_new\*.eeg');
Subject_List = {dirlist(~[dirlist.isdir]).name}';

save_fig_pt = 'C:\Users\ctgil\Documents\Probablistic_Learning\ProcessedData\Plot\Behv_performance';
save_data_pt = 'C:\Users\ctgil\Documents\Probablistic_Learning\ProcessedData\Behavioral_Performance';

%% loop through all subjects
for sbj = 1 :size(Subject_List,1)
    load([pt filesep Subject_List{sbj}(1:end-4) '.mat'])%load sbj data
    
    %keep relevant markers
    cond_info = cfg.event(:,1:6);
    cond_info(:,7) = cfg.event(:,15);
    
    %determine response time
    rt = (cfg.trl(:,4)-cfg.trl(:,5))/1000;
    
    %define block ends
    Blk_ends = [60,120,180,240,300,360];
    
    %split trial data up into blocks
    for i = 1 : 6
        if i == 1
            Blk{i} = cond_info(1:Blk_ends(i),:);
        else
            Blk{i} = cond_info(Blk_ends(i-1)+1:Blk_ends(i),:);
        end
    end
    
    %track total number of optimal responses per stim type for each sbj
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
    
    %track trial numbers for all non-responses(i.e. which trials had no resp)
    non_resp_trials = [find(isnan(cond_info(:,7)))];
    
    %loop through each block
    for blk = 1:6
        ab_blk_count = 0;
        cd_blk_count = 0;
        ef_blk_count = 0;
        ab_blk_opt = 0;
        cd_blk_opt = 0;
        ef_blk_opt = 0;
        
        
        %loop through trials and count the number of each stim condition
        %and the number of optimal responses per stim condition
        for trl = 1:60
            if Blk{blk}(trl,1) == 14
                ab_blk_count = ab_blk_count + 1;
                if Blk{blk}(trl,7) == 300
                    ab_blk_opt = ab_blk_opt + 1;
                end
                
            elseif Blk{blk}(trl,1) == 15
                cd_blk_count = cd_blk_count + 1;
                if Blk{blk}(trl,7) == 300
                    cd_blk_opt = cd_blk_opt + 1;
                end
            else
                ef_blk_count = ef_blk_count + 1;
                if Blk{blk}(trl,7) == 300
                    ef_blk_opt = ef_blk_opt + 1;
                end
            end
        end
        
        %calculate percent optimal responses per stim condition per block
        Blk_perf{blk}.AB_percent_corr = ab_blk_opt/ab_blk_count;
        Blk_perf{blk}.CD_percent_corr = cd_blk_opt/cd_blk_count;
        Blk_perf{blk}.EF_percent_corr = ef_blk_opt/ef_blk_count;
        
        %calculate the total number of optimal responses across all blocks
        ab_opt_total = ab_opt_total + ab_blk_opt;
        cd_opt_total = cd_opt_total + cd_blk_opt;
        ef_opt_total = ef_opt_total + ef_blk_opt;
        
        %calculate the number of optimal responses across the last 3 blocks
        %for each stim condition in order to determine if sbj achieves
        %threshold performance
        if blk > 3
            ab_thresh_count = ab_thresh_count + ab_blk_count;
            cd_thresh_count = cd_thresh_count + cd_blk_count;
            ef_thresh_count = ef_thresh_count + ef_blk_count;
            
            ab_thresh_corr = ab_thresh_corr + ab_blk_opt;
            cd_thresh_corr = cd_thresh_corr + cd_blk_opt;
            ef_thresh_corr = ef_thresh_corr + ef_blk_opt;
        end
    end
    
    %create data struct for saving behavioral performance
    Sbj_Behv_data.Sbj = Subject_List{sbj}(1:end-4);
    Sbj_Behv_data.Blk_performance = Blk_perf;
    Sbj_Behv_data.All_Blk_avg.ab_avg_opt_resp = ab_opt_total/120;
    Sbj_Behv_data.All_Blk_avg.cd_avg_opt_resp = cd_opt_total/120;
    Sbj_Behv_data.All_Blk_avg.ef_avg_opt_resp = ef_opt_total/120;
    Sbj_Behv_data.thresh_performance.AB = ab_thresh_corr/ab_thresh_count;
    Sbj_Behv_data.thresh_performance.CD = cd_thresh_corr/cd_thresh_count;
    Sbj_Behv_data.thresh_performance.EF = ef_thresh_corr/ef_thresh_count;
    Sbj_Behv_data.RT = rt;
    Sbj_Behv_data.non_resp_trials = [find(isnan(Sbj_Behv_data.RT))];
    Sbj_Behv_data.early_resp_trials = [find(Sbj_Behv_data.RT < .1)]; %trials in which RT was less than 100ms
    

    
    %Determine whether threshold was met for each stim condition
    if ((ab_thresh_corr/ab_thresh_count) >= .6)
        Sbj_Behv_data.AB_thresh_met = 1;
    else
        Sbj_Behv_data.AB_thresh_met = 0;
    end
    
    if ((cd_thresh_corr/cd_thresh_count) >= .6)
        Sbj_Behv_data.CD_thresh_met = 1;
    else
        Sbj_Behv_data.CD_thresh_met = 0;
    end
    
    if ((ef_thresh_corr/ef_thresh_count) >= .6)
        Sbj_Behv_data.EF_thresh_met = 1;
    else
        Sbj_Behv_data.EF_thresh_met = 0;
    end
    
    
    %assign subject to behavioral performance group based on the number of
    %stim conditions in which thresh performance was met
    if ((Sbj_Behv_data.AB_thresh_met + Sbj_Behv_data.CD_thresh_met + Sbj_Behv_data.EF_thresh_met) == 0 || (Sbj_Behv_data.AB_thresh_met + Sbj_Behv_data.CD_thresh_met + Sbj_Behv_data.EF_thresh_met) == 1)
        Sbj_Behv_data.Behav_performance_group = 1;
    elseif (Sbj_Behv_data.AB_thresh_met + Sbj_Behv_data.CD_thresh_met + Sbj_Behv_data.EF_thresh_met) == 2
        Sbj_Behv_data.Behav_performance_group = 2;
    else
        Sbj_Behv_data.Behav_performance_group = 3;
    end
    
    
    %save data
    if save_sbjPerformance_data
        save([save_data_pt filesep 'Behv_Performance_data_' Subject_List{sbj}(1:end-4)],'Sbj_Behv_data');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plotting Performance
    
    if plot_figs == 1
        X=[Sbj_Behv_data.Blk_performance{1,1}.AB_percent_corr,Sbj_Behv_data.Blk_performance{1,1}.CD_percent_corr,Sbj_Behv_data.Blk_performance{1,1}.EF_percent_corr;...
            Sbj_Behv_data.Blk_performance{1,2}.AB_percent_corr,Sbj_Behv_data.Blk_performance{1,2}.CD_percent_corr,Sbj_Behv_data.Blk_performance{1,2}.EF_percent_corr;...
            Sbj_Behv_data.Blk_performance{1,3}.AB_percent_corr,Sbj_Behv_data.Blk_performance{1,3}.CD_percent_corr,Sbj_Behv_data.Blk_performance{1,3}.EF_percent_corr;...
            Sbj_Behv_data.Blk_performance{1,4}.AB_percent_corr,Sbj_Behv_data.Blk_performance{1,4}.CD_percent_corr,Sbj_Behv_data.Blk_performance{1,4}.EF_percent_corr;...
            Sbj_Behv_data.Blk_performance{1,5}.AB_percent_corr,Sbj_Behv_data.Blk_performance{1,5}.CD_percent_corr,Sbj_Behv_data.Blk_performance{1,5}.EF_percent_corr;...
            Sbj_Behv_data.Blk_performance{1,6}.AB_percent_corr,Sbj_Behv_data.Blk_performance{1,6}.CD_percent_corr,Sbj_Behv_data.Blk_performance{1,6}.EF_percent_corr;...
            Sbj_Behv_data.thresh_performance.AB,Sbj_Behv_data.thresh_performance.CD,Sbj_Behv_data.thresh_performance.EF];
        
        fig = figure;
        bar(X)
        ylabel('Percent Correct');
        xlabel('Block numer');
        ylim([0 1])
        title([Subject_List{sbj}(size(Subject_List{sbj},2)-8:end-4) ' - Behv Performance']);
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
            savefig(fig,[save_fig_pt filesep 'percent_corr_bar_graphs' filesep sprintf('%s_behv_performance.fig',Subject_List{sbj}(1:end-4))]);
        end
        close all
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% display group stats

pt2 = 'C:\Users\ctgil\Documents\Probablistic_Learning\ProcessedData\Behavioral_Performance';

group1_count = 0; %combination of subjects who reached thresh performance on either 0 or 1 stim conditions.
group2_count = 0; %number of sbjs who reached thresh performance on 2 of the stim conditions
group3_count = 0; %number of sbjs who reached thresh performance on 3 of the stim conditions

AB_count = 0; %number of sbjs who reached thresh performance on AB stim condition
CD_count = 0; %number of sbjs who reached thresh performance on CD stim condition
EF_count = 0; %number of sbjs who reached thresh performance on EF stim condition

group1 = [];
group2 = [];
group3 = [];

for sbj = 1 :size(Subject_List,1)
    load([pt2 filesep ['Behv_Performance_data_' Subject_List{sbj}(1:end-4) '.mat']])%load sbj data
    
    if Sbj_Behv_data.Behav_performance_group == 0 || Sbj_Behv_data.Behav_performance_group == 1
        group1_count = group1_count + 1;
        group1 = [group1 sbj];
    end
    if Sbj_Behv_data.Behav_performance_group == 2
        group2_count = group2_count + 1;
        group2 = [group2 sbj];
    end
    if Sbj_Behv_data.Behav_performance_group == 3
        group3_count = group3_count + 1;
        group3 = [group3 sbj];
    end
    
    if Sbj_Behv_data.AB_thresh_met == 1
        AB_count = AB_count + 1;
    end
    if Sbj_Behv_data.CD_thresh_met == 1
        CD_count = CD_count + 1;
    end
    if Sbj_Behv_data.EF_thresh_met == 1
        EF_count = EF_count + 1;
    end
end

disp(['Group 1 count = ' num2str(group1_count)])
disp(['Group 2 count = ' num2str(group2_count)])
disp(['Group 3 count = ' num2str(group3_count)])

disp([num2str(AB_count) ' subjects achieved threshold performance in AB stim condition.'])
disp([num2str(CD_count) ' subjects achieved threshold performance in CD stim condition.'])
disp([num2str(EF_count) ' subjects achieved threshold performance in EF stim condition.'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting RT

if plot_RTs == 1
    pt2 = 'C:\Users\ctgil\Documents\Probablistic_Learning\ProcessedData\Behavioral_Performance';
    
    RTs_Group1 = zeros(360,group1_count);
    RTs_Group2 = zeros(360,group2_count);
    RTs_Group3 = zeros(360,group3_count);
   
    G1_slopes = zeros(1,group1_count);
    G2_slopes = zeros(1,group2_count);
    G3_slopes = zeros(1,group3_count);
    G1_ints = zeros(1,group1_count);
    G2_ints = zeros(1,group2_count);
    G3_ints = zeros(1,group3_count);
    G1_RT_MNs = zeros(1,group1_count);
    G2_RT_MNs = zeros(1,group2_count);
    G3_RT_MNs = zeros(1,group3_count);
    
    Group1_count = 0;
    Group2_count = 0;
    Group3_count = 0;
    
    %Individual RT Plots
    for sbj = 1 :size(Subject_List,1)
        load([pt2 filesep ['Behv_Performance_data_' Subject_List{sbj}(1:end-4) '.mat']]);%load sbj data
        
             
        %Group RT data by Performance group
        if Sbj_Behv_data.Behav_performance_group == 1
            Group1_count = Group1_count + 1;
            RTs_Group1(:,Group1_count) = Sbj_Behv_data.RT;
        elseif Sbj_Behv_data.Behav_performance_group == 2
            Group2_count = Group2_count + 1;
            RTs_Group2(:,Group2_count) = Sbj_Behv_data.RT;
        else
            Group3_count = Group3_count + 1;
            RTs_Group3(:,Group3_count) = Sbj_Behv_data.RT;
        end
        
        fig1 = figure;
        scatter(1:length(Sbj_Behv_data.RT),Sbj_Behv_data.RT)
        title([Subject_List{sbj}(size(Subject_List{sbj},2)-8:end-4) ' - RT'])
        xlabel('Trial')
        ylabel('Time [s]')
        ylim([0 4])
        hl = lsline;
        hl.Color = 'r';
        B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
        Int = B(1);
        Slope = B(2);
        mn_RT = nanmean(Sbj_Behv_data.RT);
        
        dim = [0.2 0.5 0.3 0.3];
        str = {['Slope = ' num2str(Slope)],['\mu = ' num2str(mn_RT)]};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','red');
        
        
        Sbj_Behv_data.RT_slope = Slope;
        Sbj_Behv_data.mn_RT = mn_RT;
        
        %save data
        save([save_data_pt filesep 'Behv_Performance_data_' Subject_List{sbj}(1:end-4)],'Sbj_Behv_data');
        
        
        if Sbj_Behv_data.Behav_performance_group == 1
            G1_slopes(1,Group1_count) = Slope;
            G1_ints(1,Group1_count) = Int;
            G1_RT_MNs(1,Group1_count) = mn_RT;
        elseif Sbj_Behv_data.Behav_performance_group == 2
            G2_slopes(1,Group2_count) = Slope;
            G2_ints(1,Group2_count) = Int;
            G2_RT_MNs(1,Group2_count) = mn_RT;
        else
            G3_slopes(1,Group3_count) = Slope;
            G3_ints(1,Group3_count) = Int;
            G3_RT_MNs(1,Group3_count) = mn_RT;
        end
        
        if save_Sbj_RTs == 1
            savefig(fig1,[save_fig_pt filesep 'RTs' filesep sprintf('%s_RTs.fig',Subject_List{sbj}(1:end-4))]);
        end
        close all
    end
    
    %RT plots by Group
    if save_group_RTs == 1

        %Group 1
        fig1 = figure(1);
        for ii = 1:size(RTs_Group1,2)
            scatter(1:length(Sbj_Behv_data.RT),RTs_Group1(:,ii),'k')
            hold on;
        end
        title(['Group 1 - RT'])
        xlabel('Trial')
        ylabel('Time [s]')
        hline = refline(mean(G1_slopes),mean(G1_ints));
        hline.Color = 'r';
      
        dim = [0.14 0.61 0.3 0.3];
        str = {['Slope = ' num2str(mean(G1_slopes))],['Int =' num2str(mean(G1_ints))],['\mu = ' num2str(mean(G1_RT_MNs))]};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','red');
        
        savefig(fig1,[save_fig_pt filesep 'RTs' filesep 'Groups' filesep 'Group1_RTs.fig']);
        
        %Group 2
        fig2 =figure(2);
        for ii = 1:size(RTs_Group2,2)
            scatter(1:length(Sbj_Behv_data.RT),RTs_Group2(:,ii),'k')
            hold on;
        end
        title(['Group 2 - RT'])
        xlabel('Trial')
        ylabel('Time [s]')
        hline = refline(mean(G2_slopes),mean(G2_ints));
        hline.Color = 'r';
      
        dim = [0.14 0.61 0.3 0.3];
        str = {['Slope = ' num2str(mean(G2_slopes))],['Int =' num2str(mean(G2_ints))],['\mu = ' num2str(mean(G2_RT_MNs))]};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','red');
        
        savefig(fig2,[save_fig_pt filesep 'RTs' filesep 'Groups' filesep 'Group2_RTs.fig']);
        
        %Group 3
        fig3 = figure(3);
        for ii = 1:size(RTs_Group3,2)
            scatter(1:length(Sbj_Behv_data.RT),RTs_Group3(:,ii),'k')
            hold on;
        end
        title(['Group 3 - RT'])
        xlabel('Trial')
        ylabel('Time [s]')
        hline = refline(mean(G3_slopes),mean(G3_ints));
        hline.Color = 'r';
      
        dim = [0.14 0.61 0.3 0.3];
        str = {['Slope = ' num2str(mean(G3_slopes))],['Int =' num2str(mean(G3_ints))],['\mu = ' num2str(mean(G3_RT_MNs))]};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','red');
        
        savefig(fig3,[save_fig_pt filesep 'RTs' filesep 'Groups' filesep 'Group3_RTs.fig']);
    end
    
    close all
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create table with subjects behavioral data


for sbj = 1 :size(Subject_List,1)
        load([pt2 filesep ['Behv_Performance_data_' Subject_List{sbj}(1:end-4) '.mat']]);%load sbj data
        non_resp_count(sbj) = numel(Sbj_Behv_data.non_resp_trials);
        early_resp_count(sbj) = numel(Sbj_Behv_data.early_resp_trials);
        non_resp_trials{sbj} = Sbj_Behv_data.non_resp_trials;
        early_resp_trials{sbj} = Sbj_Behv_data.early_resp_trials;
        Perf_group(sbj) = Sbj_Behv_data.Behav_performance_group;
        mn_RT(sbj) = Sbj_Behv_data.mn_RT;
        RT_slope(sbj) = Sbj_Behv_data.RT_slope;
end

T = table(Subject_List,Perf_group',non_resp_count',non_resp_trials',early_resp_count',early_resp_trials',mn_RT',RT_slope');%,non_resp_count',early_resp_count',mn_RT,RT_slope');
T.Properties.VariableNames{'Var2'} = 'Perf_Group';
T.Properties.VariableNames{'Var3'} = 'Non_Resp_count';
T.Properties.VariableNames{'Var4'} = 'Non_Resp_trls';
T.Properties.VariableNames{'Var5'} = 'Early_Resp_count';
T.Properties.VariableNames{'Var6'} = 'Early_Resp_trls';
T.Properties.VariableNames{'Var7'} = 'mn_RT';
T.Properties.VariableNames{'Var8'} = 'RT_Slope';
