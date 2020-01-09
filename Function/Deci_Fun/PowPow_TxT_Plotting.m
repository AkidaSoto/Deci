function PowPow_TxT_Plotting(Deci,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by CTGill
% last edited 1/6/20
%
% PowPow_TxT_Plotting plots Power-Power correlations averaged across selected
% frequecys.
%
% Power-Power Spearman correlation coefficients are averaged across specified
% frequency bins and then plotted, giving a seperate [time x trial]
% plot for each frequency bin, for each subject, condition, and lock.
%
% If Plot_Grand_Average is selected, the average across all subjects for
% each plot is computed and plotted.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data

data{length(Deci.SubjectList),length(Deci.Plot.Extra.PowPow_TxT_Plotting.locks),length(Deci.Plot.Extra.PowPow_TxT_Plotting.cond),length(Deci.Plot.Extra.PowPow_TxT_Plotting.freq)} = [];

for  subject = 1:length(Deci.SubjectList)
    for locks = 1:length(Deci.Plot.Extra.PowPow_TxT_Plotting.locks)
        for conditions = 1:length(Deci.Plot.Extra.PowPow_TxT_Plotting.cond)
            for freqs = 1:length(Deci.Plot.Extra.PowPow_TxT_Plotting.freq)
                if Deci.Plot.Extra.PowPow_TxT_Plotting.Laplace
                    load([Deci.Folder.Analysis filesep 'PowPow_TxT' filesep 'Laplacian' filesep Deci.SubjectList{subject}  filesep Deci.Plot.Extra.PowPow_TxT_Plotting.locks{locks} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions} filesep [Deci.Plot.Extra.PowPow_TxT_Plotting.chan1{1} '_' Deci.Plot.Extra.PowPow_TxT_Plotting.chan2{1}] filesep Deci.Plot.Extra.TxTPlotting.freq{freqs} '.mat'],'PowPow_corr_data');
                else
                    load([Deci.Folder.Analysis filesep 'PowPow_TxT' filesep Deci.SubjectList{subject}  filesep Deci.Plot.Extra.PowPow_TxT_Plotting.locks{locks} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions} filesep [Deci.Plot.Extra.PowPow_TxT_Plotting.chan1{1} '_' Deci.Plot.Extra.PowPow_TxT_Plotting.chan2{1}] filesep Deci.Plot.Extra.TxTPlotting.freq{freqs} '.mat'],'PowPow_corr_data');
                end
                data{subject,locks,conditions,freqs} = PowPow_corr_data;
            end
        end
    end
end

%% plots for individual subjects
if Deci.Plot.Extra.PowPow_TxT_Plotting.save_individual_sbj_plots
    for  subject = 1:length(Deci.SubjectList)
        for locks = 1:length(Deci.Plot.Extra.PowPow_TxT_Plotting.locks)
            for conditions = 1:length(Deci.Plot.Extra.PowPow_TxT_Plotting.cond)
                
                if Deci.Plot.Extra.PowPow_TxT_Plotting.Laplace
                    mkdir([Deci.Folder.Plot filesep 'PowPow_Corr' filesep 'Laplacian' filesep Deci.SubjectList{subject} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.locks{locks} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions} filesep [Deci.Plot.Extra.PowPow_TxT_Plotting.chan1{1} '_' Deci.Plot.Extra.PowPow_TxT_Plotting.chan2{1}]]);
                else
                    mkdir([Deci.Folder.Plot filesep 'PowPow_Corr' filesep Deci.SubjectList{subject} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.locks{locks} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions} filesep [Deci.Plot.Extra.PowPow_TxT_Plotting.chan1{1} '_' Deci.Plot.Extra.PowPow_TxT_Plotting.chan2{1}]]);
                end
                    
                for freqs = 1:length(Deci.Plot.Extra.PowPow_TxT_Plotting.freq)
                    
                    fig = figure;
                    contourf([1:size(data{subject,locks,conditions,freqs}.rho,2)],data{subject,locks,conditions,freqs}.time,data{subject,locks,conditions,freqs}.rho,'LineStyle','none');
                    ylabel('Time [s]');
                    xlabel('Trials');
                    title({['PowPow_Corr - ' Deci.SubjectList{subject}(size(Deci.SubjectList{subject},2)-4:end) ' - ' Deci.Plot.Extra.PowPow_TxT_Plotting.freq{freqs}]; [Deci.Plot.Extra.PowPow_TxT_Plotting.locks{locks} ' - ' Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions}(1:2) Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions}(4:end) ' - ' [Deci.Plot.Extra.PowPow_TxT_Plotting.chan1{1} '\_' Deci.Plot.Extra.PowPow_TxT_Plotting.chan2{1}]]});
                    colorbar
                    
                    if Deci.Plot.Extra.PowPow_TxT_Plotting.Laplace
                        savefig(fig,[Deci.Folder.Plot filesep 'PowPow_Corr' filesep 'Laplacian' filesep Deci.SubjectList{subject} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.locks{locks} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions} filesep [Deci.Plot.Extra.PowPow_TxT_Plotting.chan1{1} '_' Deci.Plot.Extra.PowPow_TxT_Plotting.chan2{1}] filesep sprintf('%s_PowPow_Corr.fig',Deci.Plot.Extra.PowPow_TxT_Plotting.freq{freqs})]);
                    else
                        savefig(fig,[Deci.Folder.Plot filesep 'PowPow_Corr' filesep Deci.SubjectList{subject} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.locks{locks} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions} filesep [Deci.Plot.Extra.PowPow_TxT_Plotting.chan1{1} '_' Deci.Plot.Extra.PowPow_TxT_Plotting.chan2{1}] filesep sprintf('%s_PowPow_Corr.fig',Deci.Plot.Extra.PowPow_TxT_Plotting.freq{freqs})]);
                    end
                    
                    close all;
                end
            end
        end
    end
    
end

%% Grand Average Plots
if Deci.Plot.Extra.PowPow_TxT_Plotting.GrandAverage
    
    for locks = 1:length(Deci.Plot.Extra.PowPow_TxT_Plotting.locks)
        for conditions = 1:length(Deci.Plot.Extra.PowPow_TxT_Plotting.cond)
            for freqs = 1:length(Deci.Plot.Extra.PowPow_TxT_Plotting.freq)
                
                if Deci.Plot.Extra.PowPow_TxT_Plotting.Laplace
                    mkdir([Deci.Folder.Plot filesep 'PowPow_Corr' filesep 'Laplacian' filesep 'All_Sbj_Avg' filesep Deci.Plot.Extra.PowPow_TxT_Plotting.locks{locks} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions} filesep [Deci.Plot.Extra.PowPow_TxT_Plotting.chan1{1} '_' Deci.Plot.Extra.PowPow_TxT_Plotting.chan2{1}]]);
                else
                    mkdir([Deci.Folder.Plot filesep 'PowPow_Corr' filesep 'All_Sbj_Avg' filesep Deci.Plot.Extra.PowPow_TxT_Plotting.locks{locks} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions} filesep [Deci.Plot.Extra.PowPow_TxT_Plotting.chan1{1} '_' Deci.Plot.Extra.PowPow_TxT_Plotting.chan2{1}]]);
                end
                
                max_trialNum = 0;
                for subject = 1:length(Deci.SubjectList)
                    trial_num = size((data{subject,locks,conditions,freqs}.rho),2);
                    
                    if (trial_num > max_trialNum)
                        max_trialNum = trial_num;
                    end
                end
                
                rho_Sum = zeros(101,max_trialNum);
                for subject = 1:length(Deci.SubjectList)
                    temp_rho = data{subject,locks,conditions,freqs}.rho;
                    new_num_trials = max_trialNum;
                    num_col = size(PowPow_corr_data.rho,1);
                    [x, y] = meshgrid(1:size(temp_rho,2),1:size(temp_rho,1));
                    [xq, yq] = meshgrid(linspace(1,size(temp_rho,2),new_num_trials),linspace(1,size(temp_rho,1),num_col));
                    interp_rho{subject} = interp2(x,y,temp_rho,xq,yq);
                    rho_Sum = rho_Sum + interp_rho{subject};
                end
                
                rho_grand_avg = rho_Sum/length(Deci.SubjectList);
                
                
                fig = figure;
                contourf([1:size(rho_grand_avg,2)],PowPow_corr_data.time,rho_grand_avg,'LineStyle','none');
                ylabel('Time [s]');
                xlabel('Trials');
                title({['PowPow\_Corr - All Subject Avg - ' Deci.Plot.Extra.PowPow_TxT_Plotting.freq{freqs}]; [Deci.Plot.Extra.PowPow_TxT_Plotting.locks{locks} ' - ' Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions}(1:2) Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions}(4:end) ' - ' [Deci.Plot.Extra.PowPow_TxT_Plotting.chan1{1} '\_' Deci.Plot.Extra.PowPow_TxT_Plotting.chan2{1}]]});
                c = colorbar;
                c.Label.String = 'Spearman rho';
                c.Label.Rotation = -90;
                c.Label.Position = [3, .06, 0];

                if Deci.Plot.Extra.PowPow_TxT_Plotting.Laplace
                    savefig(fig,[Deci.Folder.Plot filesep 'PowPow_Corr' filesep 'Laplacian' filesep 'All_Sbj_Avg' filesep filesep Deci.Plot.Extra.PowPow_TxT_Plotting.locks{locks} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions} filesep filesep [Deci.Plot.Extra.PowPow_TxT_Plotting.chan1{1} '_' Deci.Plot.Extra.PowPow_TxT_Plotting.chan2{1}] filesep sprintf('All_Sbj_Avg_%s_PowPow_Corr.fig',Deci.Plot.Extra.PowPow_TxT_Plotting.freq{freqs})]);
                else
                    savefig(fig,[Deci.Folder.Plot filesep 'PowPow_Corr' filesep 'All_Sbj_Avg' filesep filesep Deci.Plot.Extra.PowPow_TxT_Plotting.locks{locks} filesep Deci.Plot.Extra.PowPow_TxT_Plotting.cond{conditions} filesep filesep [Deci.Plot.Extra.PowPow_TxT_Plotting.chan1{1} '_' Deci.Plot.Extra.PowPow_TxT_Plotting.chan2{1}] filesep sprintf('All_Sbj_Avg_%s_PowPow_Corr.fig',Deci.Plot.Extra.PowPow_TxT_Plotting.freq{freqs})]);
                end
                
                close all;
            end
        end
    end
end

end
