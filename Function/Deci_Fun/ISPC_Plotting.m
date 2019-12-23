function ISPC_Plotting(Deci,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by CTGill
% last edited 12/21/19
%
% ISPC_Plotting plots ISPC_time averaged across selected frequecys(see
% MikeXCohen, 2014).
%
% ISPC is first computed in sliding time segments within a trial and is 
% then computed for each trial. Subsequently, ISPC is averaged across 
% selected frequency bins and then plotted, giving a seperate [time x trial]
% plot for each frequency bin, for each subject, condition, and lock. 
%
% If Plot_Grand_Average is selected, the average across all subjects for
% each plot is computed and plotted. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots for individual subjects
data{length(Deci.SubjectList),length(Deci.Plot.Extra.ISPC_Plotting.locks),length(Deci.Plot.Extra.ISPC_Plotting.cond)} = [];

for  subject = 1:length(Deci.SubjectList)
    for locks = 1:length(Deci.Plot.Extra.ISPC_Plotting.locks)
        for conditions = 1:length(Deci.Plot.Extra.ISPC_Plotting.cond)
            
            mkdir([Deci.Folder.Plot filesep 'ISPC' filesep Deci.SubjectList{subject} filesep Deci.Plot.Extra.ISPC_Plotting.locks{locks} filesep Deci.Plot.Extra.ISPC_Plotting.cond{conditions} filesep [Deci.Analysis.Connectivity.chan1{1} '_' Deci.Analysis.Connectivity.chan2{1}]]);
            
            load([Deci.Folder.Analysis filesep 'ISPC' filesep Deci.SubjectList{subject}  filesep Deci.Plot.Extra.ISPC_Plotting.locks{locks} filesep Deci.Plot.Extra.ISPC_Plotting.cond{conditions} filesep [Deci.Analysis.Connectivity.chan1{1} '_' Deci.Analysis.Connectivity.chan2{1}] '.mat'],'ISPC_data');
            data{subject,locks,conditions} = ISPC_data;
                       
            for fh = 1:length(Deci.Analysis.Connectivity.freq)
                switch Deci.Analysis.Connectivity.freq{fh}
                    case 'theta'
                        HF{fh} = [4 8];
                    case 'beta'
                        HF{fh} = [12.5 30];
                    case 'alpha'
                        HF{fh} =[8 12.5];
                    case 'lowgamma'
                        HF{fh} =[30 60];
                end

                cfg = [];
                cfg.frequency = [HF{fh}];
                ISPC_freq_band_avg = ft_selectdata(cfg,ISPC_data);
                ISPC_freq_band_avg = squeeze(mean(ISPC_data.phase_sync,1));
                
                fig = figure;
                contourf([1:size(ISPC_freq_band_avg,2)],ISPC_data.times2save,ISPC_freq_band_avg,'LineStyle','none');
                ylabel(['Time [s]']);
                xlabel(['Trials']);
                title({['ISPC - ' Deci.SubjectList{subject}(size(Deci.SubjectList{subject},2)-4:end) ' - ' Deci.Analysis.Connectivity.freq{fh}]; [Deci.Analysis.LocksTitle{locks} ' - ' Deci.Plot.Extra.ISPC_Plotting.cond{conditions}(1:2) Deci.Plot.Extra.ISPC_Plotting.cond{conditions}(4:end) ' - ' [Deci.Analysis.Connectivity.chan1{1} '\_' Deci.Analysis.Connectivity.chan2{1}]]});
                colorbar
               
                savefig(fig,[Deci.Folder.Plot filesep 'ISPC' filesep Deci.SubjectList{subject} filesep Deci.Plot.Extra.ISPC_Plotting.locks{locks} filesep Deci.Plot.Extra.ISPC_Plotting.cond{conditions} filesep [Deci.Analysis.Connectivity.chan1{1} '_' Deci.Analysis.Connectivity.chan2{1}] filesep sprintf('%s_ISPC.fig',Deci.Analysis.Connectivity.freq{fh})]);
                
                close all;
            end           
        end
    end
end

%% Grand Average Plots 
if Deci.Plot.Extra.ISPC_Plotting.GrandAverage
    
    for locks = 1:length(Deci.Plot.Extra.ISPC_Plotting.locks)
        for conditions = 1:length(Deci.Plot.Extra.ISPC_Plotting.cond)
            
            mkdir([Deci.Folder.Plot filesep 'ISPC' filesep 'All_Sbj_Avg' filesep Deci.Plot.Extra.ISPC_Plotting.locks{locks} filesep Deci.Plot.Extra.ISPC_Plotting.cond{conditions} filesep [Deci.Analysis.Connectivity.chan1{1} '_' Deci.Analysis.Connectivity.chan2{1}]]);
            
            max_trialNum = 0;
            for subject = 1:length(Deci.SubjectList)
                trial_num = size((data{subject,locks,conditions}.phase_sync),3);
                
                if (trial_num > max_trialNum)
                    max_trialNum = trial_num;
                end
            end
            
            for subject = 1:length(Deci.SubjectList)
                temp_sync = data{subject,locks,conditions}.phase_sync;
                new_num_trials = max_trialNum;
                [x, y, z] = meshgrid(1:size(temp_sync,2),1:size(temp_sync,1),1:size(temp_sync,3));
                [xq, yq, zq] = meshgrid(linspace(1,size(temp_sync,2),size(temp_sync,2)),linspace(1,size(temp_sync,1),size(temp_sync,1)),linspace(1,size(temp_sync,3),new_num_trials));
                interp_sync{subject} = interp3(x,y,z,temp_sync,xq,yq,zq);
            end

%%%%%%%%%%%%%%%%%%%%%%%% Left off here

%             temp_avg_pow = 0;
%             for  i = 1:length(Deci.SubjectList)
%                 temp_avg_pow = temp_avg_pow + interp_pow{i};
%             end
%             avg_pow = temp_avg_pow/length(Deci.SubjectList);
            
            
            fig = figure(1);
            imagesc(avg_pow);set(gca,'YDir','normal');
            view([90 -90]);
            colorbar;
            xlabel('Time [ms]');
            ylabel('Trial');
            title({['All Sbj Avg - ' Deci.Plot.Extra.TxTPlotting.locks{locks}]; [Deci.Plot.Extra.TxTPlotting.cond{conditions}(1:2) Deci.Plot.Extra.TxTPlotting.cond{conditions}(4:end) '- ' Deci.Plot.Extra.TxTPlotting.channels{chans} '- ' Deci.Plot.Extra.TxTPlotting.freq{freqs}]});
            
            savefig(fig,[Deci.Folder.Plot filesep 'Txt_avg_FreqBand_pow' filesep filesep 'All_Sbj_Avg' filesep Deci.Plot.Extra.TxTPlotting.locks{locks} filesep Deci.Plot.Extra.TxTPlotting.cond{conditions} filesep Deci.Plot.Extra.TxTPlotting.channels{chans} filesep sprintf('All_Sbj_Avg_%s_powspctrm.fig',params.freq{freqs})]);
            
            close all;
        end
    end
    
end

end
