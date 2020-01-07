function Artifactor(Deci)

for subject_list = 1:length(Deci.SubjectList)
    disp('----------------------');
    if Deci.Art.do
        display(['Starting Artifactor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
        display(' ')
        pause(.05);
        
        data = [];
        cfg = [];
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
        
        if isfield(data,'condinfo')  %replacer starting 12/22, lets keep for ~4 months
            data.postart.locks = data.condinfo{1};
            data.postart.events = data.condinfo{2};
            data.postart.trlnum = data.condinfo{3};
            
            data.locks = data.preart{1};
            data.events = data.preart{2};
            data.trlnum = data.preart{3};
            
            data = rmfield(data,'condinfo');
            data = rmfield(data,'preart');
        end
        
       locks = data.locks;
       events = data.events;
       trlnum = data.trlnum;
       postart.locks = data.postart.locks;
       postart.events = data.postart.events;
       postart.trlnum = data.postart.trlnum;
       
        %% Do ICA
        disp('---ICA---');
        display(' ')
        feedback = 'no';
        
        cfg.feedback = feedback;
        cfg.demean     = 'no';
        
        evalc('data_ica = ft_componentanalysis(cfg, data)');
        data_ica     = rmfield(data_ica,'cfg');
        
        figure;
        cfg.component = [1:20];
        cfg.viewmode = 'component';
        
        clear cfg.method
        cfg.channel = 'all';
        
        cfg.component = [];
        
        comps = [data_ica.trial{:}];
        eyes = [data.trial{:}];
        eyechan = eyes(ismember(data.label,Deci.ICA.eog),:);
        
        for eye = 1:size(eyechan,1)
            for comp = 1:size(comps,1)
                [compcorr, p] = corrcoef(eyechan(eye,:),comps(comp,:));
                corr(eye,comp,1) = compcorr(1,2);
                corr(eye,comp,2) = p(1,2);
            end
            
            component{eye} = find(abs(corr(eye,:,1)) >= Deci.ICA.cutoff);
        end
        
        possiblecorrupt = max(abs(exp(zscore(cfg.unmixing'))),[],1) > 1100;
        possiblecorrupt = possiblecorrupt | [max(abs(1./exp(zscore(cfg.unmixing'))),[],1) > 1100];
        
        if ~Deci.ICA.Automatic
            disp('Manual ICA Rejection')
            cfg.component = [1:20];
            cfg.viewmode = 'component';
            cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
            
            cfg.channelcolormap = zeros(4,3);
            cfg.channelcolormap(3,:) = [1 0 0];
            cfg.channelcolormap(2,:) = [0 0 1];
            cfg.channelcolormap(1,:) = [1 0 1];
            
            cfg.colorgroups = ones(20,1)+3;
            cfg.colorgroups(unique([component{:}]),1) = cfg.colorgroups(unique([component{:}]),1) - 1;
            cfg.colorgroups(possiblecorrupt,1) = cfg.colorgroups(possiblecorrupt,1) - 2;
            
            disp(['Found ' num2str(length(find(cfg.colorgroups == 1))) ' eye-correlated components']);
            
            cfg.channel = 'all'; 
            fakeUI = figure;
            select_labels(fakeUI,[],sort(data_ica.label));
            fakeUI.Visible =  'off';
            evalc('ft_databrowser(cfg,data_ica)');
            suptitle(Deci.SubjectList{subject_list});
            waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
            
            if isempty(fakeUI.UserData)
                cfg.component = [];
            else
                cfg.component = find(ismember(data_ica.label,fakeUI.UserData));
            end
            close(fakeUI)
            corr = [];
        else
            disp('Automated ICA Rejection')
            cfg.component = unique([component{:}])';
        end
        
        disp(['Rejecting components: ' num2str(cfg.component')]);
        cfg.demean = 'yes';
        evalc('data = ft_rejectcomponent(cfg, data_ica)');
        
        %% Manual Trial Rejection

        if Deci.Art.Manual_Trial_Rejection
            cfg =[];
            cfg.method = 'trial';
            tcfg.toilim = [abs(nanmax(locks,[],2)/1000)+Deci.Art.crittoilim(1) abs(nanmin(locks,[],2)/1000)+Deci.Art.crittoilim(2)];

            cfg.viewmode = 'vertical';
            evalc('artf = ft_databrowser(cfg,ft_redefinetrial(tcfg,data))');
            
            evalc('datacomp_rej = ft_rejectartifact(artf,ft_redefinetrial(tcfg,data))');
            
            postart.locks = postart.locks(ismember(postart.trlnum,find(datacomp_rej.saminfo)),:);
            postart.events = postart.events(ismember(postart.trlnum,find(datacomp_rej.saminfo)),:);
            postart.trlnum = postart.trlnum(ismember(postart.trlnum,find(datacomp_rej.saminfo)));
%             
%             cfg = [];
%             cfg.trials = postart.trlnum;
%             data_copy = ft_selectdata(cfg,data);

            display(' ')
            display('---Manual Trial Rejection Applied---')
            display(['Rejected ' num2str(length(find(~logical(datacomp_rej.saminfo)))) ' trials'])
            display(['Remaining ' num2str(length(postart.trlnum)) ' trials'])
            display(['Remaining ' num2str([length(postart.trlnum)/length(trlnum)]*100) '% trials'])
            display(' ')
            pause(.05);
        end

        %% Artifact Reject
        cfg =[];
        cfg.method = 'summary';
        cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
        cfg.eog = Deci.Art.eog;
        cfg.keepchannel = 'yes';
        tcfg.toilim = [abs(nanmax(locks,[],2)/1000)+Deci.Art.crittoilim(1) abs(nanmin(locks,[],2)/1000)+Deci.Art.crittoilim(2)];
        cfg.channel = 'all';

        evalc('data_rej = ft_rejectvisual(cfg,ft_redefinetrial(tcfg,data))');

        postart.locks = postart.locks(ismember(postart.trlnum,find(data_rej.saminfo)),:);
        postart.events = postart.events(ismember(postart.trlnum,find(data_rej.saminfo)),:);
        postart.trlnum = postart.trlnum(ismember(postart.trlnum,find(data_rej.saminfo)));
        
        display(' ')
        disp('---Trial Summary Rejection Applied---')
        disp(['Rejected ' num2str(length(find(~logical(data_rej.saminfo)))) ' trials'])
        disp(['Remaining ' num2str(length(postart.trlnum)) ' trials'])
        disp(['Remaining ' num2str([length(postart.trlnum)/length(trlnum)]*100) '% trials'])
        display(' ')
        pause(.05);
        
        if ~isempty(Deci.Art.RT)
            RT_Art = [locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))] < Deci.Art.RT.minlength;
            RT_Nans = isnan([locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))]);
            
            postart.locks = postart.locks(ismember(postart.trlnum,find(~RT_Art)),:);
            postart.events = postart.events(ismember(postart.trlnum,find(~RT_Art)),:);
            postart.trlnum = postart.trlnum(ismember(postart.trlnum,find(~RT_Art)));
            
            display(' ')
            disp('---Reaction Time Rejection Applied---')
            disp(['Rejected ' num2str(length(find(RT_Art))) ' trials'])
            disp(['Remaining ' num2str(length(postart.trlnum)) ' trials'])
            disp(['Remaining ' num2str([length(postart.trlnum)/length(trlnum)]*100) '% trials'])
            display(' ')
            pause(.05);
        end
       
        
        if Deci.Art.AddComponents
%             cfg =[];
%             cfg.viewmode = 'vertical';
%             
%             scfg.trials = condinfo{3};
%             
%             data_ica     = rmfield(ft_componentanalysis(ica, data),'cfg');
%             data_comp = ft_selectdata(scfg,data_ica);
%             
%             tcfg.toilim = [abs(nanmax(condinfo{1},[],2)/1000)+Deci.Art.crittoilim(1) abs(nanmin(condinfo{1},[],2)/1000)+Deci.Art.crittoilim(2)];
%             
%             artf = ft_databrowser(cfg,ft_redefinetrial(tcfg,data_comp));
%             
%             datacomp_rej = ft_rejectartifact(artf,ft_redefinetrial(tcfg,data_comp));
%             
%             condinfo{1} = condinfo{1}(logical(datacomp_rej.saminfo),:);
%             condinfo{2} = condinfo{2}(logical(datacomp_rej.saminfo),:);
%             if length(condinfo) > 2
%                 condinfo{3} = condinfo{3}(logical(datacomp_rej.saminfo));
%             end
        end
        
        data.locks = locks;
        data.events = events;
        data.trlnum = trlnum;
        data.postart = postart;
        mkdir([Deci.Folder.Artifact])
        save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data','-v7.3')
        data = rmfield(data,'trial');
        save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} '_info'],'data','-v7.3')
        
    else
        disp('Skipping Artifactor');
        
        mkdir([Deci.Folder.Artifact])
        
        load([Deci.Folder.Preproc filesep Deci.SubjectList{subject_list} '.mat']);
        
        if isfield(data,'condinfo')  %replacer starting 12/22, lets keep for ~4 months
            data.postart.locks = data.condinfo{1};
            data.postart.events = data.condinfo{2};
            data.postart.trlnum = data.condinfo{3};
            
            data.locks = data.preart{1};
            data.events = data.preart{2};
            data.trlnum = data.preart{3};
            
            data = rmfield(data,'condinfo');
            data = rmfield(data,'preart');
        end
        
        data = rmfield(rmfield(data,'unmixing'),'topolabel');
        save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data','-v7.3')
        data = rmfield(data,'trial');
        save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} '_info'],'data','-v7.3')
    end
    
    
    disp('----------------------');
end