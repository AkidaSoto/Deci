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
        
        
        %cfg.detrend = 'yes';
        bcfg.bpfreq = Deci.ICA.bpfreq;
        bcfg.bpfilter = 'no';
        bcfg.detrend = 'yes';
        bcfg.demean = 'yes';
        
        tempdata = ft_preprocessing(bcfg,data);
        
        %% Do ICA
        
        if Deci.ICA.do
            
            disp('---ICA---');
            display(' ')
            feedback = 'no';
            
            tempcfg = cfg;
            
            cfg.feedback = feedback;
            cfg.demean     = 'no';
            

            
            evalc('data_ica = ft_componentanalysis(cfg, tempdata)');
            data_ica     = rmfield(data_ica,'cfg');
            
            figure;
            cfg.component = [1:length(data_ica.label)];
            cfg.viewmode = 'component';
            
            clear cfg.method
            cfg.channel = 'all';
            
            cfg.component = [];
            
            comps = [data_ica.trial{:}];
            all = [tempdata.trial{:}];
            eyechan = all(ismember(tempdata.label,Deci.ICA.eog),:);
            
            for eye = 1:size(eyechan,1)
                for comp = 1:size(comps,1)
                    [compcorr, p] = corrcoef(eyechan(eye,:),comps(comp,:));
                    corr(eye,comp,1) = compcorr(1,2);
                    corr(eye,comp,2) = p(1,2);
                end
                
                ecomponent{eye} = find(abs(corr(eye,:,1)) >= Deci.ICA.cutoff);
            end
            
            refchan = all(ismember(tempdata.label,{'TP9' 'TP10'}),:);
            
            for ref = 1:size(refchan,1)
                for comp = 1:size(comps,1)
                    [compcorr, p] = corrcoef(refchan(ref,:),comps(comp,:));
                    corr(ref,comp,1) = compcorr(1,2);
                    corr(ref,comp,2) = p(1,2);
                end
                
                rcomponent{ref} = find(abs(corr(ref,:,1)) >= Deci.ICA.cutoff);
            end
            
            possiblecorrupt = max(abs(exp(zscore(cfg.unmixing'))),[],1) > 1100;
            possiblecorrupt = possiblecorrupt | [max(abs(1./exp(zscore(cfg.unmixing'))),[],1) > 1100];
            
            if ~Deci.ICA.Automatic
                disp('Manual ICA Rejection')
                cfg = [];
                cfg.component = [1:length(data_ica.label)];
                cfg.viewmode = 'vertical';
                %cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
                
                cfg.channelcolormap = zeros(8,3);
                
                cfg.channelcolormap(8,:) = [1 1 1];
                cfg.channelcolormap(7,:) = [0 1 1];
                cfg.channelcolormap(6,:) = [1 0 1];
                cfg.channelcolormap(5,:) = [1 1 0];
                cfg.channelcolormap(4,:) = [0 0 1];
                cfg.channelcolormap(3,:) = [0 1 0];
                cfg.channelcolormap(2,:) = [1 0 0];
                cfg.channelcolormap(1,:) = [0 0 0];
                
                cfg.colorgroups = ones(length(data_ica.label),1)+7;
                
                for col = 1:length(cfg.colorgroups)
                    eye  = ismember(col,unique([ecomponent{:}]));
                    ref =  ismember(col,unique([rcomponent{:}]));
                    cor =  ismember(col,find(possiblecorrupt));
                    
                    cfg.colorgroups(col) = find(ismember(cfg.channelcolormap,[eye ref cor],'rows'));
                    
                end

                disp(['Found ' num2str(length(find(cfg.colorgroups == 1))) ' eye-correlated components']);
                
                cfg.channel = 'all';
                fakeUI = figure;
                select_labels(fakeUI,[],sort(cellfun(@(c)c(10:end), data_ica.label, 'un', 0)));
                fakeUI.Visible =  'off';

%                 figure;
%                 ft_topoplotIC(cfg, data_ica)
                
                evalc('ft_databrowser(cfg,data_ica)');
                suptitle(Deci.SubjectList{subject_list});
                waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
                
                if isempty(fakeUI.UserData)
                    cfg.component = [];
                else
                    cfg.component = find(ismember(cellfun(@(c)c(10:end),data_ica.label, 'un', 0),fakeUI.UserData));
                end
                close(fakeUI)
                corr = [];
            else
                disp('Automated ICA Rejection')
                cfg.component = unique([component{:}])';
            end
            
            evalc('data = ft_componentanalysis(tempcfg, data)');
            
            disp(['Rejecting components: ' num2str(cfg.component')]);
            cfg.demean = 'yes';
            evalc('data = ft_rejectcomponent(cfg, data)');
        end


        %% Manual Trial Rejection
        
        
        
        if Deci.Art.Manual_Trial_Rejection
            cfg =[];
            cfg.method = 'trial';
            tcfg.toilim = [abs(nanmax(locks,[],2)/1000)+Deci.Art.crittoilim(1) abs(nanmin(locks,[],2)/1000)+Deci.Art.crittoilim(2)];
            
            cfg.viewmode = 'vertical';
            evalc('artf = ft_databrowser(cfg,ft_redefinetrial(tcfg,tempdata))');
            
            evalc('datacomp_rej = ft_rejectartifact(artf,ft_redefinetrial(tcfg,tempdata))');
            
            
            postart.locks = postart.locks(ismember(postart.trlnum,trlnum(logical(datacomp_rej.saminfo))),:);
            postart.events = postart.events(ismember(postart.trlnum,trlnum(logical(datacomp_rej.saminfo))),:);
            postart.trlnum = postart.trlnum(ismember(postart.trlnum,trlnum(logical(datacomp_rej.saminfo))));
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
        
        if ~isempty(Deci.Art.More)
            evalc('data = ft_preprocessing(Deci.Art.More,data)');
            disp('Additional Preprocessing');
        end
        
        %% Artifact Reject
        
        Deci.Art = Exist(Deci.Art,'ATR',true);
        
        if Deci.Art.ATR
            cfg =[];
            cfg.method = 'summary';
            cfg.layout    = Deci.Layout.eye; % specify the layout file that should be used for plotting
            cfg.eog = Deci.Art.eog;
            cfg.keepchannel = 'no';
            tcfg.toilim = [abs(nanmax(locks,[],2)/1000)+Deci.Art.crittoilim(1) abs(nanmin(locks,[],2)/1000)+Deci.Art.crittoilim(2)];
            cfg.channel = 'all';
            
            evalc('data_rej = ft_rejectvisual(cfg,ft_redefinetrial(tcfg,data))');
            
            
            if Deci.Art.ShowArt
                Summary_artifacts = ~data_rej.saminfo;
                SAcfg.trials = Summary_artifacts;
                SAcfg.viewmode = 'vertical';
                evalc('savedtrls = ft_databrowser(SAcfg,ft_selectdata(SAcfg,ft_redefinetrial(tcfg,data)))');
                
                evalc('datacomp_saved = ft_rejectartifact(savedtrls,ft_selectdata(SAcfg,ft_redefinetrial(tcfg,data)))');
                
                savedtrls = find(Summary_artifacts);
                savedtrls =  savedtrls(~datacomp_saved.saminfo);
                Summary_artifacts(savedtrls) = 0;
                data_rej.saminfo = ~Summary_artifacts;
                
            end

            
            postart.locks = postart.locks(ismember(postart.trlnum,trlnum(logical(data_rej.saminfo))),:);
            postart.events = postart.events(ismember(postart.trlnum,trlnum(logical(data_rej.saminfo))),:);
            postart.trlnum = postart.trlnum(ismember(postart.trlnum,trlnum(logical(data_rej.saminfo))));
            
            display(' ')
            disp('---Trial Summary Rejection Applied---')
            disp(['Rejected ' num2str(length(find(~ismember(postart.trlnum,trlnum(data_rej.saminfo))))) ' trials'])
            %disp(['Remaining ' num2str(length(postart.trlnum)) ' trials'])
            disp (['Remaining ' num2str([length(find(ismember(postart.trlnum,trlnum(data_rej.saminfo))))]) ' trials'])
            disp(['Remaining ' num2str([length(find(ismember(postart.trlnum,trlnum(data_rej.saminfo))))/length(postart.trlnum)]*100) '% trials'])
            display(' ')
            pause(.05);
            
        end


        %% Interpolation

if Deci.Art.ATR && any(~ismember(data.label(~ismember(data.label,data_rej.label)),cfg.eog))
    rej_chan = data.label(~ismember(data.label,data_rej.label));
    interp_chan = rej_chan(~ismember(rej_chan,cfg.eog));

    Deci.Art.interp.badchannel = interp_chan;

        if isfield(Deci.Art,'interp')
            Deci.Art.interp.method = 'spline';
            load('elec1010_neighb.mat','neighbours');
            Deci.Art.interp.neighbours = neighbours;
            
            
            if exist([Deci.SubjectList{subject_list} '.bvct']) == 2
                [elec.label, elec.elecpos] = CapTrakMake([Deci.Folder.Raw  filesep Deci.SubjectList{subject_list} '.bvct']);
            else
                elec = ft_read_sens('standard_1020.elc');
            end
            Deci.Art.interp.elec = elec;
            display('Laplace Interpolation Applied')
            
            nonrepairs.channel = data.label(~ismember(data.label,elec.label));
            nonrepairs = ft_selectdata(nonrepairs,data);
            [data_interp] = ft_channelrepair(Deci.Art.interp, data);
            
            data = ft_appenddata([],nonrepairs,data_interp);
            clear nonrepairs
        end
end
%% RT
        
        if ~isempty(Deci.Art.RT)
            
            if Deci.Art.RT.dominlength
                RT_Art = [locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))] < Deci.Art.RT.minlength;
                RT_Nans = isnan([locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))]);
            
                display(' ')
                disp('---Minimum Reaction Time Rejection Applied---')
                disp(['Rejected ' num2str(length(find(~ismember(postart.trlnum,trlnum(~RT_Art))))) ' min length trials'])
                disp(['Remaining ' num2str([length(find(ismember(postart.trlnum,trlnum(~RT_Art))))]) ' trials'])
                disp(['Remaining ' num2str([length(find(ismember(postart.trlnum,trlnum(~RT_Art))))/length(trlnum)]*100) '% trials'])
                display(' ')
                pause(.05);
            
                postart.locks = postart.locks(ismember(postart.trlnum,trlnum(~RT_Art)),:);
                postart.events = postart.events(ismember(postart.trlnum,trlnum(~RT_Art)),:);
                postart.trlnum = postart.trlnum(ismember(postart.trlnum,trlnum(~RT_Art)));
            end 
           
           
            % auto reject max length trials
           if Deci.Art.RT.domaxlength
           RT_Art = [locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))] > Deci.Art.RT.maxlength;
           RT_Nans = isnan([locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))]);
             
           
              display(' ')
             disp('---Maximum Reaction Time Rejection Applied---')
             disp(['Rejected ' num2str(length(find(~ismember(postart.trlnum,trlnum(~RT_Art))))) ' max length trials'])
             disp(['Remaining ' num2str([length(find(ismember(postart.trlnum,trlnum(~RT_Art))))]) ' trials'])
             disp(['Remaining ' num2str([length(find(ismember(postart.trlnum,trlnum(~RT_Art))))/length(trlnum)]*100) '% trials'])
             display(' ')
             pause(.05);
             
             postart.locks = postart.locks(ismember(postart.trlnum,trlnum(~RT_Art)),:);
             postart.events = postart.events(ismember(postart.trlnum,trlnum(~RT_Art)),:);
             postart.trlnum = postart.trlnum(ismember(postart.trlnum,trlnum(~RT_Art)));
           end
           
            %%% auto-reject trials based on 2std    
           if Deci.Art.RT.dotwostd
            RT_Art = [locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))] > mean([locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))])+2*std([locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))]) ...
                    | [locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))] < mean([locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))])-2*std([locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))]) ;
            
            RT_Nans = isnan([locks(:,Deci.Art.RT.locks(2)) - locks(:,Deci.Art.RT.locks(1))]);
            
            display(' ')
            disp('---Reaction Time Rejection Applied---')
            disp(['Rejected ' num2str(length(find(~ismember(postart.trlnum,trlnum(~RT_Art))))) ' 2*std length trials'])
            disp(['Remaining ' num2str([length(find(ismember(postart.trlnum,trlnum(~RT_Art))))]) ' trials'])
            disp(['Remaining ' num2str([length(find(ismember(postart.trlnum,trlnum(~RT_Art))))/length(trlnum)]*100) '% trials'])
            display(' ')
            pause(.05);
            
            postart.locks = postart.locks(ismember(postart.trlnum,trlnum(~RT_Art)),:);
            postart.events = postart.events(ismember(postart.trlnum,trlnum(~RT_Art)),:);
            postart.trlnum = postart.trlnum(ismember(postart.trlnum,trlnum(~RT_Art)));
           end
           
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
        
        
        
%         
%         if Deci.Art.ATR && any(~ismember(data.label(~ismember(data.label,data_rej.label)),cfg.eog))
%             rej_chan = data.label(~ismember(data.label,data_rej.label));
%             interp_chan = rej_chan(~ismember(rej_chan,cfg.eog));
%             
%             TempDeci = Deci;
%             TempDeci.Art.interp.missingchannel = interp_chan;
%             TempDeci.PP.More.channel = data.label(~ismember(data.label,interp_chan));
%             TempDeci.SubjectList = Deci.SubjectList(subject_list);
%             TempDeci.Step = 2;
%             TempDeci.Proceed = 0;
%             TempDeci.PCom               = false;                                                      % Activates Parallel Computing for PP and Analysis only
%             TempDeci.GCom               = false;
%             TempDeci.DCom               = false;
%             TempDeci.ICA.RankReduction = length(interp_chan);
%             
%             Deci_Backend(TempDeci);
%             
%             TempDeci.Step = 3;
%             Deci_Backend(TempDeci);
%         else
            
            data.locks = locks;
            data.events = events;
            data.trlnum = trlnum;
            data.postart = postart;
            mkdir([Deci.Folder.Artifact])
            
            Deci.Art = Exist(Deci.Art,'AppendNewArtifacts',false);
            
            if Deci.Art.AppendNewArtifacts
                olddata =  load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data','info');
                
                
                postart.locks = olddata.info.postart.locks(ismember(olddata.info.postart.trlnum,postart.trlnum),:);
                postart.events = olddata.info.postart.events(ismember(olddata.info.postart.trlnum,postart.trlnum),:);
                postart.trlnum = olddata.info.postart.trlnum(ismember(olddata.info.postart.trlnum,postart.trlnum));
                
                data = olddata.data;
                data.postart = postart;
                info = rmfield(data,'trial');
                
                save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data','info','-v7.3')
                
            else
                
                info = rmfield(data,'trial');
                save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data','info','-v7.3')
            end
            %save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} '_info'],'data','-v7.3')
%         end
%         

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
        
        if isfield(data,'unmixing')
        data = rmfield(rmfield(data,'unmixing'),'topolabel');
        end
        
        Deci.Art = Exist(Deci.Art,'AppendNewArtifacts',false);
        
        if Deci.Art.AppendNewArtifacts
        olddata =  load([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'info','-v7.3');
            
                  
        postart.locks = olddata.info.postart.locks(ismember(olddata.info.postart.trlnum,postart.trlnum),:);
        postart.events = olddata.info.postart.events(ismember(olddata.info.postart.trlnum,postart.trlnum),:);
        postart.trlnum = olddata.info.postart.trlnum(ismember(olddata.info.postart.trlnum,postart.trlnum));
        
        data = olddata;
        data.postart = postart;
        
        save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data','info','-v7.3')
        
        else
        
        info = rmfield(data,'trial');
        save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list}],'data','info','-v7.3')
        end
        %save([Deci.Folder.Artifact filesep Deci.SubjectList{subject_list} '_info'],'data','-v7.3')
end


disp('----------------------');
end