function Deci_Corr(Deci,info,freq,params)
%% Load Data
freq.dimord = 'rpt_chan_freq_time';


for brains = 1:length(params.Brain)

    switch params.Brain{brains}
        case 'Magnitude'
            brain = abs(freq.fourierspctrm).^2;

            if params.bsl.do
                Bsl = brain;
                toi = freq.time >= round(params.bsl.time(1),4) & freq.time <= round(params.bsl.time(2),4);

                Bsl = nanmean(Bsl(:,:,:,toi),4);
                Bsl = repmat(Bsl,[1 1 1 size(brain,4)]);

                switch params.bsl.type
                    case 'none'
                    case 'absolute'
                        brain =  brain - Bsl;
                    case 'relative'
                        brain=  brain ./ Bsl;
                    case 'relchange'
                        brain = ( brain - Bsl) ./ Bsl;
                    case 'db'
                        brain = 10*log10( brain ./ Bsl);
                end

            end

        case 'Phase'
            brain = ang2rad(angle(freq.fourierspctrm));
    end

    toi = freq.time;

    if ~strcmpi(params.type,'regress')
        for behaviors = 1:length(params.Behavior)

            behavior = load([Deci.Folder.Analysis filesep 'Extra' filesep params.Behavior{behaviors} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.CondTitle{info.Cond}]);

            Type = fieldnames(behavior);

            behavior = behavior.(Type{1});

            for foi = 1:length(freq.freq)

                for ti = 1:length(toi)

                    b_time = brain(:,:,foi,toi(ti) == toi);

                    if ~isequal(size(b_time),size(behavior))
                        behavior = behavior';
                    end

                    parameter{behaviors} = behavior;

                    switch params.Brain{brains}
                        case 'Magnitude'

                            if strcmpi(params.type,'spearman')
                                params = Exist(params,'logsig',false);
                                if params.logsig
                                    b_time = logsig(b_time);
                                end

                                params = Exist(params,'type',[]);

                                [R(1,foi,ti), P(1,foi,ti)] = corr(b_time,parameter{behaviors},'Type','Spearman');
                            elseif strcmpi(params.type,'linear')

                                [y] = polyfit(b_time,parameter{behaviors},1);
                                R(1,foi,ti) = y(1);
                                P(1,foi,ti) = y(2);

                            else
                                [r,p] = corrcoef(b_time,parameter{behaviors});
                                R(1,foi,ti) = r(1,2);
                                P(1,foi,ti) = p(1,2);
                            end


                        case 'Phase'
                            [R(1,foi,ti),P(1,foi,ti)] =  circ_corrcl(b_time, parameter{behaviors});

                    end
                end
            end
        end

    else

        for behaviors = 1:length(params.Behavior)

            behavior = load([Deci.Folder.Analysis filesep 'Extra' filesep params.Behavior{behaviors} filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.CondTitle{info.Cond}]);

            Type = fieldnames(behavior);

            behavior = behavior.(Type{1});

            if ~isequal(size(b_time,1),size(behavior,1))
                behavior = behavior';
            end

            parameter{behaviors} = behavior;
        end

        for foi = 1:length(freq.freq)

            for ti = 1:length(toi)

                b_time = brain(:,:,foi,toi(ti) == toi);

                switch params.Brain{brains}
                    case 'Magnitude'
                        mdl = fitlm(parameter,b_time);
                    case 'Phase'
                        [beta,R2,p] = CircularRegression(parameter,b_time);
                end

            end
        end



    end

    %         if any(any(arrayfun(@(c) iscomplex(c),R)))
    %            k = 0;
    %         end

    extracorr.label = freq.label;
    extracorr.freq = freq.freq;
    extracorr.time = freq.time(toi);
    extracorr.dimord =  'chan_freq_time';

    extracorr.powspctrm = R;
    R = extracorr;

    extracorr.powspctrm = P;
    P = extracorr;

    mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Corr' filesep params.Brain{brains} '_' params.Behavior{behaviors}  filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond}])
    save([Deci.Folder.Analysis filesep 'Extra' filesep 'Corr' filesep params.Brain{brains} '_' params.Behavior{behaviors}  filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.LocksTitle{info.Lock} filesep Deci.Analysis.CondTitle{info.Cond} filesep info.Channels{info.ChanNum}],'R','P');
    clear R P

end

end

%%
