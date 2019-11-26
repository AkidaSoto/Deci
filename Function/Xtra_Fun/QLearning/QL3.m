function [out] = QL3(Deci,info,dat,params)
%params.States = [20 21 23 24];
%params.Actions = [31 32];
%params.Reward  = [51 52];
%params.Value  = {[20 10] [10 0] [0 -10] [-10 -20]};


maxt = max(sum(ismember(dat.preart{2},Deci.Analysis.Conditions{info.Cond}),2));
info.alltrials = find(sum(ismember(dat.preart{2},Deci.Analysis.Conditions{info.Cond}),2) == maxt);

%ignore all locks that are missing
minamountofnans = min(mean(isnan(dat.preart{1}(info.alltrials,:)),2));
info.allnonnans = mean(isnan(dat.preart{1}(info.alltrials,:)),2) == minamountofnans;% & ~isnan(mean(condinfo{2},2));

trialtypes = dat.preart{2}(info.alltrials(info.allnonnans),:);


%Basic Data Maintence for finding block numbers
if any(any(trialtypes < 0))
    trialtypes(:,find(trialtypes(1,:) < 1)) = ceil(trialtypes(:,find(trialtypes(1,:) < 1)));
    blocknumbers = sort(unique(trialtypes(:,find(trialtypes(1,:) < 1))),'descend');
else
    blocknumbers = -1;
end

Values = find(ismember(params.States,trialtypes));

%Sorting data into blcks
for blk = 1:length(blocknumbers)
    
    blkmrk = sum(ismember(trialtypes,[blocknumbers(blk)]),2) ==  1;
    
    Actmrk = logical(sum(ismember(trialtypes(blkmrk,:),[params.Actions]),1));
    Actions{blk} = trialtypes(blkmrk,Actmrk);
    
    Rewmrk = logical(sum(ismember(trialtypes(blkmrk,:),[params.Reward]),1));
    Rewards{blk} = trialtypes(blkmrk,Rewmrk);

    for rew = 1:length(params.Reward)
        Rewards{blk}(Rewards{blk} == params.Reward(rew)) = params.Value{Values}(rew);
    end

end

% Loop for running models 10 times with randomized starting conditinos

try
for init = 1:params.Reps
>>>>>>> a7e8eaa4403c906a1e5345eaacc930e82c6361df
        
        Fit1.LB = [min([params.Value{:}]) min([params.Value{:}]) 0 0 1e-6];
        Fit1.UB = [max([params.Value{:}]) max([params.Value{:}]) 1 1 30];
        Fit1.init =rand(1,length(Fit1.LB)).*(Fit1.UB-Fit1.LB)+Fit1.LB;
        
        Fit2.LB = [0 0 1e-6];
        Fit2.UB = [1 1 30];
        Fit2.init =rand(1,length(Fit2.LB)).*(Fit2.UB-Fit2.LB)+Fit2.LB;
        
        Fit3.LB = [0 1e-6];
        Fit3.UB = [1 3];
        
        Fit3.init =rand(1,length(Fit3.LB)).*(Fit3.UB-Fit3.LB)+Fit3.LB;
        q(init,:) = params.Start;

        [Value{1,init}] = ...
            fmincon(@(x)  FeedbackQ([x(1) x(2)],Actions,Rewards,x(3),x(4),x(5)),...
            [Fit1.init],[],[],[],[],[Fit1.LB],[Fit1.UB],[],...
            optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
        
        [LLE(1,init),qs{1,init},P{1,init},PE{1,init},qblock{1,init},Peblock{1,init}] = FeedbackQ([Value{1,init}(1) Value{1,init}(2)],Actions,Rewards, Value{1,init}(3),Value{1,init}(4),Value{1,init}(5));
        
        [Value{2,init}] = ...
            fmincon(@(x) FeedbackQ(q(init,:),Actions,Rewards,x(1),x(2),x(3)),...
            [Fit2.init],[],[],[],[],[Fit2.LB],[Fit2.UB],[],...
            optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
        
        [LLE(2,init),qs{2,init},P{2,init},PE{2,init},qblock{2,init},Peblock{2,init}] = FeedbackQ(q(init,:),Actions,Rewards, Value{2,init}(1),Value{2,init}(2),Value{2,init}(3));
         

        [Value{3,init}] = ...
            fmincon(@(x) SimpleQ(q(init,:),Actions,Rewards,x(1),x(2)),...
            [Fit3.init],[],[],[],[],[Fit3.LB],[Fit3.UB],[],...
            optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
        
        [LLE(3,init),qs{3,init},P{3,init},PE{3,init},qblock{3,init},Peblock{3,init}] = SimpleQ(q(init,:),Actions,Rewards, Value{3,init}(1),Value{3,init}(2));
        
        LLE2(:,init) = aicbic(-LLE(:,init),[5 3 2]);
        
        
end
catch
    error(['init ' num2str(Fit1.init) ' ' num2str(Fit2.init) ' ' num2str(Fit3.init) ' ']);
end


[Best,I] = min(LLE,[],2);
[Best2,I2] = min(LLE2,[],2);


for m = 1:size(Best,1)
    %out{1,m,blk} = Value{m,I(m,1)};
    QL{m} = Value{m,I2(m,1)};
    Q{m} = qs{m,I2(m,1)};
    Pe{m} = PE{m,I2(m,1)};
    PseudoR{m} = 1 - [-Best(m)/[length(cat(1,Actions{:}))*log(.05)]];
    AIC{m} = Best2(m);
    
    Q_all{m} = qblock{m};
    Pe_all{m} = Peblock{m};
    
    if Deci.Analysis.ApplyArtReject
        Q{m} = Q{m}(find(ismember(dat.preart{3}(info.alltrials(info.allnonnans)),dat.condinfo{3})));
        Pe{m} = Pe{m}(find(ismember(dat.preart{3}(info.alltrials(info.allnonnans)),dat.condinfo{3})));
    end
    
end



out = {QL Q Pe PseudoR AIC};


warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'QL' filesep Deci.SubjectList{info.subject_list} ])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'QL' filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.CondTitle{info.Cond}],'QL');

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Q' filesep Deci.SubjectList{info.subject_list}])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'Q' filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.CondTitle{info.Cond}],'Q');

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Pe' filesep Deci.SubjectList{info.subject_list}  filesep  ])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'Pe' filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.CondTitle{info.Cond}],'Pe');

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Q_all' filesep Deci.SubjectList{info.subject_list}])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'Q_all' filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.CondTitle{info.Cond}],'Q_all');

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'Pe_all' filesep Deci.SubjectList{info.subject_list}  filesep  ])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'Pe_all' filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.CondTitle{info.Cond}],'Pe_all');

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'PseudoR' filesep Deci.SubjectList{info.subject_list} ])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'PseudoR' filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.CondTitle{info.Cond}],'PseudoR');

mkdir([Deci.Folder.Analysis filesep 'Extra' filesep 'AIC' filesep Deci.SubjectList{info.subject_list} ])
save([Deci.Folder.Analysis filesep 'Extra' filesep 'AIC' filesep Deci.SubjectList{info.subject_list} filesep Deci.Analysis.CondTitle{info.Cond}],'AIC');

    function [LLE,qout,P,PE,qblock,PEblock] = SimpleQ(q,a,r,alp,beta)
        
        P = [];
        PE = [];
        qout = [];
        qin = q;
        
        for block = 1:length(a)
            
            actions = params.Actions;
            qblock{block} = [];
            PEblock{block} = [];
            q = qin;
            
            for Act = 1:length(a{block})
                
                which = find(a{block}(Act) == actions);
                
                qout(end+1) = q(Act,which);
                qblock{block}(end+1) =  qout(end);
                
                PE(end+1) = (r{block}(Act) - q(Act,which));
                PEblock{block}(end+1) =  PE(end);
               
                P(end+1) = exp(q(Act,which)/beta)/[exp(q(Act,which)/beta) + exp(q(Act,find(a{block}(Act) ~= actions))/beta)];
                
                if isnan(P(end))
                    P(end) = q(Act,which) > q(Act,find(a{block}(Act) ~= actions));
                end
                
                q(Act+1,which) = q(Act,which) + alp*PE(end);
                q(Act+1,find(a{block}(Act) ~= actions))  = q(Act,find(a{block}(Act) ~= actions));
            end
        end
        
        LLE = -sum(P);
    end

    function [LLE,qout,P,PE,qblock,PEblock] = FeedbackQ(q,a,r,alpP,alpN,beta)
        
         P = [];
         PE = [];
         qout = [];
         qin = q;
         
         
         for block = 1:length(a)
             
             actions = params.Actions;
             qblock{block} = [];
             PEblock{block} = [];
             q = qin;

             for Act = 1:length(a{block})
                 
                 which = find(a{block}(Act) == actions);
                 
                 qout(end+1) = q(Act,which);
                 qblock{block}(end+1) =  qout(end);
                                  
                 PE(end+1) = (r{block}(Act) - q(Act,which));
                 PEblock{block}(end+1) =  PE(end);
                 
                 P(end+1) = exp(q(Act,which)/beta)/[exp(q(Act,which)/beta) + exp(q(Act,find(a{block}(Act) ~= actions))/beta)];
                 
                 if isnan(P(end))
                     P(end) = q(Act,which) > q(Act,find(a{block}(Act) ~= actions));
                 end
                 
                 
                 if PE(end) > 0
                     q(Act+1,which) = q(Act,which) + alpP*PE(end);
                 else
                     q(Act+1,which) = q(Act,which) + alpN*PE(end);
                 end
                 q(Act+1,find(a{block}(Act) ~= actions))  = q(Act,find(a{block}(Act) ~= actions));

                 
                 
             end
         end
        
        LLE = -sum(P);
        
    end

    function [LLE,qout,P,PE,qblock,PEblock] = InteractiveQ(q,a,r,alpP,alpN,beta)
        P = [];
        PE = [];
        qout = [];
        qin = q;
        
        
        for block = 1:length(a)
            
            actions = params.Actions;
            qblock{block} = [];
            PEblock{block} = [];
            q = qin;
            
            for Act = 1:length(a{block})
                
                which = find(a{block}(Act) == actions);
                
                PE(end+1) = (r{block}(Act) -q(Act,which) );
               
                if PE(end) > 0
                    q(Act+1,which) = q(Act,which) + alpP*PE(end);
                else
                    q(Act+1,which) = q(Act,which) + alpN*PE(end);
                end
                
                q(Act+1,find(a{block}(Act) ~= actions))  = q(Act,find(a{block}(Act) ~= actions));
                
                P(end+1) = exp(beta*q(Act+1,which))/sum(exp(beta*q(Act+1,:)));
                qout(end+1) = q(Act+1,which);
                
                qblock{block}(end+1) =  qout(end);
                PEblock{block}(end+1) =  PE(end);
            end
        end
        
        LLE = -nansum(log(P));
    end

    function [LLE,qout,P,PE,qblock,PEblock] = ModelQ(q,a,r,alpP,alpN,betaP,betaR)
        
        actions = sort(unique(a));
        
        for Act = 1:length(a)
            
            which = find(a(Act) == actions);
            
            PE(Act) = (r(Act) - q(Act,which));
            
            if PE(Act) > 0
                q(Act+1,which) = q(Act,which) + alpP*PE(Act);
            else
                q(Act+1,which) = q(Act,which) + alpN*PE(Act);
            end
            
            q(Act+1,find(a(Act) ~= actions))  = q(Act,find(a(Act) ~= actions));
            
            if r(Act) > 0
                P(Act) = exp(betaR*q(Act+1,which))/sum(exp(betaR*q(Act+1,:)));
            elseif r(Act) < 0
                P(Act) = exp(betaP*q(Act+1,which))/sum(exp(betaP*q(Act+1,:)));
            end
            
            qout(Act) = q(Act+1,which);

        end
        
         LLE = -nansum(log(P));
    end
end