function [out] = dc_RLmodel(Deci,info,dat,params)
% dc_RLmodel is ran once per block per state (Reward Conditions vs Punishment
% Conditions) per subject.
% 3 types of models are currently 
% 1 - Intial expected value are left to freely vary
        %- 5 total Free Parameters (Q0BTE,Q0WTE,alphaP,alphaN,beta)
% 2 - Intial expected value are left to params.start, contains seperate
% alpha for negative and positive PE
        %- 3 total Free Parameters (alphaP,alphaN,beta)
% 3- Intial expected value are left to params.start, no seperate alpha
        %- 2 total Free Parameters (alpha,beta)

% Model fit works by using fmincon to find the local minima (highest probability of state | action)
%from the randomized starting free parameter. The model is rurun multiple
%times to find the global minima.
        
%out - outgoing parameters of best fit, cell array 1x5
        %QL - Free Parameters
        %Q - Expected value
        %Pe - Prediction Error
        %PseudoR - Model Fitness
        %AIC - Model comparison

%Deci - Configuration File parameters
        %.Folder.Analysis - Folder where all analysis data is saved
        %.SubjectList - List of subject
        %.Analysis.CondTitle - List of conditions

%info - Unique parameters for this instance of function run
        %.subject_list - current subject index
        %.Cond - current condition index

%dat - Output of step 3 from Deci, fieldtrip data file
        %Artifact trial rejected
        % Function is usually called from Analyzor line 172 so you can read
        % how Analyzor loads an example Step 3 datafile.
        
%params - Free variables depending on model
        %.State - Event code for all possible conditions
        %.Actions - Event code for all possible Actions
        %.Reward - Event code for all possible Rewards
        %.Value - Value for each rewardxstate combo, cell format
        %.Reps - number of iterations to run randomized starting free
        %parameters
        %.Start - Starting Intial Expected Value
        %.ModelNum - number array of models to run


%Ex:
%params.States = [20 21 23 24];
%params.Actions = [31 32];
%params.Reward  = [51 52];
%params.Value  = {[20 10] [10 0] [0 -10] [-10 -20]};

trialtypes = dat.events;

%Basic Data Maintence for finding block numbers
% Deci segments blocks by negative numbers
if any(any(trialtypes < 0))
    trialtypes(:,find(trialtypes(1,:) < 1)) = ceil(trialtypes(:,find(trialtypes(1,:) < 1)));
    blocknumbers = sort(unique(trialtypes(:,find(trialtypes(1,:) < 1))),'descend');
else
    blocknumbers = -1;
end

Values = find(ismember(params.States,trialtypes));

%Sorting data into blocks
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

% Loop for running models N times with randomized starting conditinos

LLE = nan([3 params.Reps]);
LLE2 = nan([3 params.Reps]);

try
for init = 1:params.Reps

     q(init,:) = params.Start;
    if ismember(1,params.ModelNum)
        Fit1.LB = [min([params.Value{:}]) min([params.Value{:}]) 0 0 1e-6];
        Fit1.UB = [max([params.Value{:}]) max([params.Value{:}]) 1 1 30];
        Fit1.init =rand(1,length(Fit1.LB)).*(Fit1.UB-Fit1.LB)+Fit1.LB;
        
        [Value{1,init}] = ...
            fmincon(@(x)  FeedbackQ([x(1) x(2)],Actions,Rewards,x(3),x(4),x(5)),...
            [Fit1.init],[],[],[],[],[Fit1.LB],[Fit1.UB],[],...
            optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
        
        [LLE(1,init),qs{1,init},P{1,init},PE{1,init},qblock{1,init},Peblock{1,init}] = FeedbackQ([Value{1,init}(1) Value{1,init}(2)],Actions,Rewards, Value{1,init}(3),Value{1,init}(4),Value{1,init}(5));
    
     LLE2(1,init) = aicbic(-LLE(1,init),[5]);
    end
        
    if ismember(2, params.ModelNum)
        Fit2.LB = [0 0 1e-6];
        Fit2.UB = [1 1 30];
        Fit2.init =rand(1,length(Fit2.LB)).*(Fit2.UB-Fit2.LB)+Fit2.LB;
        
                [Value{2,init}] = ...
            fmincon(@(x) FeedbackQ(q(init,:),Actions,Rewards,x(1),x(2),x(3)),...
            [Fit2.init],[],[],[],[],[Fit2.LB],[Fit2.UB],[],...
            optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
        
        [LLE(2,init),qs{2,init},P{2,init},PE{2,init},qblock{2,init},Peblock{2,init}] = FeedbackQ(q(init,:),Actions,Rewards, Value{2,init}(1),Value{2,init}(2),Value{2,init}(3));
         
     LLE2(2,init) = aicbic(-LLE(2,init),[3]);    
    end
        
    
    if ismember(3, params.ModelNum)
        Fit3.LB = [0 1e-6];
        Fit3.UB = [1 3];
        
        Fit3.init =rand(1,length(Fit3.LB)).*(Fit3.UB-Fit3.LB)+Fit3.LB;
       

        [Value{3,init}] = ...
            fmincon(@(x) SimpleQ(q(init,:),Actions,Rewards,x(1),x(2)),...
            [Fit3.init],[],[],[],[],[Fit3.LB],[Fit3.UB],[],...
            optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
        
        [LLE(3,init),qs{3,init},P{3,init},PE{3,init},qblock{3,init},Peblock{3,init}] = SimpleQ(q(init,:),Actions,Rewards, Value{3,init}(1),Value{3,init}(2));
    
     LLE2(3,init) = aicbic(-LLE(3,init),[2]);
    end
             
        
end
catch
    %sometimes the model crashes from the given randomized starting free
    %parameers, so we try to catch it. Still unsure why, but I usually just
    %rerun the code again if it happens.
    error(['init ' num2str(Fit1.init) ' ' num2str(Fit2.init) ' ' num2str(Fit3.init) ' ']);
end


[Best,I] = nanmin(LLE,[],2);
[Best2,I2] = nanmin(LLE2,[],2);

model = find(max(Best2) == Best2);
I2 = I2(max(Best2) == Best2);


QL = Value{model,I2};
Q = qs{model,I2};
Pe = PE{model,I2};
PseudoR = 1 - [-Best2(model)/[length(cat(1,Actions{:}))*log(.05)]];
AIC = Best2(model);

Q_all = qblock{I2};
Pe_all = Peblock{I2};



%out = {QL Q Pe PseudoR AIC};

out = nan(size(Pe));




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