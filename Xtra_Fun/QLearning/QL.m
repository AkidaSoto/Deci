function QL(Deci,info,freq,params)
%params.States = [20 21 23 24];
%params.Actions = [31 32];
%params.Reward  = [51 52];
%params.Value  = {[20 10] [10 0] [0 -10] [-10 -20]};

trialtypes = freq.condinfo{2}(info.alltrials,:);


%Basic Data Maintence for finding block numbers
if any(any(trialtypes < 0))
    blocknumbers = trialtypes(:,end);
else
    blocknumbers = -1;
end

Values = any(ismember(trialtypes',params.States)');

%Sorting data into blcks
for blk = 1:length(blocknumbers)
    
    
    Actions{blk} = zeros([size(trialtypes,1) 1]);
    Rewards{blk} = zeros([size(trialtypes,1) 1]);
    
    for dec = 1:length(params.Actions)
        maxt = max(sum(ismember(trialtypes,[params.Actions(dec) blocknumbers(blk)]),2));
        Actions{blk}(sum(ismember(trialtypes,[params.Actions(dec) blocknumbers(blk)]),2) == maxt) = params.Actions(dec);
        Rewards{blk}(sum(ismember(trialtypes,[params.Actions(dec) blocknumbers(blk)]),2) == maxt) = trialtypes(sum(ismember(trialtypes,[params.Actions(dec) blocknumbers(blk)]),2) == maxt,find(any(ismember(trialtypes',params.Reward)')));
    end
    
    for rew = 1:length(params.Reward)
        Rewards{blk}(Rewards{blk} == params.Reward(rew)) = params.Value{Values}(rew);
    end
    
     Rewards{blk} =  Rewards{blk}(Actions{blk} ~= 0);
     Actions{blk} =  Actions{blk}(Actions{blk} ~= 0);
end

% Loop for running models 10 times with randomized starting conditinos
for init = 1:10
    Fit1.LB = [0 1e-6];
    Fit1.UB = [1 30];
    
    Fit1.init =rand(1,length(Fit1.LB)).*(Fit1.UB-Fit1.LB)+Fit1.LB;
    
    Fit2.LB = [0 0 1e-6];
    Fit2.UB = [1 1 30];
    Fit2.init =rand(1,length(Fit2.LB)).*(Fit2.UB-Fit2.LB)+Fit2.LB;
    
    Fit3.LB = [0 0 0 0 1e-6];
    Fit3.UB = [1 1 1 1 30];
    Fit3.init =rand(1,length(Fit3.LB)).*(Fit3.UB-Fit3.LB)+Fit3.LB;
    
    Fit4.LB = [0 0 1e-6 1e-6];
    Fit4.UB = [1 1 30 30];
    Fit4.init =rand(1,length(Fit4.LB)).*(Fit4.UB-Fit4.LB)+Fit4.LB;
    
    %q = ones([1 length(unique(Actions{blk}))]);
    q = zeros([1 length(unique(Actions{blk}))]);
    %q = rand([1 length(unique(Actions{blk}))]);
    
    % run through each model with the randomized starting values
    for blk = 1:length(blocknumbers)
        
            [res(1),LLE(1)] = ...
                fmincon(@(x) SimpleQ(q,Actions{blk},Rewards{blk},x(1),x(2)),...
                [Fit1.init],[],[],[],[],[Fit1.LB],[Fit1.UB],[],...
                optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
            
            [res(2),LLE(2)] = ...
                fmincon(@(x) FeedbackQ(q,Actions{blk},Rewards{blk},x(1),x(2),x(3)),...
                [Fit2.init],[],[],[],[],[Fit2.LB],[Fit2.UB],[],...
                optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
            
            [res(3),LLE(3)] = ...
                fmincon(@(x) InteractiveQ(q,Actions{blk},Rewards{blk},x(1),x(2),x(3),x(4),x(5)),...
                [Fit3.init],[],[],[],[],[Fit3.LB],[Fit3.UB],[],...
                optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
            
            
            [res(4),LLE(4)] = ...
                fmincon(@(x) ModelQ(q,Actions{blk},Rewards{blk},x(1),x(2),x(3),x(4)),...
                [Fit4.init],[],[],[],[],[Fit4.LB],[Fit4.UB],[],...
                optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
            
            
            Result.Eta(init,1) = res1(1);
            Result.Beta(init,1) = res1(2);
            Result.Lik(init,1) = lik1;
            
            
            Result.Eta1(init,2) = res2(1);
            Result.Eta2(init,2) = res2(2);
            Result.Beta(init,2) = res2(3);
            Result.Lik(init,2) = lik2;
            
            Result.Eta1(init,3) = res3(1);
            Result.Eta2(init,3) = res3(2);
            Result.Beta(init,3) = res3(3);
            Result.Lik(init,3) = lik3;
            
            Result.Eta(init,4) = res4(1);
            Result.Beta(init,4) = res4(2);
            Result.Lik(init,4) = lik4;
            
            Result.Lik(init,:)  % to view progress so far

    end
end

% I have not reviewed the code beyond this point. you may find that you can
% just scrap everything past this point if needed

for l = 1:size(Result.Lik,2)
    BestFit = find(Result.Lik(:,l) == min(Result.Lik(:,l)),1,'first');
    
    if l == 1
        [~,qs,P] = SimpleQ(q,s,a,r, Result.Eta(BestFit,l), Result.Beta(BestFit,l));
    elseif l == 2
        [~,qs,P] = FeedbackQ(q,s,a,r, Result.Eta1(BestFit,l),Result.Eta2(BestFit,l), Result.Beta(BestFit,l));
    elseif l == 3
        [~,qs,P] = InteractiveQ(q,s,a,r, Result.Eta1(BestFit,l),Result.Eta2(BestFit,l), Result.Beta(BestFit,l));
    elseif l == 4
        [~,qs,P] = ModelQ(q,s,a,r, Result.Eta(BestFit,l), Result.Beta(BestFit,l));
    end
    
    Pb = Simu(P);
    
    trew = [];
    
    for t = 1:1000
        rew = [];
        for i = 16:23
            rew(i,:) = Pb(find(s == i));
        end
        rew = rew(16:23,:);
        
        trew(:,:,t) = rew;
    end
    hold on
    plot(mean(mean(trew,3),1))
end

legend(['Actual Data'],...
    ['SimpleQ (' num2str(Result.Lik(BestFit,1)) ')' ], ...
    ['FeedbackQ(' num2str(Result.Lik(BestFit,2)) ')' ],...
    ['InteractiveQ (' num2str(Result.Lik(BestFit,3)) ')' ],...
    ['ModelQ (' num2str(Result.Lik(BestFit,4)) ')' ]);


    function [LLE,qout,P] = SimpleQ(q,a,r,alp,beta)
        
        %[qout,PE,P,LLE] =
        
        lik = 0;
        
        actions = sort(unique(a));
        
        for Act = 1:length(a)
            
            which = find(a(Act) == actions);
            
            PE(Act) = (r(Act) - q(Act,which));
            q(Act+1,which) = q(Act,which) + alp*PE(Act);
            q(Act+1,find(a(Act) ~= actions))  = q(Act,find(a(Act) ~= actions));
            
            P(Act) = exp(beta*q(Act+1,which))/sum(exp(beta*q(Act+1,:)));
            qout(Act) = q(Act+1,which);
            
            
        end
        
        LLE = -nansum(log(P));
        
    end

    function [LLE,qout,P] = FeedbackQ(q,a,r,alpP,alpN,beta)
        
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
            
            P(Act) = exp(beta*q(Act+1,which))/sum(exp(beta*q(Act+1,:)));
            qout(Act) = q(Act+1,which);
            
        end
        
        LLE = -nansum(log(P));
        
    end

    function [LLE,qout,P] = InteractiveQ(q,a,r,alpPP,alpNP,alpPR,alpNR,beta)
        
        actions = sort(unique(a));
        
        for Act = 1:length(a)
            
            which = find(a(Act) == actions);
            
            PE(Act) = (r(Act) - q(Act,which));
            
            if PE(Act) > 0 && [r(Act) > 0 || [r(Act) == 0 && mean(unique(r)) < 0]]
                q(Act+1,which) = q(Act,which) + alpPR*PE(Act);
            elseif PE(Act) > 0 && r(Act) < 0 || [r(Act) == 0 && mean(unique(r)) > 0]
                q(Act+1,which) = q(Act,which) + alpPP*PE(Act);
            elseif PE(Act) < 0 && r(Act) > 0 || [r(Act) == 0 && mean(unique(r)) < 0]
                q(Act+1,which) = q(Act,which) + alpNR*PE(Act);
            elseif PE(Act) < 0 && r(Act) < 0 || [r(Act) == 0 && mean(unique(r)) < 0]
                q(Act+1,which) = q(Act,which) + alpNP*PE(Act);
            end
            
            q(Act+1,find(a(Act) ~= actions))  = q(Act,find(a(Act) ~= actions));
            
            P(Act) = exp(beta*q(Act,which))/sum(exp(beta*q(Act,:)));
            qout(Act) = q(Act+1,which);
            
        end
        
        LLE = -nansum(log(P));
    end

    function [LLE,qout,P] = ModelQ(q,a,r,alpP,alpN,betaP,betaR)
        
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

    function B = Simu(A)
        
        for act = 1:length(A)
            if ~isnan(A(act))
                
                B(act) = rand < A(act);
                
            else
                B(act) = 0;
            end
        end
        
    end
end