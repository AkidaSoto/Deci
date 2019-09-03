function out = QL(Deci,info,freq,params)
%params.States = [20 21 23 24];
%params.Actions = [31 32];
%params.Reward  = [51 52];
%params.Value  = {[20 10] [10 0] [0 -10] [-10 -20]};

trialtypes = freq.condinfo{2}(info.alltrials,:);


%Basic Data Maintence for finding block numbers
if any(any(trialtypes < 0))
    blocknumbers = sort(unique(trialtypes(:,end)),'descend');
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
    Fit1.UB = [1 3];
    
    Fit1.init =rand(1,length(Fit1.LB)).*(Fit1.UB-Fit1.LB)+Fit1.LB;
    
    Fit2.LB = [0 0 1e-6];
    Fit2.UB = [1 1 3];
    Fit2.init =rand(1,length(Fit2.LB)).*(Fit2.UB-Fit2.LB)+Fit2.LB;
    
    Fit3.LB = [0 0 0 0 1e-6];
    Fit3.UB = [1 1 1 1 3];
    Fit3.init =rand(1,length(Fit3.LB)).*(Fit3.UB-Fit3.LB)+Fit3.LB;
    
    Fit4.LB = [0 0 1e-6 1e-6];
    Fit4.UB = [1 1 3 3];
    Fit4.init =rand(1,length(Fit4.LB)).*(Fit4.UB-Fit4.LB)+Fit4.LB;
    
    q = ones([1 length(unique(params.Actions))]) * params.Start;
    %q = zeros([1 length(unique(Actions{blk}))]);
    %q = rand([1 length(unique(Actions{blk}))]);
    
    % run through each model with the randomized starting values
    for blk = 1:length(blocknumbers)
        
            [Value{1,init,blk}] = ...
                fmincon(@(x) SimpleQ(q,Actions{blk},Rewards{blk},x(1),x(2)),...
                [Fit1.init],[],[],[],[],[Fit1.LB],[Fit1.UB],[],...
                optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
            
            [LLE(1,init,blk),qs,P{1,init,blk}] = SimpleQ(q,Actions{blk},Rewards{blk}, Value{1,init,blk}(1),Value{1,init,blk}(2));
            
            
            [Value{2,init,blk}] = ...
                fmincon(@(x) FeedbackQ(q,Actions{blk},Rewards{blk},x(1),x(2),x(3)),...
                [Fit2.init],[],[],[],[],[Fit2.LB],[Fit2.UB],[],...
                optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
            
            [LLE(2,init,blk),qs,P{2,init,blk}] = FeedbackQ(q,Actions{blk},Rewards{blk}, Value{2,init,blk}(1),Value{2,init,blk}(2),Value{2,init,blk}(3));
            
            [Value{3,init,blk}] = ...
                fmincon(@(x) InteractiveQ(q,Actions{blk},Rewards{blk},x(1),x(2),x(3),x(4),x(5)),...
                [Fit3.init],[],[],[],[],[Fit3.LB],[Fit3.UB],[],...
                optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
            
            [LLE(3,init,blk),qs,P{3,init,blk}] = InteractiveQ(q,Actions{blk},Rewards{blk}, Value{3,init,blk}(1),Value{3,init,blk}(2),Value{3,init,blk}(3),Value{3,init,blk}(4),Value{3,init,blk}(5));
            
%             [Value] = ...
%                 fmincon(@(x) ModelQ(q,Actions{blk},Rewards{blk},x(1),x(2),x(3),x(4)),...
%                 [Fit4.init],[],[],[],[],[Fit4.LB],[Fit4.UB],[],...
%                 optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
%             
%             [LLE(4,init,blk),qs,P] = ModelQ(q,Actions{blk},Rewards{blk}, Value(1), Value(2),Value(3),Value(4))
%             
            LLE2(:,init,blk) = aicbic(-LLE(:,init,blk),[2 3 5]);
    end
end

[Best,I] = min(LLE,[],2);
[Best2,I2] = min(LLE2,[],2);

for blk = 1:size(Best,3)
    for m = 1:size(Best,1)
       %out{1,m,blk} = Value{m,I(m,1,blk),blk};
       out{m,blk} = Value{m,I2(m,1,blk),blk};
    end
end

% 
% for blk = 1:size(Best,3)
%     
%     figure;
%     
%     
%     plot(Actions{blk} == params.Actions(1))
%     
%     Pb = [];
%     for m = 1:size(Best,1)
%         for t = 1:1000
%             Pb(t,:) = Simu(P{m,I(m,1,blk),blk});
%         end
%         hold on
%         plot(mean(Pb,1))
%     end
%     
%     legend(['Actual Data'],...
%     ['SimpleQ (' num2str(Best(1,blk)) ')' ], ...
%     ['FeedbackQ(' num2str(Best(2,blk)) ')' ],...
%     ['InteractiveQ (' num2str(Best(3,blk)) ')' ]);
% 
% title(['block ' num2str(blk)])
% 
% end



    function [LLE,q,P] = SimpleQ(q,a,r,alp,beta)
        
        %[qout,PE,P,LLE] =
        
        lik = 0;
        
        actions = sort(unique(a));
        
        P = [];
        
        for Act = 1:length(a)
            
            which = find(a(Act) == actions);
            
            PE(Act) = (r(Act) -q(Act,which) );
            q(Act+1,which) = q(Act,which) + alp*PE(Act);
            q(Act+1,find(a(Act) ~= actions))  = q(Act,find(a(Act) ~= actions));
            
            P(Act) = exp(beta*q(Act+1,which))/sum(exp(beta*q(Act+1,:)));
            qout(Act) = q(Act+1,which);
            
            
        end
        
        LLE = -nansum(log(P));
        
    end

    function [LLE,q,P] = FeedbackQ(q,a,r,alpP,alpN,beta)
        
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

    function [LLE,q,P] = InteractiveQ(q,a,r,alpPP,alpNP,alpPR,alpNR,beta)
        
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

    function [LLE,q,P] = ModelQ(q,a,r,alpP,alpN,betaP,betaR)
        
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