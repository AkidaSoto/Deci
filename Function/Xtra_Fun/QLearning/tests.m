
load('C:\Users\Researcher\Documents\GitHub\Nona\Experiments\Multimodal_Exp1_Data\Multimodal_EEG_009BU.mat')
%load('C:\Users\Researcher\Documents\GitHub\Nona\Experiments\Multimodal_Exp1_Data\Mulitmodal_EEG_058CG.mat')
events = StandardizeEventMarkers(events);

Start = find(ismember({events.value},num2str(102)));
End = find(ismember({events.value},num2str(104)));

Block = {events(Start(1):End(1)).value};

TrialStart = find(ismember(Block,num2str(9)));
TrialEnd = find(ismember(Block,num2str(10)));

s = [];
a = [];
r = [];

for j =  1:length(TrialStart)
    
    ismember(strsplit(num2str(20:23)),Block(TrialStart(j):TrialEnd(j)));
    ismember(strsplit(num2str(60:63)),Block(TrialStart(j):TrialEnd(j)));
    ismember(strsplit(num2str([50 52])),Block(TrialStart(j):TrialEnd(j)));
    
    if any(ismember(strsplit(num2str(20:23)),Block(TrialStart(j):TrialEnd(j)))) && any(ismember(strsplit(num2str(60:63)),Block(TrialStart(j):TrialEnd(j)))) && any(ismember(strsplit(num2str([50 52])),Block(TrialStart(j):TrialEnd(j))))
        
        
        Trial = Block(TrialStart(j):TrialEnd(j));
        s(end+1) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str(20:23)))});
        
        a(end+1) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str(60:63)))});
        
        r(end+1) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str([50 52])))});
        
    elseif any(ismember(strsplit(num2str(20:23)),Block(TrialStart(j):TrialEnd(j)))) && any(ismember(strsplit(num2str(64:67)),Block(TrialStart(j):TrialEnd(j)))) && any(ismember(strsplit(num2str([50 52])),Block(TrialStart(j):TrialEnd(j))))
        
        Trial = Block(TrialStart(j):TrialEnd(j));
        s(end+1) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str(20:23)))});
        a(end+1) = nan;
        r(end+1) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str([50 52])))});
        
    end
    
    ismember(strsplit(num2str(16:19)),Block(TrialStart(j):TrialEnd(j)));
    ismember(strsplit(num2str(64:67)),Block(TrialStart(j):TrialEnd(j)));
    ismember(strsplit(num2str([50 52])),Block(TrialStart(j):TrialEnd(j)));
    
    if any(ismember(strsplit(num2str(16:19)),Block(TrialStart(j):TrialEnd(j)))) && any(ismember(strsplit(num2str(64:67)),Block(TrialStart(j):TrialEnd(j)))) && any(ismember(strsplit(num2str([50 52])),Block(TrialStart(j):TrialEnd(j))))
        
        
        Trial = Block(TrialStart(j):TrialEnd(j));
        s(end+1) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str(16:19)))});
        
        a(end+1) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str(64:67)))});
        
        r(end+1) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str([50 52])))});
    elseif any(ismember(strsplit(num2str(16:19)),Block(TrialStart(j):TrialEnd(j)))) && any(ismember(strsplit(num2str(60:63)),Block(TrialStart(j):TrialEnd(j))))&& any(ismember(strsplit(num2str([50 52])),Block(TrialStart(j):TrialEnd(j))))
        
        Trial = Block(TrialStart(j):TrialEnd(j));
        s(end+1) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str(16:19)))});
        a(end+1) = nan;
        r(end+1) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str([50 52])))});
        
    end
    
    
    if ~any(ismember(strsplit(num2str(60:67)),Block(TrialStart(j):TrialEnd(j))))
        Trial = Block(TrialStart(j):TrialEnd(j));
        s(end+1) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str(16:23)))});
        a(end+1) = nan;
        r(end+1) = nan;
    end
    
end

rew = [];

for i = 16:23
    rew(i,:) = r(find(s == i)) == 50;
end
rewd = rew(16:23,:);

r = double(r == 50);

rewd = rew(16:23,:);

%r = [r - 51]*-1;

Fit1.LB = [0 1e-6];
Fit1.UB = [1 30];



Fit2.LB = [0 -1 1e-6];
Fit2.UB = [1 1 30];




for init = 1:10
    
    q = zeros(length(unique(s)),length(unique(a(~isnan(a)))));
%     q = [rand(length(unique(s)),length(unique(a(~isnan(a)))))/2];
%     q =  ones(length(unique(s)),length(unique(a(~isnan(a)))));

Fit1.init =rand(1,length(Fit1.LB)).*(Fit1.UB-Fit1.LB)+Fit1.LB;
Fit2.init =rand(1,length(Fit2.LB)).*(Fit2.UB-Fit2.LB)+Fit2.LB;

[res1,lik1] = ...
    fmincon(@(x) SimpleQ(q,s,a,r,x(1),x(2)),...
    [Fit1.init],[],[],[],[],[Fit1.LB],[Fit1.UB],[],...
    optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));

[res2,lik2] = ...
    fmincon(@(x) ComplexQa(q,s,a,r,x(1),x(2),x(3)),...
    [Fit2.init],[],[],[],[],[Fit2.LB],[Fit2.UB],[],...
    optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));

[res3,lik3] = ...
    fmincon(@(x) ComplexQb(q,s,a,r,x(1),x(2)),...
    [Fit1.init],[],[],[],[],[Fit1.LB],[Fit1.UB],[],...
    optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));



        Result.Eta(init,1) = res1(1);
        Result.Beta(init,1) = res1(2);
        Result.Lik(init,1) = lik1;

        Result.Eta1(init,2) = res2(1);
        Result.Eta2(init,2) = res2(2);
        Result.Beta(init,2) = res2(3);
        Result.Lik(init,2) = lik2;

            Result.Eta(init,3) = res3(1);
            Result.Beta(init,3) = res3(2);
            Result.Lik(init,3) = lik3;

    Result.Lik(init,:)  % to view progress so far
    
end


figure;

plot(mean(rewd,1),'LineWidth',4);

for l = 1:size(Result.Lik,2)
BestFit(l) = find(Result.Lik(:,l) == min(Result.Lik(:,l)),1,'first');

if l == 1
[~,qs,P] = SimpleQ(q,s,a,r, Result.Eta(BestFit(l),l), Result.Beta(BestFit(l),l));
elseif l == 2
[~,qs,P] = ComplexQa(q,s,a,r, Result.Eta1(BestFit(l),l),Result.Eta2(BestFit(l),l), Result.Beta(BestFit(l),l));
elseif l == 3
[~,qs,P] = ComplexQb(q,s,a,r, Result.Eta(BestFit(l),l), Result.Beta(BestFit(l),l));
end

Pb = Simu(P);

trew = [];

for t = 1:1
    rew = [];
    for i = 16:23
        rew(i,:) = Pb(find(s == i));
    end
    rew = rew(16:23,:);
    
    trew(:,:,t) = rew;
end
hold on


plot(mean(mean(trew,3),1),'LineWidth',2)

end

legend(['Actual Data'],...
       ['SimpleQ (' num2str(Result.Lik(BestFit(1),1)) ')' ], ...
       ['AQ(' num2str(Result.Lik(BestFit(2),2)) ')' ],...
       ['BQ (' num2str(Result.Lik(BestFit(3),3)) ')' ]);



function [LLE,qout,P] = ComplexQa(q,s,a,r,alp,alp2,beta)

%[qout,PE,P,LLE] =

lik = 0;

s = s - min(s) + 1;
a = a - min(a) + 1;

prsa = repmat({[]},[length(unique(s)) length((unique(a(~isnan(a)))))]);
APE  = repmat({[0]},[length(unique(s)) length((unique(a(~isnan(a)))))]);

betasa = repmat(beta,[length(unique(s)) length((unique(a(~isnan(a)))))]);


for Stim = 1:length(s)
    
    if ~isnan(a(Stim))
        for Act = unique(a(~isnan(a)))
            
            if a(Stim) == Act
                
                PE(s(Stim), Act) = (r(Stim) - q(s(Stim),Act));
                q(s(Stim), Act) = q(s(Stim), Act) +  alp*PE(s(Stim),a(Stim));
                
                prsa{s(Stim), Act}(end+1) = r(Stim);
                APE{s(Stim), Act}(end+1) = abs(nanmean(prsa{s(Stim), Act}) - q(s(Stim), Act));
                
                dAPE = APE{s(Stim), Act}(end-1) - APE{s(Stim), Act}(end);
                
%                 if dAPE > 0
%                    betasa(s(Stim), Act) =  betasa(s(Stim), Act)*2;
%                 else 
%                    betasa(s(Stim), Act) =  betasa(s(Stim), Act)/2;
%                 end

                if dAPE > 0
                   alp =  alp^[1 - alp2];
                else 
                   alp =  alp^[1 + alp2];
                end
                
%                 if dAPE > 0
%                     betasa(s(Stim), Act) =  betasa(s(Stim), Act)*2;
%                 else
%                     betasa(s(Stim), Act) =  betasa(s(Stim), Act)/2;
%                 end
%                 
                
                P(Stim) = exp(betasa(s(Stim), Act)*q(s(Stim),Act))/sum(exp(betasa(s(Stim), Act)*q(s(Stim),:)));
                
                if isnan(P(Stim))
                   P(Stim) = 1;
                end
                
                qout(Stim) = q(s(Stim), Act);
            else
                PE(s(Stim), Act) = 0;
            end
            
            
             
            
        end
    else
        
        P(Stim) = nan;
    end
end

LLE = -nansum(log(P),2);

end

function [LLE,qout,P] = SimpleQ(q,s,a,r,alp,beta)

%[qout,PE,P,LLE] =

lik = 0;

for Stim = 1:length(s)
    
    if ~isnan(a(Stim))
        for Act = unique(a(~isnan(a)))
            
            if a(Stim) == Act
                
                PE(s(Stim)-15, Act-59) = (r(Stim) - q(s(Stim)-15,Act-59));
                q(s(Stim)-15, Act-59) = q(s(Stim)-15, Act-59) + alp*PE(s(Stim)-15,a(Stim)-59);
                
                P(Stim) = exp(beta*q(s(Stim)-15,Act-59))/sum(exp(beta*q(s(Stim)-15,:)));
                qout(Stim) = q(s(Stim)-15, Act-59);
            else
                PE(s(Stim)-15, Act-59) = 0;
            end
            
            
             
            
        end
    else
        
        P(Stim) = nan;
    end
end

LLE = -nansum(log(P),2);

end


function [LLE,qout,P] = ComplexQb(q,s,a,r,alp,beta)

%[qout,PE,P,LLE] =

lik = 0;

s = s - min(s) + 1;
a = a - min(a) + 1;

prsa = repmat({[]},[length(unique(s)) length((unique(a(~isnan(a)))))]);
APE  = repmat({[0]},[length(unique(s)) length((unique(a(~isnan(a)))))]);

betasa = repmat(beta,[length(unique(s)) length((unique(a(~isnan(a)))))]);

alphasa = repmat(alp,[length(unique(s)) length((unique(a(~isnan(a)))))]);

for Stim = 1:length(s)
    
    if ~isnan(a(Stim))
        for Act = unique(a(~isnan(a)))
            
            if a(Stim) == Act
                
                PE(s(Stim), Act) = (r(Stim) - q(s(Stim),Act));
                q(s(Stim), Act) = q(s(Stim), Act) + alp*PE(s(Stim),a(Stim));
                
                prsa{s(Stim), Act}(end+1) = r(Stim);
                APE{s(Stim), Act}(end+1) = abs(nanmean(prsa{s(Stim), Act}) - q(s(Stim), Act));
                
                dAPE = APE{s(Stim), Act}(end-1) - APE{s(Stim), Act}(end);
                
                if dAPE > 0
                   betasa(s(Stim), Act) =  betasa(s(Stim), Act)*2;
                else 
                   betasa(s(Stim), Act) =  betasa(s(Stim), Act)/2;
                end
                
                P(Stim) = exp(betasa(s(Stim), Act)*q(s(Stim),Act))/sum(exp(betasa(s(Stim), Act)*q(s(Stim),:)));
                
                if isnan(P(Stim))
                   P(Stim) = 1;
                end
                
                qout(Stim) = q(s(Stim), Act);
            else
                PE(s(Stim), Act) = 0;
            end
            
            
             
            
        end
    else
        
        P(Stim) = nan;
    end
end

LLE = -nansum(log(P),2);

end



function B = Simu(A)

for a = 1:length(A)
    if ~isnan(A(a))
      
       B(a) = rand < A(a);
        
    else
       B(a) = 0;
    end
end

end

