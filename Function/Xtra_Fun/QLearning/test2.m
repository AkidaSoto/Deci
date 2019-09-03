

load('C:\Users\Researcher\Documents\GitHub\Nona\Experiments\ProbStim_Data\ProbStim_Stim_024ST');

events = StandardizeEventMarkers(events);

Start = find(ismember({events.value},strsplit(num2str([3]))));
End = find(ismember({events.value},strsplit(num2str([7]))));



s = [];
a = [];
r = [];

rew = [];   

for b = 1:length(Start)
    
    Block = {events(Start(b):End(b)).value};
    
    
    TrialStart = find(ismember(Block,num2str(5)));
    TrialEnd = find(ismember(Block,num2str(6)));

    
    for j =  1:length(TrialStart)
        
        Trial = Block(TrialStart(j):TrialEnd(j));
        stim = ismember(strsplit(num2str(14:16)),Block(TrialStart(j):TrialEnd(j)));
        format = ismember(strsplit(num2str(41:42)),Block(TrialStart(j):TrialEnd(j)));
        side = ismember(strsplit(num2str(37:38)),Block(TrialStart(j):TrialEnd(j)));
        rewq = ismember(strsplit(num2str([31 34])),Block(TrialStart(j):TrialEnd(j)));
        
        chooseA = [];
        if [format(1) && side(1)] || [format(2) && side(2)]
            chooseA = 3 + find(stim);
        else
            chooseA = find(stim);
        end
        
        if  any(stim) && any(rewq) && ~isempty(chooseA)
            s(b,j) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str(14:16)))});
            
            a(b,j) = chooseA;
            
            r(b,j) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str([31 34])))});
        else 
            s(b,j) = str2num(Trial{ismember(Block(TrialStart(j):TrialEnd(j)),strsplit(num2str(14:16)))});
            a(b,j) = nan;
            
            r(b,j) = nan;
            
        end
        
        
        
    end
    
 
    
for i = 14:16
    rew(b,i-13,:) = r(b,s(b,:) == i) == 34;
    acts(b,i-13,:) = a(b,s(b,:) == i) == i+3-13;
    sti(b,s(b,:) == i) = i;
    racts(b,i-13,:) = arrayfun(@(c) mean(acts(b,i-13,1:c)),[1:size(acts,3)]);
end




mean(acts(b,1,:))
rewd = squeeze(mean(rew,1));
end

pacts = squeeze(mean(racts,1))';

plot(pacts);
legend({'AB' 'CD' 'EF'});

alp = .2;
beta = .5;

for set = 1:size(a,1)
        q = squeeze(zeros(size(acts(set,:,:))))';
        
        SimpleQ(q,s(set,:),a(set,:),r(set,:),alp,beta);

end

function [LLE,qout,P] = SimpleQ(q,s,a,r,alp,beta)

s = s - min(s) + 1;
a = a - min(a) + 1;
r = r - min(r) + 1;

lik = 0;

for Stim = 1:length(s)
    
    if ~isnan(a(Stim))

            if a(Stim) -3 < 1
                
                PE(s(Stim), Act) = (r(Stim) - q(s(Stim),a(Stim)));
                q(s(Stim), Act) = q(s(Stim), Act) + alp*PE(s(Stim),a(Stim));
                
                P(Stim) = exp(beta*q(s(Stim),Act))/sum(exp(beta*q(s(Stim),:)));
                qout(Stim) = q(s(Stim), Act);
            end

    else
        
        P(Stim) = nan;
    end
end

LLE = -nansum(log(P),2);

end




