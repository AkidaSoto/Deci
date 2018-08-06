function trl = timeestfun(cfg)

Row{1} = [{'28'} {'17'} {'35'} {'27'}];
Row{2} = [{'37'} {'17'} {'35'} {'27'}];
Row{3} = [{'28'} {'17'} {'35'} {'27'}];
Row{4} = [{'37'} {'17'} {'35'} {'27'}];

Minust{1} = '28';
Minust{2} = '28';
Minust{3} = '37';
Minust{4} = '37';

Lock = '35';
ReLock = '27';

Learned = '33';

startstop = {'10' '11'};
bstartstop = {'5' '6'};
sstime = [-2 3];

event = ft_read_event(cfg.dataset);
hdr   = ft_read_header(cfg.dataset);


trl = [];
retrl = [];

event =  StandardizeEventMarkers(event);

if rem(length(find(ismember({event.value},startstop)')),2) ~= 0
    error('invalid number of start and end pairings, check data!')
end

if rem(length(find(ismember({event.value},bstartstop)')),2) ~= 0
    event = event(1:find(ismember({event.value},'6'),1,'last'));
end

startstopseg = reshape(find(ismember({event.value},startstop)'),[2 length(find(ismember({event.value},startstop)'))/2]);
bstartstopseg = reshape(find(ismember({event.value},bstartstop)'),[2 length(find(ismember({event.value},bstartstop)'))/2]);
Learned = find(ismember({event.value},Learned));

for b = 1:length(Learned)
    
    blocks = find(Learned(b) > bstartstopseg(1,:) & Learned(b) < bstartstopseg(2,:));
    
    trials = find(all(bstartstopseg(1,blocks) < startstopseg & bstartstopseg(2,blocks) >startstopseg,1));
    
    for t = trials(1:end-1)
        value = {event(startstopseg(1,t):startstopseg(2,t)).value};
        sample = {event(startstopseg(1,t):startstopseg(2,t)).sample};
        
        
        for i = 1:length(Row)
            
            condvalue = all(ismember(Row{i},value));
            Minus = 0;
            
            if t ~= 1
                value2 = {event(startstopseg(1,t-1):startstopseg(2,t-1)).value};
                
                condvalue2 = all(ismember(Minust{i},value2));
                
                if condvalue2
                    Minus = 1 ;
                end
                
            end
            
            
            if Minus && condvalue
                condlock = sample{ismember(value,Lock)};
                
                condrelock = sample{ismember(value,ReLock)};
                
                begsample = condlock + sstime(1)*hdr.Fs;
                endsample = condrelock + sstime(2)*hdr.Fs;
                retrl(end +1,:) = [condrelock - condlock]/hdr.Fs;
                offset        = sstime(1)*hdr.Fs;
                
                trl(end + 1,:) = [begsample endsample offset i b];
                
            end
            
        end
        
    end
    
    
end

[~,k] = sort(trl(:,4));
retrl = retrl(k);
mkdir([cfg.Raw filesep 'Redefine']);
save([cfg.Raw filesep 'Redefine' filesep cfg.Subject], 'retrl');
