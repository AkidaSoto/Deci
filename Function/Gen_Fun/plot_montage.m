function plot_montage(table,type)

tits = cellstr(unique(table.OutcomeofInterest));

chans = cellfun(@(c) strsplit(c,', '),cellstr(table.Montage),'UniformOutput',false);
chans = [cellfun(@(c) c(1),chans) cellfun(@(c) c(2),chans)];

for t = 1:length(tits)
    
    ES = table.EffectSize(ismember(table.OutcomeofInterest,tits{t}));
    FR = table.FrequencyHz(ismember(table.OutcomeofInterest,tits{t}));
    FR = cellfun(@num2str,mat2cell(FR,ones([length(FR) 1]),1),'un',0);
    CH = chans(ismember(table.OutcomeofInterest,tits{t}),:);
    

    [B,I,O] = unique([array2table(FR) array2table(CH)],'rows');
    ind = find(~ismember(1:size(CH,1)',I));
    
    if strcmp(type,'mean')
        dd = ismember(1:length(CH)',I);
        x = ind(find(ES(ind) == max(ES(ind))));
        ES(x) = mean(ES(ind));
        
        dd(x) = true;
        CH = CH(dd,:);
        FR = FR(dd);
        ES = ES(dd);
        
    elseif strcmp(type,'max')
        x = ind(find(ES(ind) == max(ES(ind))));
        dd = ismember(1:length(CH)',I);
        dd(x) = true;
        CH = CH(dd,:);
        FR = FR(dd);
        ES = ES(dd);
    end
    
    for f = 1:length(FR) 
        
        fr = str2num(FR{f});
        
        if fr <= 4
            F(f) = 3;
        elseif fr > 4 && fr <= 8 
            F(f) = 6;
        elseif fr > 8 && fr <= 12
            F(f) = 10;
        elseif fr > 12 && fr <= 30
            F(f) = 24;
        elseif fr > 30 && fr <= 80
            F(f) = 55;
        end
        
        
    end
    
    uu = table2array(unique(array2table(CH),'rows'));
    stringuu = strcat(uu(:,1),uu(:,2));
    
    free = [3 6 10 24 55];
    ff = nan(length(uu),1,5);
    
    for ch = 1:size(CH,1)
        stringch = strcat(CH(ch,1),CH(ch,2));
        
        
        ff(ismember(stringuu,stringch),1,ismember(free,F(ch))) = ES(ch);
        
    end
    
    data.powspctrm= ff+2;
    data.labelcmb = uu;
    data.freq = free;
    
    cfg.foi(1) = min(data.freq);
    cfg.foi(2) = max(data.freq);
    cfg.layout = 'easycap_supraorbital.mat';
    cfg.colorparam = 'freq';
    cfg.widthparam = 'powspctrm';
    
    %data.powspctrm = arrayfun(@(c) rand*6,repmat(1,[4 1 2]));
    %data.freq = [1 3];
    %data.method = 'max';
    %data.label = {'FCz' 'Pz' 'CP3' 'F3'};
    %data.labelcmb = {'FCz' 'Pz'; 'FCz' 'CP3'; 'FCz' 'F3'; 'FCz' 'F4'};
    

    
    
    
    %Exist(data,'label');
    
    data.dimord = 'chan_chan_freq';
    %cfg.inputfile = tits{t};
    ft_topoplotCC_freqs(cfg,data);
    title(tits{t})
    
end

end