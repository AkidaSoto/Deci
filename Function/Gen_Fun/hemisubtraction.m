function dataout = hemisubtraction(cfg,datain)

cfg = Exist(cfg,'parameter','powspctrm');


labels = datain.label;

%odd
elecs_odd = cell2mat(cellfun(@(c) any(rem(str2num(c(isstrprop(c,'digit'))),2)),  labels,'un',0));
elecs_even = cell2mat(cellfun(@(c) any(~rem(str2num(c(isstrprop(c,'digit'))),2)),  labels,'un',0));


dataout = datain;


dataout.(cfg.parameter) = dataout.(cfg.parameter)(elecs_even,:,:) - dataout.(cfg.parameter)(elecs_odd,:,:);

end