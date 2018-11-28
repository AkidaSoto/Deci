function dataout = hemifieldflip(datain)

labels = datain.label;

%odd
elecs_odd = cell2mat(cellfun(@(c) any(rem(c(isstrprop(c,'digit')),2)),  labels,'un',0));
elecs_even = cell2mat(cellfun(@(c) any(~rem(c(isstrprop(c,'digit')),2)),  labels,'un',0));

labels(elecs_odd) = cellfun(@(c) [c(isstrprop(c,'alpha')) num2str(str2num(c(isstrprop(c,'digit')))+1) ], labels(elecs_odd),'un',0);
labels(elecs_even) = cellfun(@(c) [c(isstrprop(c,'alpha')) num2str(str2num(c(isstrprop(c,'digit')))-1) ], labels(elecs_even),'un',0);

dataout = datain;
dataout.label = labels;
end