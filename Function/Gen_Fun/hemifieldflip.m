function dataout = hemifieldflip(datain,varargin)

if ~isempty(varargin)
    labelname = varargin{1};
else
    labelname = 'label';
end

labels = datain.(labelname);

%odd
elecs_odd = cell2mat(cellfun(@(c) any(rem(str2num(c(isstrprop(c,'digit'))),2)),  labels,'un',0));
elecs_even = cell2mat(cellfun(@(c) any(~rem(str2num(c(isstrprop(c,'digit'))),2)),  labels,'un',0));

labels(elecs_odd) = cellfun(@(c) [c(isstrprop(c,'alpha')) num2str(str2num(c(isstrprop(c,'digit')))+1) ], labels(elecs_odd),'un',0);
labels(elecs_even) = cellfun(@(c) [c(isstrprop(c,'alpha')) num2str(str2num(c(isstrprop(c,'digit')))-1) ], labels(elecs_even),'un',0);

dataout = datain;
dataout.(labelname) = labels;
end