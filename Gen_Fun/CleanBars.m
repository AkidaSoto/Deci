function CleanBars(data,varargin)

A =gca;
if length(find(size(data) ~= 1)) == 1
    
    
    bar(A,[zeros([1 size(data,2)]);data]);
    arrayfun(@(c) set(c,'YData',c.YData(2)),A.Children)
    
else
    bar(A,data);
end

hold on

if ~isempty(varargin{1})
    ngroups = size(data, 1);
    nbars = size(data, 2);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, data(:,i), varargin{1}(:,i), 'k', 'linestyle', 'none');
    end
end

end