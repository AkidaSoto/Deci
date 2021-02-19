function h = CleanBars(data,varargin)

A =gca;
if length(find(size(data) ~= 1)) == 1
    
    
    h = bar(A,[zeros([1 size(data,2)]);data]);
    arrayfun(@(c) set(c,'YData',c.YData(2)),h)
    
else
    h = bar(A,data);
end

hold on

if ~isempty(varargin)
    ngroups = size(data, 1);
    nbars = size(data, 2);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        
        if ~isempty(varargin{1})
        errorbar(x, data(:,i), varargin{1}(:,i), 'k', 'linestyle', 'none');
        end
        
        if length(varargin) >= 2
            if ~isempty(varargin{2})
                randylim = (2*1-1) * groupwidth / (2*nbars);
                rander = rand([length(varargin{2}(:,i)) 1])-.5;
                randy =  rander* [randylim];
                
                scatter(x*ones([size(varargin{2},1) 1])+randy, varargin{2}(:,i),'k')
            end
        end
        
        if length(varargin) >= 3
             if ~isempty(varargin{3})
                 if varargin{3}(i) == 1
                 scatter(x,data(i)*1.1,'k*')
                 end
             end
        end
        
    end
end

end