function [Deci, Subjects] = dc_math(Deci,Subjects,info)   

display(' ')
    display(['Doing ' num2str(length(Deci.Plot.Math)) ' Maths'] )

%     %reduce number of math needed 
%     
%     for figs = find(Deci.Plot.Figures)
%         for draw = 1:length(Deci.Plot.Draw{figs})
%             
%             if Deci.Plot.Draw{figs}(draw) > length(Deci.Analysis.Conditions)
%                 
%                 arginstr = sprintf('x%i,', 1:length(varargin));
%                 
%                 
%             end
%             
%             
%         end
%     end
%     
%     while 
%     
%     
%       
%         arginstr = arginstr(1:end-1); % remove the trailing ','
%         eval(sprintf('operation = @(x) %s;',regexprep( cfg.operation,'x(\d*)','x{$1}')));
%         
%         
        
    for conds = 1:length(Deci.Plot.Math)
        for subj = 1:size(Subjects,1)
            scfg.parameter = info.parameter;
            
            if contains(Deci.Plot.Math{conds},'N')
                pos = find(Deci.Plot.Math{conds} == 'N');
                Deci.Plot.Math{conds} = replace(Deci.Plot.Math{conds},'N',num2str(length(size(Subjects{subj,1}.powspctrm))+1));
                
            end
            
            scfg.operation = Deci.Plot.Math{conds};

%                         if any(cellfun(@isempty,{Subjects{subj,:}}))
%                 empties = cellfun(@isempty,{Subjects{subj,:}});
%                 nempties = ~cellfun(@isempty,{Subjects{subj,:}});
% 
%                 for emp = find(empties)
%                     Subjects{subj,emp} = Subjects{subj,nempties(1)};
%                     Subjects{subj,emp}.(info.parameter)(:) = nan; 
%                 end
%             end

            evalc('MathData{subj} = ft_math(scfg,Subjects{subj,:})');
        end 
        
        Subjects(:,size(Subjects,2)+1) = MathData;
    end
    
end