Dirs = CleanDir('C:\Users\User\Desktop\Shrey\SRT\ProcessedData\Artifact');

for dir = 1:length(Dirs)
   load(['C:\Users\User\Desktop\Shrey\SRT\ProcessedData\Artifact' filesep Dirs{dir}])
   condinfo = data.condinfo;
   preart = data.preart;
   
   if size([condinfo{1} condinfo{2}],1) ~= length(unique([condinfo{1} condinfo{2}],'rows'))
      error('not unique') 
   end
   
   for trl = 1:size(condinfo{1},1)
       trlnum(trl) = find(ismember([preart{1} preart{2}(:,1:end-1)],[condinfo{1}(trl,:) condinfo{2}(trl,1:end-1)],'rows'));
   end
   condinfo{1} = [condinfo{1} condinfo{3}];
   condinfo{3} = trlnum'; 
   
   preart{1} = [preart{1} preart{3}];
   preart{3} = [1:length(preart{3})]';

   data.condinfo = condinfo;
   data.preart = preart;
   
   
   save(['C:\Users\User\Desktop\Shrey\SRT\ProcessedData\Artifact' filesep Dirs{dir}],'data','-v7.3');
end

