
index = ismember(layout.label,Deci.Plot.Freq.Channel)
layout.pos = layout.pos(index,:);
layout.width = layout.width(index);
layout.height = layout.height(index);
layout.label = layout.label(index);

save('C:\Users\User\Documents\GitHub\Deci\Gen_Fun\Collab_eye.mat','layout');

index = ~ismember(layout.label,{'VEM','HEM'});
layout.pos = layout.pos(index,:);
layout.width = layout.width(index);
layout.height = layout.height(index);
layout.label = layout.label(index);

save('C:\Users\User\Documents\GitHub\Deci\Gen_Fun\Collab_noeye.mat','layout');