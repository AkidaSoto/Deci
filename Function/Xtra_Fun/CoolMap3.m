function cm = CoolMap3

magenta = [150,59,0];
Silver = [255,255,244];
Teal = [0,123,9];

for i = 1:3
fh(:,i) = linspace(magenta(i),Silver(i),128);
end

for i = 1:3
lh(:,i) = linspace(Silver(i),Teal(i),128);
end

cm = cat(1,fh,lh);
cm = cm/255;