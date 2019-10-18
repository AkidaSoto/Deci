function cm = CoolMap

magenta = [0,96,200];
Silver = [249,239,255];
Teal = [169,0,0];

for i = 1:3
fh(:,i) = linspace(magenta(i),Silver(i),128);
end

for i = 1:3
lh(:,i) = linspace(Silver(i),Teal(i),128);
end

cm = cat(1,fh,lh);
cm = cm/255;