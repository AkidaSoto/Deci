function cm = CoolMap2

magenta = [140,0,79];
Silver = [238,236,251];
Teal = [0,129,80];

for i = 1:3
fh(:,i) = linspace(magenta(i),Silver(i),128);
end

for i = 1:3
lh(:,i) = linspace(Silver(i),Teal(i),128);
end

cm = cat(1,fh,lh);
cm = cm/255;