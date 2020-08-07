function cm = CoolMap

black1 = [0 29 63];
magenta = [255,0,255];
Silver = [0,255,0];
Teal = [0,0,255];
black2 = [55 0 0];

for i = 1:3
h1(:,i) = linspace(black1(i),magenta(i),64);
end

for i = 1:3
h2(:,i) = linspace(magenta(i),Silver(i),64);
end

for i = 1:3
h3(:,i) = linspace(Silver(i),Teal(i),64);
end

for i = 1:3
h4(:,i) = linspace(Teal(i),black2(i),64);
end

cm = cat(1,h1,h2,h3,h4);
cm = cm/255;