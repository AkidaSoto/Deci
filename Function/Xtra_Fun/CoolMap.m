function cm = CoolMap

black1 = [0 29 63];
magenta = [0,96,200];
Silver = [249,239,255];
Teal = [169,0,0];
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