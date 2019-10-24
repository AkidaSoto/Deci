function cm = SixMap

black = [18 2 26];
magenta = [112,18,101];
pink = [175 53 125];
orange = [209 112 71];
gold = [205 172 76];
white = [255 255 255];




for i = 1:3
h1(:,i) = linspace(black(i),magenta(i),51);
end

for i = 1:3
h2(:,i) = linspace(magenta(i),pink(i),51);
end

for i = 1:3
h3(:,i) = linspace(pink(i),orange(i),51);
end

for i = 1:3
h4(:,i) = linspace(orange(i),gold(i),51);
end

for i = 1:3
h5(:,i) = linspace(gold(i),white(i),52);
end

cm = cat(1,h1,h2,h3,h4,h5);
cm = cm/255;