function out =  cooljet

out = rgb2hsv(jet);

out(:,2) = out(:,2)*3;
out(:,3) = out(:,3)*2;

out(out(:,2) > 1,2) = 1;
out(out(:,3) > 1,3) = 1;

out = hsv2rgb(out);

end