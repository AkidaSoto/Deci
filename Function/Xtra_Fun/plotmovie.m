pow = squeeze(mean(AvgData{1}.powspctrm,[1 3],'omitnan'));
powdata = AvgData{1};
powdata.dimord = 'chan_time';
powdata = rmfield(powdata,'freq');
powdata.avg = pow;

figure;
clear F
for toi = 1:8:size(powdata.avg,2)
    cfg.xlim = [powdata.time(toi) powdata.time(toi)];
    cfg.colormap = 'jet';
    ft_topoplotER(cfg,powdata);
    drawnow
    if toi == 1
       F(1) =  getframe(gcf);
%        gif('theta.gif')
    else
        F(end+1) = getframe(gcf);
%         gif
    end
end

fig = figure;
movie(fig,F,12);
