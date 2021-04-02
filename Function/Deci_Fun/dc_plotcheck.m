
function [Deci,info] = dc_plotcheck(Deci,Dims)

for Dim = 1:length(Dims)
    if Deci.Plot.(Dims{Dim}).do
        if isequal(Deci.Plot.(Dims{Dim}).Channel,'Reinhart-All')
            Deci.Plot.(Dims{Dim}).Channel = dc_getchans('noeyes');
        end
        
        if isequal(Deci.Plot.(Dims{Dim}).Channel,'Reinhart-All_eyes')
            Deci.Plot.(Dims{Dim}).Channel = dc_getchans('all');
        end
        
        if isequal(Deci.Plot.(Dims{Dim}).Channel,'RNET')
            Deci.Plot.(Dims{Dim}).Channel = dc_getchans('RNET');
        end
        
        if isequal(Deci.Plot.(Dims{Dim}).Channel,'PD-All')
            Deci.Plot.(Dims{Dim}).Channel = dc_getchans('PD-All');
        end
        
        Tois{Dim} = Deci.Plot.(Dims{Dim}).Toi;
        Fois{Dim} = Deci.Plot.(Dims{Dim}).Foi;
        Chois{Dim} = Deci.Plot.(Dims{Dim}).Channel;
    end
end

Tois = sort([Tois{:}]);
Tois = [Tois(1) Tois(end)];

Fois = sort([Fois{:}]);
Fois = [Fois(1) Fois(end)];

Chois = unique([Chois{:}]);

if length(Deci.Plot.Topo.Channel) == 1 && Deci.Plot.Topo.do
    error('You asked for Topo but requested one channel!')
end

info.Chois = Chois;
info.Fois = Fois;
info.Tois = Tois;

end