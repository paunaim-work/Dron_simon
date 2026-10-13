function [bat, cal, percent] = bateria_actualitzar(bat, mov, p, dom)
%BATERIA_ACTUALITZAR  Descompta l'energia d'un pas i comprova si el dron
%   ha d'acabar la missio.
%
%   [bat, cal, percent] = bateria_actualitzar(bat, mov, p, dom)
%     mov     : moviment fet en aquest pas (zeros = s'ha quedat quiet)
%     p       : posicio actual del dron (m)
%     cal     : true si s'ha de tornar a la base: l'energia que queda,
%               menys la d'un pas mes, no arriba per tornar al punt
%               d'enlairament (o, sense reserva de retorn, s'ha esgotat)
%     percent : bateria util que queda (%)
cal = false;
percent = NaN;
if isempty(bat)
    return;
end
if any(mov ~= 0)
    bat.E = bat.E - bat.Pm * dom.tViatge - bat.Ph * dom.tEspera;
else
    bat.E = bat.E - bat.Ph * (dom.tViatge + dom.tEspera);
end
ePas = bat.Pm * dom.tViatge + bat.Ph * dom.tEspera;
eRetorn = 0;
if bat.retorn
    eRetorn = bat.Pm * norm(p - bat.base) / bat.v;
end
cal = bat.E - ePas <= eRetorn;
percent = max(100 * bat.E / bat.E0, 0);
end
