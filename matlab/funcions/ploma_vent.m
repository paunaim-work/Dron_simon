function [u, v] = ploma_vent(pl, x, y)
%PLOMA_VENT  Velocitat del vent (m/s) als punts (x, y) per a la ploma
%   complexa: vent uniforme amb meandreig i rafegues + remolins de
%   divergencia nul.la derivats d'una funcio de corrent psi:
%       u = d(psi)/dy,   v = -d(psi)/dx

U  = pl.U * (1 + pl.rafega);
th = pl.theta0 + pl.dTheta;
u  = U * cos(th) * ones(size(x));
v  = U * sin(th) * ones(size(x));

if ~isempty(pl.a)
    arg = (x(:) - pl.desp(1)) * pl.kx' + (y(:) - pl.desp(2)) * pl.ky' + pl.fase';
    cs  = cos(arg);
    u = u + reshape(cs * (pl.a .* pl.ky), size(x));
    v = v - reshape(cs * (pl.a .* pl.kx), size(x));
end
end
