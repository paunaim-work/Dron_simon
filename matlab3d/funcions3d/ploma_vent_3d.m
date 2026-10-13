function [u, v] = ploma_vent_3d(pl, x, y, z)
%PLOMA_VENT_3D  Vent horitzontal (m/s) als punts (x, y, z) de la ploma
%   complexa. El vent mitja creix amb l'alcada (perfil potencial
%   U(z) = U (z/zRef)^p, tipic de la capa limit atmosferica), meandreja i
%   te rafegues; a sobre s'hi afegeixen remolins horitzontals de
%   divergencia nul.la (u = d(psi)/dy, v = -d(psi)/dx).

perfil = (max(z, 1) / pl.zRef).^pl.pPerfil;
U  = pl.U * (1 + pl.rafega) * perfil;
th = pl.theta0 + pl.dTheta;
u  = U * cos(th);
v  = U * sin(th);

if ~isempty(pl.a)
    arg = (x(:) - pl.desp(1)) * pl.kx' + (y(:) - pl.desp(2)) * pl.ky' + pl.fase';
    cs  = cos(arg);
    u = u + reshape(cs * (pl.a .* pl.ky), size(x));
    v = v - reshape(cs * (pl.a .* pl.kx), size(x));
end
end
