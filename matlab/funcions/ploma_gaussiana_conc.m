function c = ploma_gaussiana_conc(ddx, ddy, p)
%PLOMA_GAUSSIANA_CONC  Concentracio de la ploma gaussiana estacionaria.
%
%   c = ploma_gaussiana_conc(ddx, ddy, p) avalua la concentracio en punts
%   desplacats (ddx, ddy) (m) respecte de la font. p conte: Q, u, H, z,
%   estabilitat, sigma0 i theta (direccio cap on bufa el vent, rad).
%
%            Q            n^2              (z-H)^2            (z+H)^2
%   c = ------------ exp(- ------) [ exp(- -------) + exp(- -------) ]
%       2 pi u sy sz      2 sy^2           2 sz^2             2 sz^2
%
%   on s es la distancia a sotavent i n la distancia transversal. El
%   segon terme vertical es la "reflexio" al terra. A sobrevent (s < 0)
%   la concentracio es zero. Per evitar la singularitat a s = 0, la font
%   te una mida inicial sigma0: sy^2 = sy(s)^2 + sigma0^2 (idem sz).

s =  ddx * cos(p.theta) + ddy * sin(p.theta);   % a favor del vent
n = -ddx * sin(p.theta) + ddy * cos(p.theta);   % transversal

c = zeros(size(s));
m = s > -p.sigma0;                   % inclou el punt de la font
if ~any(m(:))
    return;
end
[sy, sz] = sigma_briggs(max(s(m), 0), p.estabilitat);
sy = sqrt(sy.^2 + p.sigma0^2);
sz = sqrt(sz.^2 + p.sigma0^2);

vertical = exp(-(p.z - p.H)^2 ./ (2 * sz.^2)) + exp(-(p.z + p.H)^2 ./ (2 * sz.^2));
c(m) = p.Q ./ (2 * pi * p.u * sy .* sz) .* exp(-n(m).^2 ./ (2 * sy.^2)) .* vertical;
end
