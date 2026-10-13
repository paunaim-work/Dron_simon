function c = ploma_gaussiana_conc_3d(ddx, ddy, z, H, p)
%PLOMA_GAUSSIANA_CONC_3D  Ploma gaussiana estacionaria en 3D.
%
%   c = ploma_gaussiana_conc_3d(ddx, ddy, z, H, p) dona la concentracio a
%   l'alcada z en punts desplacats horitzontalment (ddx, ddy) de la font,
%   si la font es a l'alcada H. z i H poden ser escalars o matrius de la
%   mateixa mida que ddx. p conte Q, u, estabilitat, sigma0 i theta.
%
%            Q            n^2              (z-H)^2            (z+H)^2
%   c = ------------ exp(- ------) [ exp(- -------) + exp(- -------) ]
%       2 pi u sy sz      2 sy^2           2 sz^2             2 sz^2
%
%   El segon terme vertical es la reflexio del contaminant al terra.

s =  ddx * cos(p.theta) + ddy * sin(p.theta);   % a favor del vent
n = -ddx * sin(p.theta) + ddy * cos(p.theta);   % transversal
z = z + zeros(size(s));
H = H + zeros(size(s));

c = zeros(size(s));
m = s > -p.sigma0;
if ~any(m(:))
    return;
end
sm = max(s(m), 0);
[sy, sz] = sigma_briggs(sm, p.estabilitat);
sy = sqrt(sy.^2 + p.sigma0^2);
sz = sqrt(sz.^2 + p.sigma0^2);

Hm = H(m);
if isfield(p, 'flot') && ~isempty(p.flot)
    % Densitat del gas (mode realista): la ploma puja o baixa amb la distancia
    Hm = max(Hm + desplacament_flotacio(sm / p.u, p.flot), 0);
end
vertical = exp(-(z(m) - Hm).^2 ./ (2 * sz.^2)) + exp(-(z(m) + Hm).^2 ./ (2 * sz.^2));
c(m) = p.Q ./ (2 * pi * p.u * sy .* sz) .* exp(-n(m).^2 ./ (2 * sy.^2)) .* vertical;
end
