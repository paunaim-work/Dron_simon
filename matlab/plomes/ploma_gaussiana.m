function C = ploma_gaussiana(x, y, z, Q, u, H)
% PLOMA_GAUSSIANA  Concentració d'una ploma gaussiana estacionària.
%   Font puntual a (0,0,H), vent u en la direcció x, reflexió al terra.
if nargin < 4, Q = 1; end
if nargin < 5, u = 2; end
if nargin < 6, H = 10; end
sy = 0.22 * x ./ sqrt(1 + 0.0001 * x);
sz = 0.20 * x;
C = (Q ./ (2*pi*u.*sy.*sz)) .* exp(-y.^2 ./ (2*sy.^2)) .* ...
    (exp(-(z-H).^2 ./ (2*sz.^2)) + exp(-(z+H).^2 ./ (2*sz.^2)));
C(x <= 0) = 0;
end
