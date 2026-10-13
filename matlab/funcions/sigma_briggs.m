function [sy, sz] = sigma_briggs(s, classe)
%SIGMA_BRIGGS  Coeficients de dispersio sigma_y, sigma_z (m) en funcio de
%   la distancia a sotavent s (m) i la classe d'estabilitat de
%   Pasquill-Gifford (A = molt inestable ... F = molt estable).
%   Formules de Briggs (1973) per a terreny rural.

s = max(s, 0);
switch upper(classe)
    case 'A'
        sy = 0.22 * s .* (1 + 1e-4 * s).^(-1/2);
        sz = 0.20 * s;
    case 'B'
        sy = 0.16 * s .* (1 + 1e-4 * s).^(-1/2);
        sz = 0.12 * s;
    case 'C'
        sy = 0.11 * s .* (1 + 1e-4 * s).^(-1/2);
        sz = 0.08 * s .* (1 + 2e-4 * s).^(-1/2);
    case 'D'
        sy = 0.08 * s .* (1 + 1e-4 * s).^(-1/2);
        sz = 0.06 * s .* (1 + 1.5e-3 * s).^(-1/2);
    case 'E'
        sy = 0.06 * s .* (1 + 1e-4 * s).^(-1/2);
        sz = 0.03 * s .* (1 + 3e-4 * s).^(-1);
    case 'F'
        sy = 0.04 * s .* (1 + 1e-4 * s).^(-1/2);
        sz = 0.016 * s .* (1 + 3e-4 * s).^(-1);
    otherwise
        error('Classe d''estabilitat desconeguda: %s (usa A..F)', classe);
end
end
