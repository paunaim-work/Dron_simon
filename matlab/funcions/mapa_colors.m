function cmap = mapa_colors(tipus, n)
%MAPA_COLORS  Mapes de colors de la visualitzacio.
%   'ploma'        : blanc -> groc -> taronja -> vermell -> granat
%   'probabilitat' : blanc -> blau clar -> blau fosc
if nargin < 2
    n = 256;
end
switch tipus
    case 'ploma'
        ancores = [1 1 1; 1 0.96 0.65; 1 0.72 0.25; 0.88 0.25 0.10; 0.45 0 0.10];
    case 'probabilitat'
        ancores = [1 1 1; 0.78 0.87 1; 0.40 0.58 0.95; 0.15 0.20 0.65; 0.05 0.02 0.30];
end
xa = linspace(0, 1, size(ancores, 1));
xi = linspace(0, 1, n)';
cmap = [interp1(xa, ancores(:, 1), xi), interp1(xa, ancores(:, 2), xi), interp1(xa, ancores(:, 3), xi)];
end
