function camp = Palm2ploma(fitxerEntrada, fitxerSortida, posFuga, concMax, zona, ventCap)
%PALM2PLOMA  Converteix un .nc de PALM (variable kc_PM10) al format que llegeix Ploma.llegir.
%
%   1) Mirar el fitxer (nomes escriu informacio):
%        palm2ploma('tgn_port_dispersion_3d_N03_PM10_t0.nc');
%
%   2) Convertir:
%        palm2ploma(entrada, 'ploma_complexa_3d.nc', [xF yF zF], 300, [xmin xmax ymin ymax zmin zmax], '+y');
%
%   fitxerEntrada  .nc de PALM
%   fitxerSortida  .nc que es crea (el que llegira main.m)
%   posFuga        [x y z] de la fuga, en les coordenades ORIGINALS del fitxer PALM (m)
%   concMax        concentracio maxima del camp un cop reescalat (per defecte 300,
%                  com la ploma gaussiana de Ploma.m). El sensor del dron compta
%                  1 detecció/s per unitat de concentracio, i en kg/m3 el senyal seria nul.
%   zona           [xmin xmax ymin ymax zmin zmax] per retallar (coordenades originals).
%                  Opcional, pero molt recomanable: el fitxer sencer son 750x750x202 cel.les.
%
%   ventCap        cap a on bufa el vent a PALM: '+x' (per defecte), '+y', '-x' o '-y'.
%                  main.m/Dron.m suposen vent cap a +x, aixi que el camp es GIRA (rotacio de
%                  90 o 180 graus, sense emmirallar) perque el vent quedi cap a +x.
%                  posFuga es dona sempre en les coordenades originals; la funcio la gira tambe.
%
%   El resultat te x i y desplacades perque comencin a 0 (main.m ho necessita); la z no es toca.
%   Les concentracions negatives o buides (_FillValue) es posen a 0.

nomC = 'kc_PM10';
x = ncread(fitxerEntrada, 'x')';
y = ncread(fitxerEntrada, 'y')';
z = ncread(fitxerEntrada, 'zu_3d')';

% ---- Nomes informacio ----
if nargin < 2
    fprintf('x: %g .. %g (pas %g m), %d punts\n', x(1), x(end), mean(diff(x)), numel(x));
    fprintf('y: %g .. %g (pas %g m), %d punts\n', y(1), y(end), mean(diff(y)), numel(y));
    fprintf('z: %g .. %g m, %d nivells (pas %g .. %g m)\n', z(1), z(end), numel(z), min(diff(z)), max(diff(z)));
    c = ncread(fitxerEntrada, nomC);
    c(~isfinite(c) | c < 0) = 0;
    [cmax, idx] = max(c(:));
    [ix, iy, iz] = ind2sub(size(c), idx);
    fprintf('Concentracio maxima: %g kg/m3 a x=%g, y=%g, z=%g\n', cmax, x(ix), y(iy), z(iz));
    fprintf('(el maxim normalment es a prop de la fuga, pero comprova-ho amb el que saps de la simulacio)\n');
    cx = sum(sum(sum(c, 2), 3)' .* x) / sum(c(:));
    cy = sum(sum(sum(c, 1), 3) .* y) / sum(c(:));
    fprintf('Centre de masses de la ploma: x=%g, y=%g  (la fuga hauria d''estar cap al costat d''on ve el vent)\n', cx, cy);
    dx = cx - x(ix);  dy = cy - y(iy);
    if abs(dx) >= abs(dy), dir = sprintf('%+gx', sign(dx)); else, dir = sprintf('%+gy', sign(dy)); end
    dir = strrep(strrep(dir, '+1', '+'), '-1', '-');
    fprintf('Del maxim al centre de masses: dx=%g, dy=%g  ->  el vent sembla bufar cap a %s\n', dx, dy, dir);    camp = [];
    return
end

if nargin < 4 || isempty(concMax), concMax = 300; end
if nargin < 6 || isempty(ventCap), ventCap = '+x'; end
if nargin < 5 || isempty(zona)
    zona = [x(1) x(end) y(1) y(end) z(1) z(end)];
end

% ---- Retall ----
ix = find(x >= zona(1) & x <= zona(2));
iy = find(y >= zona(3) & y <= zona(4));
iz = find(z >= zona(5) & z <= zona(6));
c = ncread(fitxerEntrada, nomC, [ix(1) iy(1) iz(1) 1], [numel(ix) numel(iy) numel(iz) 1]);
c = squeeze(c);
if ndims(c) < 3, c = reshape(c, numel(ix), numel(iy), numel(iz)); end

% ---- Neteja i escala ----
c = double(c);
c(~isfinite(c) | c < 0) = 0;
cmax = max(c(:));
if cmax <= 0, error('El camp retallat no te cap concentracio positiva.'); end
c = c * (concMax / cmax);
fprintf('Escala aplicada: x %g (maxim original %g kg/m3 -> %g)\n', concMax / cmax, cmax, concMax);

% ---- Estructura que espera Ploma.guardar / Ploma.llegir ----
x0 = x(ix(1));  y0 = y(iy(1));
xl = x(ix) - x0;                       % coordenades locals (comencen a 0)
yl = y(iy) - y0;
fx = posFuga(1) - x0;  fy = posFuga(2) - y0;
Lx = xl(end);  Ly = yl(end);
switch ventCap
    case '+x'                          % no cal girar
        camp.x = xl;  camp.y = yl;  fuga = [fx fy];
    case '+y'                          % rotacio -90: (x,y) -> (y, Lx - x)
        c = flip(permute(c, [2 1 3]), 2);
        camp.x = yl;  camp.y = fliplr(Lx - xl);  fuga = [fy, Lx - fx];
    case '-y'                          % rotacio +90: (x,y) -> (Ly - y, x)
        c = flip(permute(c, [2 1 3]), 1);
        camp.x = fliplr(Ly - yl);  camp.y = xl;  fuga = [Ly - fy, fx];
    case '-x'                          % rotacio 180: (x,y) -> (Lx - x, Ly - y)
        c = flip(flip(c, 1), 2);
        camp.x = fliplr(Lx - xl);  camp.y = fliplr(Ly - yl);  fuga = [Lx - fx, Ly - fy];
    otherwise
        error('ventCap ha de ser ''+x'', ''+y'', ''-x'' o ''-y''.');
end
camp.z = z(iz);
camp.t = 0;
camp.posFuga = [fuga(1), fuga(2), posFuga(3)];
camp.c = single(c);
camp.c(camp.c < 1e-4) = 0;

Ploma.guardar(camp, fitxerSortida);

fprintf('Guardat %s\n', fitxerSortida);
fprintf('  zona: midaX = %g, midaY = %g, z = %g .. %g m\n', camp.x(end), camp.y(end), camp.z(1), camp.z(end));
fprintf('  posFuga dins del nou fitxer: [%g %g %g]\n', camp.posFuga);
fprintf('A main.m posa escenari.midaX/midaY/midaZ a aquests valors i posInicialDron dins de la zona.\n');
end