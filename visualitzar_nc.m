function visualitzar_nc(fitxer)
%VISUALITZAR_NC  Inspecciona i dibuixa una ploma en format .nc.
%
%   Us:
%       visualitzar_nc                          % fa servir el fitxer per defecte (a sota)
%       visualitzar_nc('el_meu_fitxer.nc')
%
%   Que fa:
%     1. Mostra l'estructura del fitxer (dimensions, variables, atributs).
%     2. Calcula estadistiques de la concentracio i on es el maxim.
%     3. Dibuixa talls horitzontals i verticals de la ploma.
%     4. Desa tot el text a "resum_nc.txt" i la figura a "resum_nc.png".
%
%   Despres, enganxa el contingut de resum_nc.txt a la conversa i podre
%   ajustar la simulacio del dron a aquest fitxer.

if nargin < 1 || isempty(fitxer)
    fitxer = 'tgn_port_dispersion_3d_N03_PM10_t0.nc';
end
if ~exist(fitxer, 'file')
    error('No trobo "%s". Ha d''estar a la carpeta actual (%s) o al path.', fitxer, pwd);
end

ruta = which(fitxer);                    % ruta completa si el fitxer es al path
if ~isempty(ruta), fitxer = ruta; end

informe = 'resum_nc.txt';
if exist(informe, 'file'), delete(informe); end
diary(informe);
tanca = onCleanup(@() diary('off'));     % es tanca el diari encara que hi hagi un error

%% ============ 1. ESTRUCTURA DEL FITXER ============
fprintf('==================== FITXER ====================\n');
d = dir(fitxer);
if isempty(d)
    fprintf('Fitxer: %s  (mida desconeguda)\n', fitxer);
else
    fprintf('Fitxer: %s  (%.1f MB)\n', fitxer, d(1).bytes / 1e6);
end
info = ncinfo(fitxer);
fprintf('Format: %s\n', info.Format);

fprintf('\n--- Dimensions ---\n');
for i = 1:numel(info.Dimensions)
    di = info.Dimensions(i);
    fprintf('  %-20s longitud %d%s\n', di.Name, di.Length, ternari(di.Unlimited, '  (il.limitada)', ''));
end

fprintf('\n--- Atributs globals ---\n');
if isempty(info.Attributes)
    fprintf('  (cap)\n');
end
for i = 1:numel(info.Attributes)
    fprintf('  %s = %s\n', info.Attributes(i).Name, attr2str(info.Attributes(i).Value));
end

fprintf('\n--- Variables ---\n');
noms = {info.Variables.Name};
for i = 1:numel(info.Variables)
    v = info.Variables(i);
    dn = strjoin({v.Dimensions.Name}, ', ');
    fprintf('  %-20s [%s]  mida %s  tipus %s\n', v.Name, dn, mat2str(v.Size), v.Datatype);
    for j = 1:numel(v.Attributes)
        fprintf('        @%s = %s\n', v.Attributes(j).Name, attr2str(v.Attributes(j).Value));
    end
end

%% ============ 2. VARIABLE DE CONCENTRACIO ============
nomC = '';
for i = 1:numel(noms)                    % primer, una que digui PM10
    if contains(lower(noms{i}), 'pm10'), nomC = noms{i}; break; end
end
if isempty(nomC)
    for cand = {'c', 'conc', 'concentration', 'C', 'plume', 'tracer', 'concentracio', 'data'}
        if any(strcmp(noms, cand{1})), nomC = cand{1}; break; end
    end
end
if isempty(nomC)                         % si no, la variable mes gran amb >= 2 dimensions
    mida = zeros(1, numel(info.Variables));
    for i = 1:numel(info.Variables)
        if numel(info.Variables(i).Dimensions) >= 2
            mida(i) = prod(info.Variables(i).Size);
        end
    end
    [mx, k] = max(mida);
    if mx == 0, error('No trobo cap variable amb 2 o mes dimensions.'); end
    nomC = noms{k};
end
fprintf('\n>>> Variable de concentracio triada: "%s"\n', nomC);

iv   = find(strcmp(noms, nomC), 1);
dims = info.Variables(iv).Dimensions;
nd   = numel(dims);
c    = double(ncread(fitxer, nomC));     % ncread ja aplica scale_factor, add_offset i _FillValue (-> NaN)

%% ============ 3. ROLS DE LES DIMENSIONS ============
llistaX = {'x', 'lon', 'longitude', 'nx', 'west_east', 'xc', 'easting', 'xi_rho', 'i'};
llistaY = {'y', 'lat', 'latitude', 'ny', 'south_north', 'yc', 'northing', 'eta_rho', 'j'};
llistaZ = {'z', 'alt', 'altitude', 'height', 'lev', 'level', 'nz', 'bottom_top', 'depth', 'k'};
llistaT = {'t', 'time', 'nt', 'times'};
rol = zeros(1, nd);                      % 1=x 2=y 3=z 4=t
for i = 1:nd
    n = lower(dims(i).Name);
    if any(strcmp(llistaX, n)), rol(i) = 1;
    elseif any(strcmp(llistaY, n)), rol(i) = 2;
    elseif any(strcmp(llistaZ, n)), rol(i) = 3;
    elseif any(strcmp(llistaT, n)), rol(i) = 4;
    end
end
if any(rol == 0) || numel(unique(rol)) < nd
    rol = 1:nd;
    fprintf('AVIS: no he reconegut els noms de les dimensions; assumeixo l''ordre (x, y, z, t).\n');
end
fprintf('Dimensions de "%s": ', nomC);
etiq = {'x', 'y', 'z', 't'};
for i = 1:nd
    fprintf('%s=%s(%d) ', dims(i).Name, etiq{rol(i)}, dims(i).Length);
end
fprintf('\n');

ordre = zeros(1, 4);  extra = nd;
for r = 1:4
    k = find(rol == r, 1);
    if isempty(k), extra = extra + 1; ordre(r) = extra; else, ordre(r) = k; end
end
c4 = permute(c, ordre);                  % ara sempre [x y z t]
n4 = size(c4);  n4(end + 1:4) = 1;
nx = n4(1);  ny = n4(2);  nz = n4(3);  nt = n4(4);

% Eixos (vectors); si no n'hi ha, es posen indexos
nomDim = cell(1, 4);
cands  = {llistaX, llistaY, llistaZ, llistaT};
eix    = cell(1, 4);
unit   = {'', '', '', ''};
for r = 1:4
    k = find(rol == r, 1);
    nomDim{r} = ternari(isempty(k), etiq{r}, dims(k).Name);
    [eix{r}, unit{r}, fet] = llegirEix(fitxer, info, [{nomDim{r}}, cands{r}], n4(r));
    if ~fet
        fprintf('AVIS: l''eix %s no es un vector 1D del fitxer (o no existeix): faig servir indexos 1..%d.\n', etiq{r}, n4(r));
    end
end
x = eix{1};  y = eix{2};  z = eix{3};  t = eix{4};

%% ============ 4. ESTADISTIQUES ============
fprintf('\n==================== CONCENTRACIO ====================\n');
fprintf('Mida [nx ny nz nt] = [%d %d %d %d]\n', nx, ny, nz, nt);
fprintf('x: %g .. %g %s (pas %s)\n', x(1), x(end), unit{1}, pasEix(x));
fprintf('y: %g .. %g %s (pas %s)\n', y(1), y(end), unit{2}, pasEix(y));
fprintf('z: %g .. %g %s (pas %s), %d nivell(s)\n', z(1), z(end), unit{3}, pasEix(z), nz);
if nz <= 40, fprintf('   nivells z = %s\n', mat2str(z, 6)); end
fprintf('t: %g .. %g %s (pas %s), %d instant(s)\n', t(1), t(end), unit{4}, pasEix(t), nt);

fprintf('\nTotal valors: %d | NaN: %d (%.1f%%) | zeros: %d (%.1f%%) | negatius: %d\n', ...
        numel(c4), sum(isnan(c4(:))), 100 * mean(isnan(c4(:))), ...
        sum(c4(:) == 0), 100 * mean(c4(:) == 0), sum(c4(:) < 0));
[mx, im] = max(c4(:));
fprintf('Minim: %g | Maxim: %g | Mitjana: %g\n', min(c4(:)), mx, mean(c4(:), 'omitnan'));

idx  = round(linspace(1, numel(c4), min(numel(c4), 2e6)));
vals = c4(idx);  vals = sort(vals(isfinite(vals) & vals > 0));
if ~isempty(vals)
    pc = @(q) vals(max(1, round(q / 100 * numel(vals))));
    fprintf('Percentils dels valors > 0: p50=%g  p90=%g  p99=%g  p99.9=%g\n', pc(50), pc(90), pc(99), pc(99.9));
    fprintf('Valors > 0 mes petits: %g\n', vals(1));
    for llindar = [1e-6 1e-4 1e-2 1 10]
        fprintf('  cel.les amb c > %-8g : %.3f%%\n', llindar, 100 * mean(c4(:) > llindar));
    end
else
    fprintf('AVIS: no hi ha cap valor positiu!\n');
end

[ix, iy, iz, it] = ind2sub(n4, im);
fprintf('\n--- On es el maxim ---\n');
fprintf('Indexs (x,y,z,t) = (%d, %d, %d, %d)\n', ix, iy, iz, it);
fprintf('Posicio: x = %g, y = %g, z = %g %s, t = %g\n', x(ix), y(iy), z(iz), unit{3}, t(it));
fprintf('Valor: %g\n', mx);

if nz > 1
    fprintf('\nPerfil vertical a la columna del maxim (z -> c):\n');
    for k = 1:min(nz, 30)
        fprintf('   z = %-8g  c = %g\n', z(k), c4(ix, iy, k, it));
    end
end

% Centre de massa horitzontal del primer instant (on es concentra la ploma)
c1 = max(c4(:, :, :, 1), 0);  c1(isnan(c1)) = 0;
if sum(c1(:)) > 0
    [X, Y, Z] = ndgrid(x, y, z);
    fprintf('\nCentre de massa (instant 1): x = %g, y = %g, z = %g\n', ...
            sum(X(:) .* c1(:)) / sum(c1(:)), sum(Y(:) .* c1(:)) / sum(c1(:)), sum(Z(:) .* c1(:)) / sum(c1(:)));
end

if nt > 1
    fprintf('\nMaxim per instant (primers 20):\n');
    for k = 1:min(nt, 20)
        tmp = c4(:, :, :, k);
        [m, i2] = max(tmp(:));
        [a, b, cc] = ind2sub([nx ny nz], i2);
        fprintf('   t = %-10g max = %-12g a (x,y,z) = (%g, %g, %g)\n', t(k), m, x(a), y(b), z(cc));
    end
end

%% ============ 5. ALTRES VARIABLES (resum) ============
fprintf('\n==================== ALTRES VARIABLES ====================\n');
for i = 1:numel(info.Variables)
    v = info.Variables(i);
    if strcmp(v.Name, nomC), continue; end
    if prod(v.Size) > 5e7 || ~isnumericType(v.Datatype)
        fprintf('  %-20s (no es resumeix)\n', v.Name);
        continue;
    end
    try
        w = double(ncread(fitxer, v.Name));
        if numel(w) <= 6
            fprintf('  %-20s = %s\n', v.Name, mat2str(w(:)', 8));
        else
            fprintf('  %-20s min %g  max %g  mitjana %g  (mida %s)\n', v.Name, ...
                    min(w(:)), max(w(:)), mean(w(:), 'omitnan'), mat2str(size(w)));
        end
    catch
        fprintf('  %-20s (no s''ha pogut llegir)\n', v.Name);
    end
end

%% ============ 6. DIBUIX ============
fig = figure('Color', 'w', 'Name', ['Ploma: ' fitxer], 'Units', 'normalized', 'Position', [0.03 0.05 0.94 0.85]);
colormap(flipud(hot));
terra = max(mx * 1e-6, vals_min(vals));     % escala logaritmica: 6 ordres de magnitud
cl = log10([terra mx]);
if ~isfinite(cl(1)) || ~isfinite(cl(2)) || cl(1) >= cl(2), cl = [-6 0]; end
L = @(a) log10(max(a, terra));

subplot(2, 3, 1);
imagesc(x, y, L(c4(:, :, iz, it))');  set(gca, 'YDir', 'normal');  caxis(cl);  colorbar;  hold on;
plot(x(ix), y(iy), 'p', 'MarkerSize', 14, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
xlabel(nomDim{1});  ylabel(nomDim{2});
title(sprintf('Tall horitzontal a z = %g (log10 c)', z(iz)));

subplot(2, 3, 2);
imagesc(x, y, L(max(c4(:, :, :, it), [], 3))');  set(gca, 'YDir', 'normal');  caxis(cl);  colorbar;  hold on;
plot(x(ix), y(iy), 'p', 'MarkerSize', 14, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
xlabel(nomDim{1});  ylabel(nomDim{2});
title('Maxim en vertical (vista des de dalt)');

if nz > 1
    subplot(2, 3, 3);
    imagesc(x, z, L(squeeze(c4(:, iy, :, it)))');  set(gca, 'YDir', 'normal');  caxis(cl);  colorbar;  hold on;
    plot(x(ix), z(iz), 'p', 'MarkerSize', 14, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
    xlabel(nomDim{1});  ylabel(nomDim{3});
    title(sprintf('Tall vertical a y = %g', y(iy)));

    subplot(2, 3, 4);
    imagesc(y, z, L(squeeze(c4(ix, :, :, it)))');  set(gca, 'YDir', 'normal');  caxis(cl);  colorbar;  hold on;
    plot(y(iy), z(iz), 'p', 'MarkerSize', 14, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
    xlabel(nomDim{2});  ylabel(nomDim{3});
    title(sprintf('Tall vertical a x = %g', x(ix)));
else
    subplot(2, 3, 3);
    text(0.1, 0.5, 'Nomes hi ha un nivell z (2D)', 'FontSize', 12);  axis off;
end

subplot(2, 3, 5);
if nt > 1
    mxt = squeeze(max(max(max(c4, [], 1), [], 2), [], 3));
    mtt = squeeze(mean(mean(mean(c4, 1, 'omitnan'), 2, 'omitnan'), 3, 'omitnan'));
    semilogy(t, max(mxt, terra), 'r-', t, max(mtt, terra), 'b-');  grid on;
    legend('maxim', 'mitjana', 'Location', 'best');
    xlabel(nomDim{4});  title('Evolucio temporal');
else
    text(0.1, 0.5, 'Nomes hi ha un instant de temps', 'FontSize', 12);  axis off;
end

subplot(2, 3, 6);
if ~isempty(vals)
    histogram(log10(vals), 50);  grid on;
    xlabel('log10(c)  (nomes valors > 0)');  ylabel('cel.les');
    title('Distribucio de valors');
end

try
    saveas(fig, 'resum_nc.png');
catch
end

%% ============ 7. RECONSTRUCCIO 3D (isosuperficies) ============
if nz > 1 && mx > 0
    fig3 = figure('Color', 'w', 'Name', ['Ploma 3D: ' fitxer], 'Units', 'normalized', 'Position', [0.10 0.10 0.80 0.80]);
    % isosurface espera V(y, x, z), per aixo es permuta la matriu
    V = permute(c4(:, :, :, it), [2 1 3]);
    nivells = mx * 10 .^ (-1:-1:-5);                 % 10%, 1%, 0.1%, 0.01% i 0.001% del maxim
    nivells = nivells(nivells > terra);
    cm = flipud(hot(numel(nivells) + 2));
    hold on;
    for k = numel(nivells):-1:1                      % primer el nivell exterior (el mes tenue)
        fv = isosurface(x, y, z, V, nivells(k));
        if ~isempty(fv.vertices)
            patch(fv, 'FaceColor', cm(k + 1, :), 'EdgeColor', 'none', ...
                  'FaceAlpha', 0.12 + 0.55 * (numel(nivells) - k) / max(numel(nivells) - 1, 1));
        end
    end
    plot3(x(ix), y(iy), z(iz), 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
    plot3(x(ix) * [1 1], y(iy) * [1 1], [z(1) z(iz)], 'g-', 'LineWidth', 1.5);
    grid on;  box on;  axis tight;  view(-35, 25);
    camlight;  lighting gouraud;
    xlabel(nomDim{1});  ylabel(nomDim{2});  zlabel(nomDim{3});
    daspect([1 1 0.25]);                             % alcada exagerada x4 per veure-la millor
    title(sprintf('Ploma 3D a t = %g: isosuperficies a %s', t(it), mat2str(nivells, 3)));
    try
        saveas(fig3, 'resum_nc_3d.png');
    catch
    end
end

fprintf('\n==================== FET ====================\n');
fprintf('Text desat a resum_nc.txt i figura a resum_nc.png (carpeta %s).\n', pwd);
fprintf('Enganxa el contingut de resum_nc.txt a la conversa.\n');
end


%% ============================== AUXILIARS ==============================

function s = attr2str(v)
%ATTR2STR  Atribut (text o numeric) com a text curt.
if ischar(v) || isstring(v)
    s = char(strjoin(string(v), ' '));
    if numel(s) > 200, s = [s(1:200) '...']; end
elseif isnumeric(v) || islogical(v)
    s = mat2str(double(v(1:min(numel(v), 10))), 8);
    if numel(v) > 10, s = [s ' ...']; end
else
    s = ['<' class(v) '>'];
end
end

function r = ternari(cond, a, b)
if cond, r = a; else, r = b; end
end

function ok = isnumericType(tipus)
ok = any(strcmp(tipus, {'double', 'single', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64'}));
end

function [v, unitat, fet] = llegirEix(fitxer, info, candidats, n)
%LLEGIREIX  Primer vector 1D de longitud n amb algun dels noms candidats; si no, indexos.
noms = {info.Variables.Name};
v = 1:n;  unitat = '';  fet = false;
for i = 1:numel(candidats)
    k = find(strcmpi(noms, candidats{i}), 1);
    if ~isempty(k) && isnumericType(info.Variables(k).Datatype)
        tmp = ncread(fitxer, noms{k});
        if isvector(tmp) && numel(tmp) == n
            v = double(tmp(:))';
            fet = true;
            for j = 1:numel(info.Variables(k).Attributes)
                if strcmpi(info.Variables(k).Attributes(j).Name, 'units')
                    unitat = attr2str(info.Variables(k).Attributes(j).Value);
                end
            end
            return;
        end
    end
end
end

function s = pasEix(v)
%PASEIX  Pas entre punts d'un eix (o "irregular" si no es constant).
if numel(v) < 2, s = '-'; return; end
d = diff(v);
if max(abs(d - d(1))) <= 1e-6 * max(abs(d(1)), 1)
    s = sprintf('%g', d(1));
else
    s = sprintf('irregular, de %g a %g', min(d), max(d));
end
end

function m = vals_min(vals)
if isempty(vals), m = 1e-12; else, m = vals(1); end
end