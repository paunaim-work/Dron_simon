function T = Montecarlo(nNormal, nLent)
%MONTECARLO  Monte Carlo: crida main.m amb fugues i dron en posicions aleatories.
%   Posa'l al costat de main.m, Dron.m i Ploma.m i executa:  T = montecarlo(100, 30);
%     nNormal  simulacions per escenari (per defecte 100)
%     nLent    simulacions per als escenaris lents: 3D amb infotaxis (per defecte 30)
%   Els escenaris lents fan servir les N primeres posicions dels que no ho son.
%   Per a cada escenari (dimensio x ploma x algorisme) fa N simulacions i treu
%   el % d'exit (amb interval de confiança del 95 %) i dos temps mitjans:
%     tempsTrobades  mitja nomes de les simulacions que troben la fuga
%     tempsTotes     mitja de totes; les que fallen compten com l'autonomia sencera
if nargin < 1, nNormal = 100; end
if nargin < 2, nLent = 30; end
Nmax = max(nNormal, nLent);            % posicions sortejades per a tots els escenaris

midaX = 400;  midaY = 300;            % han de coincidir amb l'escenari de main.m
distMinima = 20;                      % distancia minima dron-fuga al sortir (m)

% ---- Posicions aleatories (les mateixes per a tots els escenaris) ----
rng(2024);
xFuga = midaX * rand(Nmax, 1);
yFuga = midaY * rand(Nmax, 1);
zFuga = 5 + 35 * rand(Nmax, 1);                       % dins de l'alcada de vol del dron (5-40 m)
xDron = xFuga + (midaX - xFuga) .* rand(Nmax, 1);     % vent avall de la fuga
yDron = midaY * rand(Nmax, 1);
zDron = 5 + 35 * rand(Nmax, 1);
for i = 1:Nmax                                        % re-sorteja si surt massa a prop
    while norm([xDron(i) - xFuga(i), yDron(i) - yFuga(i), zDron(i) - zFuga(i)]) < distMinima
        xDron(i) = xFuga(i) + (midaX - xFuga(i)) * rand;
        yDron(i) = midaY * rand;
        zDron(i) = 5 + 35 * rand;
    end
end

dimensions = {'2d', '3d'};
ploms      = {'gaussiana', 'complexa'};
algorismes = {'gradient', 'infotaxis'};

files = {};
for dimensio = dimensions
for tipus = ploms
for algorisme = algorismes
    opcions.dimensio   = dimensio{1};
    opcions.tipusPloma = tipus{1};
    opcions.algorisme  = algorisme{1};
    opcions.animacio   = false;
    opcions.silenci    = true;
    if strcmp(opcions.dimensio, '3d') && strcmp(opcions.algorisme, 'infotaxis')
        N = nLent;                    % escenari lent
    else
        N = nNormal;
    end

    trobada = false(N, 1);
    temps   = nan(N, 1);
    for i = 1:N
        opcions.llavor = 1000 + i;
        res = main([xFuga(i) yFuga(i) zFuga(i)], [xDron(i) yDron(i) zDron(i)], opcions);
        trobada(i) = res.trobada;
        temps(i)   = res.temps;
    end

    % Interval de Wilson al 95 % per al % d'exit
    k = sum(trobada);  p = k / N;  z = 1.96;
    centre = (p + z^2 / (2 * N)) / (1 + z^2 / N);
    ample  = z * sqrt(p * (1 - p) / N + z^2 / (4 * N^2)) / (1 + z^2 / N);
    pctExit = 100 * p;  icMin = 100 * (centre - ample);  icMax = 100 * (centre + ample);

    tempsTrobades = mean(temps(trobada));      % NaN si cap simulacio la troba
    tempsTotes    = mean(temps);               % els fracassos ja valen l'autonomia (res.temps)
    fprintf('%s | %-9s | %-9s | exit %5.1f %% [%4.1f-%5.1f] | temps trobades %5.0f s | temps totes %5.0f s\n', ...
            upper(opcions.dimensio), opcions.tipusPloma, opcions.algorisme, ...
            pctExit, icMin, icMax, tempsTrobades, tempsTotes);
    files(end + 1, :) = {opcions.dimensio, opcions.tipusPloma, opcions.algorisme, N, ...
                         pctExit, icMin, icMax, tempsTrobades, tempsTotes, ...
                         median(temps(trobada))}; %#ok<AGROW>
end
end
end

T = cell2table(files, 'VariableNames', ...
    {'dimensio', 'ploma', 'algorisme', 'N', 'exitPercent', 'icMin', 'icMax', ...
     'tempsTrobades', 'tempsTotes', 'tempsMedianaTrobades'});
writetable(T, 'montecarlo_resultats.csv');
end