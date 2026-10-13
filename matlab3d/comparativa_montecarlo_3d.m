%% COMPARATIVA MONTE CARLO DELS ALGORISMES DE CERCA (3D)
%
%  Repeteix moltes cerques 3D (sense animacio) per a cada combinacio de
%  ploma i algorisme i en calcula la taxa d'exit i el temps de cerca. Si
%  escenarisAleatoris = true, cada repeticio te un escenari diferent (el
%  mateix per a tots els algorismes, perque la comparacio sigui justa).
%  L'infotaxis 3D es mes costos de calcular que el 2D: compta uns
%  quants minuts.

clear; close all; clc;
carpeta = fileparts(mfilename('fullpath'));
if isempty(carpeta)
    carpeta = pwd;
end
addpath(fullfile(carpeta, 'funcions3d'));
addpath(fullfile(carpeta, '..', 'matlab', 'funcions'));

%% ========================= OPCIONS =========================
nRepeticions     = 10;                         % cerques per combinacio
plomes           = {'gaussiana', 'complexa'};
algorismes       = {'infotaxis', 'gradient'};
escenarisAleatoris = true;  % false: sempre l'escenari fix de configuracio_defecte_3d.m
%% ===========================================================

nP = numel(plomes);
nA = numel(algorismes);
taxaExit = zeros(nP, nA);
temps = cell(nP, nA);

for ip = 1:nP
    for ia = 1:nA
        tt = nan(nRepeticions, 1);
        for r = 1:nRepeticions
            cfg = configuracio_defecte_3d();
            cfg.tipusPloma = plomes{ip};
            cfg.algorisme  = algorismes{ia};
            cfg.visual.mostrar = false;
            cfg.llavor = 1000 * ip + r;
            cfg.atzar.activat = escenarisAleatoris;
            cfg.atzar.llavor  = r;          % escenari r (igual per a tots els algorismes)
            cfg.atzar.mostrar = false;
            res = simular_cerca_3d(cfg);
            if res.trobada
                tt(r) = res.temps;
            end
            fprintf('%-9s %-9s  %2d/%2d  trobada = %d  t = %5.0f s\n', plomes{ip}, ...
                    algorismes{ia}, r, nRepeticions, res.trobada, res.temps);
        end
        taxaExit(ip, ia) = mean(~isnan(tt));
        temps{ip, ia} = tt;
    end
end

%% Taula de resultats
fprintf('\n%-10s %-10s %8s %12s %12s\n', 'Ploma', 'Algorisme', 'Exit(%)', 'T mitja (s)', 'T mediana(s)');
for ip = 1:nP
    for ia = 1:nA
        tt = temps{ip, ia};
        tt = tt(~isnan(tt));
        if isempty(tt)
            tm = NaN; tmed = NaN;
        else
            tm = mean(tt); tmed = median(tt);
        end
        fprintf('%-10s %-10s %8.0f %12.0f %12.0f\n', plomes{ip}, algorismes{ia}, ...
                100 * taxaExit(ip, ia), tm, tmed);
    end
end

%% Grafiques
figure('Name', 'Comparativa Monte Carlo 3D', 'Color', 'w', 'Position', [100 100 1100 420]);
subplot(1, 2, 1);
bar(100 * taxaExit);
set(gca, 'XTick', 1:nP, 'XTickLabel', plomes);
ylabel('Cerques amb exit (%)');
legend(algorismes, 'Location', 'southoutside', 'Orientation', 'horizontal');
ylim([0 105]);
grid on;
title('Taxa d''exit (3D)');

subplot(1, 2, 2);
hold on;
colors = [0.10 0.35 0.85; 0.85 0.30 0.10];
etiq = {};
x = 0;
for ip = 1:nP
    for ia = 1:nA
        x = x + 1;
        tt = temps{ip, ia};
        tt = tt(~isnan(tt));
        if ~isempty(tt)
            plot(x + 0.15 * (rand(size(tt)) - 0.5), tt, 'o', 'Color', colors(ia, :));
            plot(x + [-0.3 0.3], median(tt) * [1 1], '-', 'Color', colors(ia, :), 'LineWidth', 3);
        end
        etiq{end + 1} = sprintf('%s-%s', plomes{ip}(1:5), algorismes{ia}(1:4)); %#ok<SAGROW>
    end
end
set(gca, 'XTick', 1:x, 'XTickLabel', etiq);
xlim([0.5 x + 0.5]);
ylabel('Temps de cerca (s)');
title('Temps de les cerques amb exit, 3D (linia = mediana)');
grid on;
