function [mov, ifx] = infotaxis_decidir(ifx, pos, dom)
%INFOTAXIS_DECIDIR  Tria el moviment que maximitza el guany esperat
%   d'informacio (eq. 1 de l'article de Vergassola et al., 2007):
%
%   dS(r -> rj) = P(rj) [-S] + [1 - P(rj)] sum_k rho_k(rj) dS_k
%
%   - Primer terme ("explotacio"): probabilitat de trobar la font just a
%     rj; en aquest cas l'entropia passa a ser zero.
%   - Segon terme ("exploracio"): si la font no es a rj, el dron rebra k
%     comptes amb probabilitat rho_k i l'entropia canviara dS_k.
%
%   Es tria, entre quedar-se quiet i els 4 veins, el moviment amb dS mes
%   negatiu (maxima reduccio esperada d'entropia = maxima informacio).
%
%   Per limitar el cost, els recomptes possibles s'agrupen en intervals
%   (ifx.bins = 0, 1, 2, 3-4, 5-7, ...): rho_k i dS_k es calculen per
%   interval. Amb pocs comptes (limit diluit) es exactament l'eq. 1.

P  = ifx.P;
S0 = ifx.S;
b  = ifx.bins;
nB = numel(b);
nMov = size(dom.moviments, 1);
dS = Inf(nMov, 1);

for j = 1:nMov
    q = pos + dom.moviments(j, :);
    if q(1) < 1 || q(2) < 1 || q(1) > dom.Ny || q(2) > dom.Nx
        continue;   % fora del domini
    end

    % Prob. de trobar la font en arribar a rj (dins del radi de deteccio)
    idx = cel_les_properes(q, dom);
    Pj = sum(P(idx));

    % Distribucio condicionada a NO trobar la font a rj
    Pp = P;
    Pp(idx) = 0;
    sp = sum(Pp(:));
    if sp <= 0
        dS(j) = -S0;
        continue;
    end
    Pp = Pp / sp;

    MU = mu_local(ifx.K, q, dom);       % comptes esperats a rj per cada r0

    % Recorrem k = 0, 1, 2, ... acumulant la llei de Poisson per intervals
    canvi = 0;
    acc   = zeros(size(MU));            % P(k < limit actual)
    term  = exp(-MU);                   % P(k) per al k actual
    k = 0;
    for ib = 1:nB
        if ib < nB
            L = zeros(size(MU));
            while k < b(ib + 1)
                L    = L + term;
                term = term .* MU / (k + 1);
                k    = k + 1;
            end
            acc = acc + L;
        else
            L = max(1 - acc, 0);        % interval obert: k >= b(end)
        end
        w   = Pp .* L;
        rho = sum(w(:));                % prob. de rebre comptes en aquest interval
        if rho > 1e-300
            canvi = canvi + rho * (entropia(w / rho) - S0);
        end
    end

    dS(j) = Pj * (-S0) + (1 - Pj) * canvi;
end

ifx.dS = dS;

% Si cap moviment aporta informacio (la posterior ja esta concentrada,
% p. ex. entropia ~ 0 lluny de la font), l'eq. 1 no distingeix entre
% moviments: aleshores el dron vola directament cap al punt mes probable.
if min(dS) > -1e-6
    [~, im] = max(P(:));
    [iy, ix] = ind2sub(size(P), im);
    d = [iy ix] - pos;
    if abs(d(2)) >= abs(d(1))
        mov = [0 sign(d(2))];
    else
        mov = [sign(d(1)) 0];
    end
    ifx.mode = 'cap al maxim';
    return;
end

% Desempat aleatori entre moviments equivalents
dSd = dS + 1e-10 * abs(S0) * rand(size(dS));
[~, jmin] = min(dSd);
mov = dom.moviments(jmin, :);
ifx.mode = 'infotaxis';
end
