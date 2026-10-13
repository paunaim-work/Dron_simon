function [mov, ifx] = infotaxis_decidir_3d(ifx, pos, dom)
%INFOTAXIS_DECIDIR_3D  Tria, entre quedar-se quiet i els 6 veins (+-x,
%   +-y, amunt, avall), el moviment que maximitza el guany esperat
%   d'informacio (eq. 1 de Vergassola et al., 2007):
%
%   dS(r -> rj) = P(rj) [-S] + [1 - P(rj)] sum_k rho_k(rj) dS_k
%
%   Vegeu infotaxis_decidir.m (versio 2D) per a l'explicacio detallada.

P  = ifx.P;
S0 = ifx.S;
b  = ifx.bins;
nB = numel(b);
nMov = size(dom.moviments, 1);
dS = Inf(nMov, 1);

for j = 1:nMov
    q = pos + dom.moviments(j, :);
    if any(q < 1) || q(1) > dom.Ny || q(2) > dom.Nx || q(3) > dom.Nz
        continue;
    end

    idx = cel_les_properes_3d(q, dom);
    Pj = sum(P(idx));                   % prob. de trobar la font en arribar a rj

    Pp = P;                             % distribucio si NO la troba a rj
    Pp(idx) = 0;
    sp = sum(Pp(:));
    if sp <= 0
        dS(j) = -S0;
        continue;
    end
    Pp = Pp / sp;

    MU = mu_local_3d(ifx.K, q, dom);

    canvi = 0;
    acc   = zeros(size(MU));
    term  = exp(-MU);
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
            L = max(1 - acc, 0);
        end
        w   = Pp .* L;
        rho = sum(w(:));
        if rho > 1e-12
            % entropia de w/rho calculada sense crear la matriu normalitzada
            wn = w(w > 0);
            Sk = log(rho) - sum(wn .* log(wn)) / rho;
            canvi = canvi + rho * (Sk - S0);
        end
    end

    dS(j) = Pj * (-S0) + (1 - Pj) * canvi;
end
ifx.dS = dS;

% Si cap moviment aporta informacio, vola cap al punt mes probable
if min(dS) > -1e-6
    [~, im] = max(P(:));
    [iy, ix, iz] = ind2sub(size(P), im);
    d = [iy ix iz] - pos;
    [~, eix] = max(abs(d));
    mov = [0 0 0];
    mov(eix) = sign(d(eix));
    ifx.mode = 'cap al maxim';
    return;
end

dSd = dS + 1e-10 * abs(S0) * rand(size(dS));
[~, jmin] = min(dSd);
mov = dom.moviments(jmin, :);
ifx.mode = 'infotaxis';
end
