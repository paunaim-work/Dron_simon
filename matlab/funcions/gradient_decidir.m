function [mov, gr] = gradient_decidir(gr, pos, dom)
%GRADIENT_DECIDIR  Tria el moviment segons el gradient estimat.
%
%   - Mode SEGUIMENT (hi ha senyal recent): es mou al vei en la direccio
%     de maxima projeccio del gradient (on la concentracio augmenta mes).
%   - Mode EXPLORACIO (no detecta res): passeig aleatori correlacionat;
%     mante la direccio amb probabilitat "persistencia" i, si no, gira.

valid = false(1, 5);
for j = 2:5
    q = pos + dom.moviments(j, :);
    valid(j) = q(1) >= 1 && q(2) >= 1 && q(1) <= dom.Ny && q(2) <= dom.Nx;
end

if gr.passosSenseSenyal <= gr.memoria && any(gr.g ~= 0)
    gr.mode = 'seguiment';
    % projeccio del gradient sobre cada direccio (dom.moviments es [diy dix])
    proj = -Inf(1, 5);
    for j = 2:5
        if valid(j)
            proj(j) = gr.g(1) * dom.moviments(j, 2) + gr.g(2) * dom.moviments(j, 1);
        end
    end
    millors = find(proj == max(proj));
    gr.dir = millors(randi(numel(millors)));
else
    gr.mode = 'exploracio';
    if rand > gr.persistencia || ~valid(gr.dir)
        cand = find(valid);
        cand(cand == gr.invers(gr.dir)) = [];     % evita fer marxa enrere
        if isempty(cand)
            cand = find(valid);
        end
        gr.dir = cand(randi(numel(cand)));
    end
end
mov = dom.moviments(gr.dir, :);
end
