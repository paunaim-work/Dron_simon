function [mov, gr] = gradient_decidir_3d(gr, pos, dom)
%GRADIENT_DECIDIR_3D  Seguiment del gradient 3D o, sense senyal, passeig
%   aleatori correlacionat (sobretot horitzontal; de tant en tant canvia
%   d'alcada amb probabilitat probVertical).

valid = false(1, 7);
for j = 2:7
    q = pos + dom.moviments(j, :);
    valid(j) = all(q >= 1) && q(1) <= dom.Ny && q(2) <= dom.Nx && q(3) <= dom.Nz;
end

if gr.passosSenseSenyal <= gr.memoria && any(gr.g ~= 0)
    gr.mode = 'seguiment';
    proj = -Inf(1, 7);
    for j = 2:7
        if valid(j)
            m = dom.moviments(j, :);      % [diy dix diz]
            proj(j) = gr.g(1) * m(2) + gr.g(2) * m(1) + gr.g(3) * m(3);
        end
    end
    millors = find(proj == max(proj));
    gr.dir = millors(randi(numel(millors)));
else
    gr.mode = 'exploracio';
    if rand > gr.persistencia || ~valid(gr.dir)
        if rand < gr.probVertical
            cand = 6:7;
        else
            cand = 2:5;
        end
        cand = cand(valid(cand));
        cand(cand == gr.invers(gr.dir)) = [];
        if isempty(cand)
            cand = find(valid);
            cand(cand == 1) = [];
        end
        gr.dir = cand(randi(numel(cand)));
    end
    % els moviments verticals d'exploracio son d'un sol pas
    if gr.dir >= 6
        mov = dom.moviments(gr.dir, :);
        gr.dir = randi([2 5]);
        return;
    end
end
mov = dom.moviments(gr.dir, :);
end
