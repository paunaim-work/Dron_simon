function gr = gradient_actualitzar(gr, pos, k, dom)
%GRADIENT_ACTUALITZAR  Estima el gradient de concentracio a partir dels
%   comptes dels 5 sensors k = [C E O N S] per diferencies finites
%   centrades:
%
%       dc/dx ~ (c_E - c_O) / (2d),    dc/dy ~ (c_N - c_S) / (2d)
%
%   La concentracio s'estima com c = (k/dt - fons) / sensibilitat, i el
%   gradient es suavitza amb un filtre exponencial (mitjana mobil):
%   g <- alfa g + (1 - alfa) g_nou, per reduir el soroll de recompte.

c = max(k(:)' / dom.dt - gr.fons, 0) / gr.sensibilitat;
gNou = [(c(2) - c(3)) / (2 * gr.d), (c(4) - c(5)) / (2 * gr.d)];

if sum(k) >= gr.llindar
    gr.passosSenseSenyal = 0;
    gr.g = gr.alfa * gr.g + (1 - gr.alfa) * gNou;
else
    gr.passosSenseSenyal = gr.passosSenseSenyal + 1;
    if gr.passosSenseSenyal > gr.memoria
        gr.g = [0 0];                 % senyal perdut: oblidem el gradient
    end
end

gr.mapaSuma(pos(1), pos(2)) = gr.mapaSuma(pos(1), pos(2)) + c(1);
gr.mapaN(pos(1), pos(2))    = gr.mapaN(pos(1), pos(2)) + 1;
end
