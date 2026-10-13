function gr = gradient_actualitzar_3d(gr, pos, k, dom)
%GRADIENT_ACTUALITZAR_3D  Estima el gradient 3D per diferencies finites
%   centrades amb els 7 sensors k = [C E O N S amunt avall] i el suavitza
%   amb un filtre exponencial.

c = max(k(:)' / dom.dt - gr.fons, 0) / gr.sensibilitat;
gNou = [c(2) - c(3), c(4) - c(5), c(6) - c(7)] / (2 * gr.d);

if sum(k) >= gr.llindar
    gr.passosSenseSenyal = 0;
    gr.g = gr.alfa * gr.g + (1 - gr.alfa) * gNou;
else
    gr.passosSenseSenyal = gr.passosSenseSenyal + 1;
    if gr.passosSenseSenyal > gr.memoria
        gr.g = [0 0 0];
    end
end

gr.mapaSuma(pos(1), pos(2), pos(3)) = gr.mapaSuma(pos(1), pos(2), pos(3)) + c(1);
gr.mapaN(pos(1), pos(2), pos(3))    = gr.mapaN(pos(1), pos(2), pos(3)) + 1;
end
