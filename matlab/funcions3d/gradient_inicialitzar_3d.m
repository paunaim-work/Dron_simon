function gr = gradient_inicialitzar_3d(cfg, dom)
%GRADIENT_INICIALITZAR_3D  Cercador per gradient en 3D. El dron porta 7
%   sensors: un al centre, quatre en creu als bracos (+-x, +-y) i dos en
%   vertical (un a sobre i un penjat a sota), a una distancia d.

g = cfg.gradient;
if isempty(g.separacio)
    gr.d = dom.dx;
else
    gr.d = g.separacio;
end
d = gr.d;
gr.offsets = [0 0 0; d 0 0; -d 0 0; 0 d 0; 0 -d 0; 0 0 d; 0 0 -d];  % C, E, O, N, S, amunt, avall

gr.alfa         = g.alfa;
gr.memoria      = g.memoria;
gr.persistencia = g.persistencia;
gr.probVertical = g.probVertical;
gr.sensibilitat = cfg.sensor.sensibilitat;
gr.fons         = cfg.sensor.fons;

% Llindar de deteccio (falsa alarma pel fons < 1 %)
lamFons = size(gr.offsets, 1) * cfg.sensor.fons * dom.tMesura;
gr.llindar = 1;
while 1 - sum(exp(-lamFons) * lamFons.^(0:gr.llindar - 1) ./ factorial(0:gr.llindar - 1)) > 0.01
    gr.llindar = gr.llindar + 1;
end

gr.g = [0 0 0];                       % gradient filtrat [dc/dx dc/dy dc/dz]
gr.passosSenseSenyal = Inf;
gr.mode = 'exploracio';
gr.dir  = randi([2 5]);               % comenca movent-se en horitzontal
gr.invers = [1 3 2 5 4 7 6];

gr.mapaSuma = zeros(dom.Ny, dom.Nx, dom.Nz);
gr.mapaN    = zeros(dom.Ny, dom.Nx, dom.Nz);
end
