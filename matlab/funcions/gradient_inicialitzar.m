function gr = gradient_inicialitzar(cfg, dom)
%GRADIENT_INICIALITZAR  Estat inicial del cercador per gradient
%   (quimiotaxi). El dron porta 5 sensors: un al centre i quatre en creu
%   a una distancia d (+x, -x, +y, -y), com si n'hi hagues un a cada
%   bra d'un quadricopter.

g = cfg.gradient;
if isempty(g.separacio)
    gr.d = dom.dx;
else
    gr.d = g.separacio;
end
gr.offsets = [0 0; gr.d 0; -gr.d 0; 0 gr.d; 0 -gr.d];   % [x y]: C, E, O, N, S

gr.alfa         = g.alfa;
gr.memoria      = g.memoria;
gr.persistencia = g.persistencia;
gr.sensibilitat = cfg.sensor.sensibilitat;
gr.fons         = cfg.sensor.fons;

% Llindar de deteccio: minim de comptes (sumant els 5 sensors) perque
% la probabilitat de falsa alarma pel fons sigui < 1 %
lamFons = size(gr.offsets, 1) * cfg.sensor.fons * dom.dt;
gr.llindar = 1;
while 1 - sum(exp(-lamFons) * lamFons.^(0:gr.llindar - 1) ./ factorial(0:gr.llindar - 1)) > 0.01
    gr.llindar = gr.llindar + 1;
end

gr.g = [0 0];                         % gradient filtrat [dc/dx dc/dy]
gr.passosSenseSenyal = Inf;
gr.mode = 'exploracio';
gr.dir  = randi([2 5]);               % direccio actual (index de dom.moviments)
gr.invers = [1 3 2 5 4];              % direccio oposada

% Mapa de memoria del dron (nomes per visualitzar el que ha mesurat)
gr.mapaSuma = zeros(dom.Ny, dom.Nx);
gr.mapaN    = zeros(dom.Ny, dom.Nx);
end
