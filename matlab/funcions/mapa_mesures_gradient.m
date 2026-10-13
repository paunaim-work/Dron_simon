function M = mapa_mesures_gradient(gr, vis)
%MAPA_MESURES_GRADIENT  Mapa (log10) de la concentracio mitjana mesurada
%   pel sensor central del dron a cada cel.la visitada.
M = log10(max(gr.mapaSuma ./ max(gr.mapaN, 1), vis.cmin));
M(gr.mapaN == 0) = vis.clim(1);
end
