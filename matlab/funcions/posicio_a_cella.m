function c = posicio_a_cella(p, dom)
%POSICIO_A_CELLA  Converteix una posicio [x y] (m) a la cel.la [iy ix].
ix = min(max(round(p(1) / dom.dx) + 1, 1), dom.Nx);
iy = min(max(round(p(2) / dom.dx) + 1, 1), dom.Ny);
c = [iy ix];
end
