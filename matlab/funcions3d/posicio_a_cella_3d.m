function c = posicio_a_cella_3d(p, dom)
%POSICIO_A_CELLA_3D  Converteix [x y z] (m) a la cel.la [iy ix iz].
ix = min(max(round(p(1) / dom.dx) + 1, 1), dom.Nx);
iy = min(max(round(p(2) / dom.dx) + 1, 1), dom.Ny);
iz = min(max(round(p(3) / dom.dz) + 1, 1), dom.Nz);
c = [iy ix iz];
end
