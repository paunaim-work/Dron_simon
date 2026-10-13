function idx = cel_les_properes_3d(q, dom)
%CEL_LES_PROPERES_3D  Index lineal de les cel.les de la zona de deteccio
%   (cilindre) al voltant de la cel.la q = [iy ix iz]: les que el dron
%   "veu" des d'alla. Vegeu crear_domini_3d.m.
c = q + dom.offDet;
c = c(c(:, 1) >= 1 & c(:, 1) <= dom.Ny & c(:, 2) >= 1 & c(:, 2) <= dom.Nx ...
      & c(:, 3) >= 1 & c(:, 3) <= dom.Nz, :);
idx = sub2ind([dom.Ny dom.Nx dom.Nz], c(:, 1), c(:, 2), c(:, 3));
end
