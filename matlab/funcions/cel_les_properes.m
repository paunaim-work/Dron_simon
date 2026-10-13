function idx = cel_les_properes(q, dom)
%CEL_LES_PROPERES  Index lineal de les cel.les a distancia <= radi de
%   deteccio de la cel.la q = [iy ix] (la zona que el dron "veu").
c = q + dom.offDet;
c = c(c(:, 1) >= 1 & c(:, 1) <= dom.Ny & c(:, 2) >= 1 & c(:, 2) <= dom.Nx, :);
idx = sub2ind([dom.Ny dom.Nx], c(:, 1), c(:, 2));
end
