function vis = actualitzar_isosuperficies(vis, C, dom)
%ACTUALITZAR_ISOSUPERFICIES  Redibuixa les superficies d'igual
%   concentracio de la ploma 3D (una de "vora" i una d'"interior").
for i = 1:numel(vis.hIso)
    f = [];
    if max(C(:)) > vis.nivells(i)
        [f, v] = isosurface(dom.X3, dom.Y3, dom.Z3, C, vis.nivells(i));
    end
    if isempty(f)
        f = [1 2 3];            % triangle invisible (sense superficie)
        v = nan(3, 3);
    end
    set(vis.hIso(i), 'Faces', f, 'Vertices', v);
end
end
