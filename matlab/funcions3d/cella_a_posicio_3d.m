function p = cella_a_posicio_3d(c, dom)
%CELLA_A_POSICIO_3D  Converteix la cel.la [iy ix iz] a [x y z] (m).
p = [dom.x(c(2)) dom.y(c(1)) dom.z(c(3))];
end
