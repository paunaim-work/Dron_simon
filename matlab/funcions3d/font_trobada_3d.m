function t = font_trobada_3d(p, posFont, dom)
%FONT_TROBADA_3D  true si el dron (a p = [x y z]) es dins de la zona de
%   deteccio de la fuga: a menys de dom.radiDet en horitzontal i a menys
%   de dom.tolAlcada en vertical.
t = norm(p(1:2) - posFont(1:2)) <= dom.radiDet + 1e-9 && ...
    abs(p(3) - posFont(3)) <= dom.tolAlcada + 1e-9;
end
