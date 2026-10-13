function vid = video_afegir(vid, fig)
%VIDEO_AFEGIR  Captura la figura i l'afegeix com a fotograma al video.
%   Tots els fotogrames han de tenir la mateixa mida: si la finestra
%   canvia de mida durant la gravacio, la imatge es retalla o s'omple de
%   blanc fins a la mida del primer fotograma.

if isempty(vid)
    return;
end

F = getframe(fig);
img = F.cdata;
if isempty(vid.mida)
    % L'MPEG-4 necessita amplada i alcada parelles
    vid.mida = 2 * floor([size(img, 1) size(img, 2)] / 2);
end
img = ajustar_mida(img, vid.mida);

if strcmp(vid.format, 'gif')
    if vid.mida(2) > 1000
        img = reduir_meitat(img);          % GIF mes lleuger
    end
    [A, mapa] = a_indexada(img);
    if vid.nFotogrames == 0
        imwrite(A, mapa, vid.fitxer, 'gif', 'LoopCount', Inf, 'DelayTime', 1 / vid.fps);
    else
        imwrite(A, mapa, vid.fitxer, 'gif', 'WriteMode', 'append', 'DelayTime', 1 / vid.fps);
    end
else
    writeVideo(vid.writer, img);
end
vid.nFotogrames = vid.nFotogrames + 1;
end

% ---------------------------------------------------------------------
function out = ajustar_mida(img, mida)
out = 255 * ones(mida(1), mida(2), 3, 'uint8');
h = min(mida(1), size(img, 1));
w = min(mida(2), size(img, 2));
out(1:h, 1:w, :) = img(1:h, 1:w, :);
end

% ---------------------------------------------------------------------
function out = reduir_meitat(img)
%REDUIR_MEITAT  Redueix la imatge a la meitat fent la mitjana de cada bloc
%   de 2x2 pixels (la mida es sempre parella).
x = double(img);
out = uint8((x(1:2:end, 1:2:end, :) + x(2:2:end, 1:2:end, :) ...
           + x(1:2:end, 2:2:end, :) + x(2:2:end, 2:2:end, :)) / 4);
end

% ---------------------------------------------------------------------
function [A, mapa] = a_indexada(img)
%A_INDEXADA  Converteix la imatge RGB a 256 colors (necessari per al GIF).
try
    [A, mapa] = rgb2ind(img, 256, 'nodither');
catch
    % Alternativa sense rgb2ind: paleta fixa de 6 x 7 x 6 = 252 colors
    nR = 6; nG = 7; nB = 6;
    [bb, gg, rr] = ndgrid(0:nB - 1, 0:nG - 1, 0:nR - 1);
    mapa = [rr(:) / (nR - 1), gg(:) / (nG - 1), bb(:) / (nB - 1)];
    img = double(img) / 255;
    r = round(img(:, :, 1) * (nR - 1));
    g = round(img(:, :, 2) * (nG - 1));
    b = round(img(:, :, 3) * (nB - 1));
    A = uint8(r * nG * nB + g * nB + b);
end
end
