function visualitzacio_finalitzar(vis)
%VISUALITZACIO_FINALITZAR  Tanca el video, si se n'esta gravant un. Abans
%   hi repeteix l'ultim fotograma 2 s perque es pugui llegir el resultat.
if ~isempty(vis.video) && ~vis.tancada && ishandle(vis.fig)
    for i = 1:round(2 * vis.video.fps)
        vis.video = video_afegir(vis.video, vis.fig);
    end
end
video_tancar(vis.video);
end
