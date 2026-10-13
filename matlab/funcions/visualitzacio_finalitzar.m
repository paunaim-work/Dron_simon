function visualitzacio_finalitzar(vis)
%VISUALITZACIO_FINALITZAR  Tanca el fitxer de video, si n'hi ha.
if ~isempty(vis.video)
    close(vis.video);
    fprintf('Video desat.\n');
end
end
