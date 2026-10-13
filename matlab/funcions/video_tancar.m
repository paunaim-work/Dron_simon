function video_tancar(vid)
%VIDEO_TANCAR  Tanca el fitxer de video (imprescindible perque el .mp4 o
%   .avi es pugui reproduir) i diu on s'ha desat.
if isempty(vid)
    return;
end
if ~isempty(vid.writer)
    try
        close(vid.writer);
    catch
    end
end
if vid.nFotogrames > 0
    fprintf('Video desat (%d fotogrames, %.0f s):\n   %s\n', vid.nFotogrames, ...
            vid.nFotogrames / vid.fps, vid.fitxer);
else
    fprintf('No s''ha gravat cap fotograma.\n');
end
end
