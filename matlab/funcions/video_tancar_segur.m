function video_tancar_segur(vid)
%VIDEO_TANCAR_SEGUR  Tanca el fitxer de video sense missatges. Es crida
%   automaticament (onCleanup) si la simulacio s'atura per un error o amb
%   Ctrl+C, perque el video gravat fins aquell moment no quedi malmes.
if ~isempty(vid) && ~isempty(vid.writer)
    try
        close(vid.writer);
    catch
    end
end
end
