function vid = video_obrir(cfg, etiqueta)
%VIDEO_OBRIR  Prepara la gravacio del video de l'animacio.
%
%   vid = video_obrir(cfg, etiqueta) retorna [] si cfg.visual.gravarVideo
%   es false. Si no, crea el fitxer a la carpeta cfg.visual.carpetaVideo
%   (per defecte, la carpeta 'videos' dins de la carpeta actual) amb un
%   nom que inclou la ploma, l'algorisme, la data i l'hora, per exemple:
%
%       videos/cerca_2d_gaussiana_infotaxis_20261013_153000.mp4
%
%   Formats (cfg.visual.formatVideo):
%     'mp4' : MPEG-4 (Windows i macOS). Si l'ordinador no el permet (p. ex.
%             MATLAB a Linux), es desa en AVI.
%     'avi' : Motion JPEG AVI (qualsevol ordinador amb MATLAB).
%     'gif' : GIF animat (tambe funciona a GNU Octave, que no te VideoWriter).

vid = [];
if ~cfg.visual.gravarVideo
    return;
end

carpeta = cfg.visual.carpetaVideo;
if isempty(carpeta)
    carpeta = fullfile(pwd, 'videos');
end
if ~exist(carpeta, 'dir')
    mkdir(carpeta);
end

nom = cfg.visual.nomVideo;
if isempty(nom)
    nom = sprintf('cerca_%s_%s_%s_%s', etiqueta, cfg.tipusPloma, cfg.algorisme, ...
                  datestr(now, 'yyyymmdd_HHMMSS'));
end
[~, nom] = fileparts(nom);              % sense extensio

vid.fps         = cfg.visual.fpsVideo;
vid.mida        = [];                   % mida dels fotogrames (es fixa amb el primer)
vid.nFotogrames = 0;
vid.writer      = [];
vid.format      = lower(cfg.visual.formatVideo);
vid.fitxer      = '';

teVideoWriter = exist('VideoWriter', 'class') == 8 || exist('VideoWriter', 'file') == 2;
if ~strcmp(vid.format, 'gif') && teVideoWriter
    if strcmp(vid.format, 'avi')
        perfils = {'Motion JPEG AVI', '.avi'};
    else
        perfils = {'MPEG-4', '.mp4'; 'Motion JPEG AVI', '.avi'};
    end
    for i = 1:size(perfils, 1)
        fitxer = fullfile(carpeta, [nom perfils{i, 2}]);
        try
            w = VideoWriter(fitxer, perfils{i, 1});
            w.FrameRate = vid.fps;
            w.Quality   = 95;
            open(w);
            vid.writer = w;
            vid.format = perfils{i, 2}(2:end);
            vid.fitxer = fitxer;
            break;
        catch err
            fprintf('No es pot gravar en format %s (%s).\n', perfils{i, 1}, err.message);
        end
    end
end

if isempty(vid.writer)
    if ~strcmp(vid.format, 'gif')
        fprintf('Es desara el video com a GIF animat.\n');
    end
    vid.format = 'gif';
    vid.fitxer = fullfile(carpeta, [nom '.gif']);
    if exist(vid.fitxer, 'file')
        delete(vid.fitxer);
    end
end
fprintf('Gravant el video a:\n   %s\n', vid.fitxer);
end
