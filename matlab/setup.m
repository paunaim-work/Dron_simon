% Afegeix les carpetes del projecte al path de MATLAB.
arrel = fileparts(mfilename('fullpath'));
addpath(genpath(arrel));
disp('Dron_simon: path configurat.');
