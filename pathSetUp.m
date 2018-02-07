globalPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Matlab\';
paths = genpath(globalPath);
addpath(paths);

matlabTelecoMasterPath = genpath('C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\Matlab');
addpath(matlabTelecoMasterPath);

% Clear variables that I don't want
clear paths
clear matlabTelecoMasterPath