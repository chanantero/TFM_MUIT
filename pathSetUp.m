globalPath = 'C:\Users\Rub�n\Google Drive\Telecomunicaci�n\M�ster 2� Curso 2015-2016\TFM MUIT\Matlab\';
paths = genpath(globalPath);
addpath(paths);

matlabTelecoMasterPath = genpath('C:\Users\Rub�n\Google Drive\Telecomunicaci�n\M�ster 2� Curso 2015-2016\Matlab');
addpath(matlabTelecoMasterPath);

% Clear variables that I don't want
clear paths
clear matlabTelecoMasterPath