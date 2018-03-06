%% Reproduction tests in laptop

%% Preamble
pathSetUp;

imagesPath = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
obj = SimulationController;

obj.NSposition = [3.35 -0.2 0]; % Assumed real position
obj.amplitude = 1;
obj.phase = 0;
obj.frequency = 440;

% Testing for laptop
obj.WFSToolObj.setNumWFSarraySources(2);
obj.WFSToolObj.setNumReceivers(1);

%% 
% Select adequate audio driver and device for reproduction and recording
obj.WFSToolObj.prepareReproduction();
obj.WFSToolObj.reproduceAndRecord('main', 'soundTime', 2); % Simple reproduction of one pulse of 2 seconds
