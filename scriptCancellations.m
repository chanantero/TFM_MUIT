%% Preamble
pathSetUp;

dataPathName = [globalPath, 'Data\'];
ID = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% System set up.
obj = SimulationController;

obj.NSposition = [3.35 -0.2 0]; % Assumed real position
obj.amplitude = 0.5;
obj.phase = 0;
obj.frequency = 440;

% Receivers
% receiverPositions = [1.5 3 0; 1.5 1.5 0];
load([dataPathName, 'acousticPathsGTAC_440.mat'])
obj.microPos = microphonePositions;

% GTAC's official acoustic paths
% acousticPath = importImpulseResponseGTAC(frequency);
% load([dataPathName, 'acousticPathsGTAC_440.mat'])
acPathWFSarrayStruct.acousticPaths = acousticPath;
acPathWFSarrayStruct.frequencies = 440;

obj.setAcousticPaths('NS', 'theoretical', 'WFS', acPathWFSarrayStruct);

%% Cancellation attempts

% Optimization options
sourceFilter = {'Loudspeakers', 'Loudspeakers'}; % It makes no sense to optimize a less real scenario, so use loudspeakers filter
maxAbsValCons = [false, false]; % We always want to be realistic about real constraints
acousticPathType = {'Current', 'Current'}; % It makes no sense to optimize with a theoric acoustic path because it depends on the parameters of the noise source, and those parameters are actually unknown. Besides, the acousic path of the loudspeakers is only known in the places where microphones have been placed.
grouping = {'Independent', 'AllTogether'};
zerosFixed = [false, false];

obj.cancel(sourceFilter, maxAbsValCons, acousticPathType, grouping, zerosFixed);
% obj.experimentalChecking(); % Reproduce for experimental checking of
% correspondence

% Visualization
    % Extract 2D map
fMap = figure;
ax = copyobj(obj.ax, fMap);
colormap(ax, 'jet')
colorbar(ax)
    % Use viewer object
s = obj.cancelResults;
% - NScoef
% - WFScoef
% - recNScoef
% - recWFSCoef
% - recCoef
% - recNSCoefExp
% - recWFSCoefExp
% - recCoefExp
    
objVis = simulationViewer(ax, obj.cancelResults);
    % Export 2D map
printfig(fMap, 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\', 'prueba2DMap', 'pdf');

%% Optimization of noise source theoric parameters
% Find the theoric parameters for the virtual noise source that, applying
% WFS cancellation with theoric acoustic path and unified optimization of loudspeakers, minimize the
% magnitude of the resulting field using the experimental acoustic path

obj.findBestVirtualSourceParameters();

