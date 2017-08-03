% Create two noise sources. One of them will be virtual. It will be the
% assumed real source to be cancelled. The other will be real and its only
% purpose is to reproduce a signal.
% The real source will remain always the same. The virtual one will change
% in amplitude, phase and position.

obj = WFSToolSimple;

obj.changeScenario(2);

obj.setNumNoiseSources(2);
obj.setVirtual([false; true]);
obj.setReal([true; false]);
obj.noiseSourceChannelMapping = [1; 0];

amplitude = 0.5;
phase = 0;
frequency = 600;

obj.amplitude = [amplitude; amplitude];
obj.phase = [phase; phase];
obj.frequency = [frequency; frequency];

obj.updateReprodPanelBasedOnVariables();

realPosition = [0.6 0 0]; % Assumed real position
obj.noiseSourcePosition = [realPosition; realPosition];

% For the current configuration. Get the experimental acoustic path
obj.reproduceAndRecordForAcousticPaths();
obj.calculateExperimentalAcousticPaths();

% Test the main configuration with prelude
obj.reproduceAndRecord('preludeAndMain', 'numRepetitions', 1);
PathName = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Matlab\Data\';
FileName = 'BasePoint.mat';

s = obj.exportInformation();
save([PathName, FileName], 's');

% Test the multiple cases
minXPos = realPosition(1) - 0.1; maxXPos = realPosition(1) + 0.1; numXPoints = 4;
minYPos = realPosition(2) - 0.1; maxYPos = realPosition(2) + 0.1; numYPoints = 4;
minZPos = 0; maxZPos = 0; numZPoints = 1;

xVec = linspace(minXPos, maxXPos, numXPoints);
yVec = linspace(minYPos, maxYPos, numYPoints);
zVec = linspace(minZPos, maxZPos, numZPoints);
[X, Y, Z] = ndgrid(xVec, yVec, zVec);
virtPos = [X(:), Y(:), Z(:)];

numPoints = size(virtPos, 1);

for p = 1:numPoints
    
    % Set virtual position
    obj.noiseSourcePosition = [realPosition; virtPos(p, :)];
    
    % Apply WFS calculation
    obj.WFScalculation();
    
    % Reproduce
    obj.reproduceAndRecord('main');
    
    % Retrieve and save information
    PathName = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Matlab\Data\';
    FileName = sprintf('%d.mat', p);
    
    s = obj.exportInformation();
    save([PathName, FileName], 's');   
end


