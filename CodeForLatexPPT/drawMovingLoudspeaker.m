%% Draw Moving Loudspeaker
% Parameters. User can touch this.

    % Noise loudspeaker
    realPosition = [-1 2 0]; % Assumed real position of noise source
    amplitude = 0.5;
    phase = 0;
    frequency = 440;
    angleNoise = 0; noiseName = 'loudspeakerSound';
    radius = 0.8;
    
    % Receivers
    receiverPositions = [1.5 3 0; 1.5 1.5 0];
    numReceivers = size(receiverPositions, 1);
   
    % Position and orientation of radiation diagrams
    loudsIndex = 6; microphIndex = 2;
    
    % Viewbox and viewport
    width = 600; height = 800;
    viewPort2viewBoxRatio = 100;
    
    % GIF file
    gifDuration = 4; % seconds
    numFrames = 120;
    
    % Folder
    pathFold = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';

%% User don't touch!

% Calculate derivated variables from user parameters

    % Fixed parameters for the WFSTool object
    numLoudspeakers = 96;
        % Predefined acoustic paths.
        theoricWFSAcPath = true;
        theoricNoiseSourceAcPath = true;

    % Generate positions and orientations of WFS array loudspeakers
    d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
    nb = 8; % Bottom and upper sides of the octogon (2 sides)
    nd = 8; % Diagonal sides of the octogon (4 sides)
    nl = 24; % Lateral side of the octogon (2 sides)
    betabd = 45; % Deviation angle between bottom/upper and diagonal sides
    [ xLoud, yLoud, angleLoud ] = octogon(d, nb, nd, nl, betabd);

    % Position and orientation of noise loudspeaker
    noiseVirtPos = movingLoudspeakerPos( 'spiral', 'radius', radius, 'numRounds', 2, 'numPoints', numFrames, 'refPos', [realPosition(1), realPosition(2)] );
    xNoise = noiseVirtPos(:, 1); yNoise = noiseVirtPos(:, 2);

    % Position and orientation of radiation diagrams
    xDiag = [xLoud(6); receiverPositions(2, 1)]; yDiag = [yLoud(6), receiverPositions(2, 2)]; angleDiag = [angleLoud(6); 0];
    diagName = {'diagramDipole', 'diagramMonopoleNoise'};

    % Viewbox and viewport
    viewBoxWidth = width/viewPort2viewBoxRatio;
    viewBoxHeight = height/viewPort2viewBoxRatio;

    rangeX = max(xLoud) - min(xLoud);
    rangeY = max(yLoud) - min(yLoud);

    marginX = (viewBoxWidth - rangeX)*0.75;
    marginY = (viewBoxHeight - rangeY)/2;

    viewBox = [min(xLoud) - marginX, min(yLoud) - marginY, viewBoxWidth, viewBoxHeight];

% Create string that will draw the WFS loudspeaker array in SVG
drawLoudspeakerStr = '<use xlink:href="symbolLibrary.svg#%s" x="%g" y="%g" transform="rotate(%g %g %g)"/>\n';
strLoud = '';
numLouds = numel(xLoud);
for k = 1:numLouds
    strLoud = [strLoud, sprintf(drawLoudspeakerStr, 'loudspeaker', xLoud(k), yLoud(k), angleLoud(k), xLoud(k), yLoud(k))];
end
strLoud = [sprintf('<g>\n'), strLoud, sprintf('</g>\n')];

% Change matlab current folder to the one that contains the image files
% (.svg, etc.) so it doesn't cause any problem
cd(pathFold);

% Set Up System
    if exist('obj', 'var') == 0 || ~isvalid(obj) || ~isvalid(obj.ax)
        obj = WFSToolSimple;
    end

    % Set noise source variables
    obj.setNumNoiseSources(2);
    obj.noiseSourceChannelMapping = [1; 0];
    obj.amplitude = [amplitude; amplitude];
    obj.phase = [phase; phase];
    obj.frequency = [frequency; frequency];
    obj.noiseSourcePosition = [realPosition; realPosition];
    obj.setVirtual([false; true]);
    obj.setReal([true; false]);

    % Set WFS array variables
    obj.setNumWFSarraySources(numLoudspeakers);

    % Set receiver variables
    obj.setNumReceivers(numReceivers);
    obj.receiverPosition = receiverPositions;

    % Set acoustic path variables
    if theoricWFSAcPath
        obj.theoricWFSacousticPath();
    else
        obj.WFSarrayAcPathStruct.acousticPaths = acPathWFSarray;
    end

    if theoricNoiseSourceAcPath
        obj.theoricNoiseSourceAcousticPath();
    else
        obj.noiseSourceAcPathStruct.acousticPaths = acPathNoiseSources;
    end

    obj.updateReprodPanelBasedOnVariables();
    obj.updateRecordPanelBasedOnVariables();
    
% Configure simulator
    xLim = [viewBox(1), viewBox(1) + viewBox(3)];
    yLim = [viewBox(2), viewBox(2) + viewBox(4)];
    obj.simulTheo.XLim = xLim;
    obj.simulTheo.YLim = yLim;
    obj.simulTheo.XnumPoints = floor(width/2);
    obj.simulTheo.YnumPoints = floor(height/2);
    obj.simulTheo.generateMeasurePoints();
    obj.simulTheo.updateTheoricAcousticPathsImage();

% Generate image sequence
name = 'movingLoudspeaker';
loudsColormap = [0 0 0; 1 1 1];
fieldColormapBasic = [1 1 1; 0 0 1];
fieldColormap = extendColormap( fieldColormapBasic, 128 );
combColormap = [loudsColormap; fieldColormap];
imagGIF = zeros(height, width, 1, numFrames);
for f = 1:numFrames
    fprintf('Frame %d/%d\n', f, numFrames)
    
    % Draw svg image and import it as bitmap
        strLoudNoise =  sprintf(drawLoudspeakerStr, noiseName, xNoise(f), yNoise(f), angleNoise, xNoise(f), yNoise(f));

        % Read file
        file = fopen([pathFold, 'WFSarrayScheme_Template.svg']);
        svgText = fread(file, Inf, '*char')';
        fclose(file);

        % Set viewport and viewbox
        svgText = strrep(svgText, '[viewPortWidth]', num2str(width));
        svgText = strrep(svgText, '[viewPortHeight]', num2str(height));
        svgText = strrep(svgText, '[viewBox]', num2str(viewBox));

        % Insert loudspeakers in the desired positions and orientations
        svgText = strrep(svgText, '[WFSarray]', strLoud);

        % Noise Loudspeakers
        svgText = strrep(svgText, '[NoiseLouds]', strLoudNoise);

        % Write file 
        destFile = fopen([pathFold, name, '.svg'], 'w', 'n', 'UTF-8');
        fwrite(destFile, svgText, 'char');
        fclose(destFile);

        % Export to PNG and EMF
        % Matlab should be in the folder with the .svg to avoid path problems
        system(['inkscape -z "', pathFold, name, '.svg"', ' -e "', pathFold, name, '.png"', ' --export-width=600', '--export-background-opacity=0']);

        [~, ~, alpha] = imread([pathFold, name, '.png'], 'png');
                
        % imWFSArray = changePixelColor(imWFSArray, [1 1 1], transparency == 0);
        %     
        % % Convert to indexed image
        % X = rgb2ind(imWFSArray, map);
        % X = X + 1;
    
        X = uint8(alpha == 0);
        loudsRGB = ind2rgb(X, loudsColormap);
    
    % Simulate field and draw it
        obj.noiseSourcePosition(2, :) = noiseVirtPos(f, :);
        obj.WFScalculation('Grouping', 'AllTogether', 'SourceFilter', 'NoFilter', 'AcousticPath', 'Theoric');
        obj.simulTheo.updateFieldImage();
        U = obj.simulTheo.fieldImage;
        U = reshape(sum(U, 2), obj.simulTheo.XnumPoints, obj.simulTheo.YnumPoints).';
        
        indIm = scaled2indexedColors(128, [0 1], abs(U));
        fieldRGB = ind2rgb(indIm, fieldColormap);
        
        fieldRGB = imresize(fieldRGB, [height width]);
    
    % Put loudspeakers over other image
    H = vision.AlphaBlender('Operation', 'Binary mask', 'MaskSource', 'Input port');    
    comb = step(H, fieldRGB, loudsRGB, alpha > 0);            
    combInd = rgb2ind(comb, combColormap) + 1;
    
    imagGIF(:, :, 1, f) = combInd;  
end

% Create GIF file
numFramesPerPeriod = numFrames;

imwrite(imagGIF, combColormap, [pathFold, name, '.gif'], ...
    'LoopCount', 100, 'DelayTime', gifDuration/numFramesPerPeriod);
   