% Based on a SVG template with variables of the form "[variableName]",
% substitute those variables for the actual values.
% It was intended for embedding SVGs inside other SVGs in order to authomatize procecess,
% but I've given up that, it presents some problems I can't solve for now.

% Set size variables
    % Variables chosen by the user
    
        % Dimensions
        gridSize = [2, 3]; % [height, width]
        index = [1 2 3; 4 5 6];
        gridCellWidth = 250;
        gridCellHeight = 200;
        frame2CellRatio_width = 0.9; frame2CellRatio_height = 0.9;
        imag2frameRatio_width = 0.9; imag2frameRatio_height = 0.9; % Proportion of the frame that the image occupies
        imagOriginalWidth = 560;
        imagOriginalHeight = 420;
    
        % Name of images
        path = 'C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\';
        name1 = 'pulseSignal.svg';
        name2 = 'DFT_channel2.svg';
        name3 = 'DFT_channel2_filter.svg';
        name4 = 'IQ.svg';
        name5 = 'correlation.svg';
        name6 = 'correlation.svg';
        names = {name1, name2, name3, name4, name5, name6};
           
    % Automatically set variables based on the user ones
    numImages = sum(index(:) ~= 0);
    
    rectangleWidth = gridCellWidth * frame2CellRatio_width; % Image frame width
    rectangleHeight = gridCellHeight * frame2CellRatio_height; % Image frame height
    imagWidth = rectangleWidth*imag2frameRatio_width;
    imagHeight = rectangleHeight*imag2frameRatio_height;
    
    frameXPos = zeros(numImages, 1);
    frameYPos = zeros(numImages, 1);
    for k = 1:numImages
        [row, column] = find(index == k);
        frameXPos(k) = (column - 1)*gridCellWidth + gridCellWidth*(1 - frame2CellRatio_width)/2;
        frameYPos(k) = (row - 1)*gridCellHeight + gridCellHeight*(1 - frame2CellRatio_height)/2;
    end
    
    imagXPos = frameXPos + (rectangleWidth - imagWidth)/2;
    imagYPos = frameYPos + (rectangleHeight - imagHeight)/2;
    
    for k = 1:numImages
        varName = ['xRect', num2str(k)];
        value = frameXPos(k);
        eval([varName, ' = value;']);
        
        varName = ['yRect', num2str(k)];
        value = frameYPos(k);
        eval([varName, ' = value;']);

        varName = ['xImag', num2str(k)];
        value = imagXPos(k);
        eval([varName, ' = value;']);
        
        varName = ['yImag', num2str(k)];
        value = imagYPos(k);
        eval([varName, ' = value;']);
        
        varName = ['imagName', num2str(k)];
        name = names{k};
        eval([varName, ' = name;']);
        
        varName = ['scaleFactorX', num2str(k)];
        value = imagWidth/imagOriginalWidth;
        eval([varName, ' = value;']);
        
        varName = ['scaleFactorY', num2str(k)];
        value = imagHeight/imagOriginalHeight;
        eval([varName, ' = value;']);
    end

% Read file
file = fopen('C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\EsquemaDeteccionSenyalesPlantilla.svg');
svgTemplate = fread(file, Inf, '*char')';
fclose(file);

% Find variable names
[start, ending, tokens] = regexp(svgTemplate, '\[(.+?)\]', 'start', 'end', 'tokens'); % Find fields
numVar = numel(start);
    % Identify variable names
    varNames = cell(numVar, 1);
    for k = 1:numVar
        varNames{k} = tokens{k}{1};
    end

% Sustitute variable names by the actual variable values
svgText = svgTemplate;
for k = 1:numVar
    svgText = strrep(svgText, ['[', varNames{k}, ']'], num2str(eval(varNames{k})));
end

% Write file
destFile = fopen('C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\EsquemaDeteccionSenyales.svg', 'w', 'n', 'UTF-8');
fwrite(destFile, svgText, 'char');
fclose(destFile);


