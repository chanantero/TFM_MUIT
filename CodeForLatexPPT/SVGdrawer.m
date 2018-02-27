classdef SVGdrawer
    
    properties
        % Positions and angles
        NSpositions
        NSangles
        WFSpositions
        WFSangles
        microPositions
        microAngles
        
        % Sizes
        microSize
        
        % Colors
        microColor % Cell string. Hexadecimal format.
        WFScolor
        NScolor
        
        % Symbols
        NSsymbol
        WFSsymbol
        microSymbol
        
        % Diagrams
        WFSdiagramIndices
        WFSdiagramNames
        NSdiagramIndices
        NSdiagramNames
        microDiagramIndices
        microDiagramNames
        
        % Other
        backgroundFileName
        viewBox
        
    end
    
    properties(Dependent)
        numNS
        numWFS
        numMicro
        backgroundFlag
        microColorUnique
        WFScolorUnique
        NScolorUnique
    end
    
    % Getters and setters
    methods
        function numNS = get.numNS(obj)
            numNS = size(obj.NSpositions, 1);
        end
        
        function numWFS = get.numWFS(obj)
            numWFS = size(obj.WFSpositions, 1);
        end
        
        function numMicro = get.numMicro(obj)
            numMicro = size(obj.microPositions, 1);
        end
        
        function backgroundFlag = get.backgroundFlag(obj)
            backgroundFlag = ~isempty(obj.backgroundFileName);
        end
        
        function microColorUnique = get.microColorUnique(obj)
            microColorUnique = numel(obj.microColor) ~= obj.numMicro;
        end
        
        function WFScolorUnique = get.WFScolorUnique(obj)
            WFScolorUnique = numel(obj.WFScolor) ~= obj.numWFS;
        end
        
        function NScolorUnique = get.NScolorUnique(obj)
            NScolorUnique = numel(obj.NScolor) ~= obj.numNS;
        end
    end
    
    methods
        function obj = SVGdrawer(varargin)
            p = inputParser;
            
            addParameter(p, 'viewBox', [-1 -1 5 8]);
            addParameter(p, 'NSpositions', [0 0]);
            addParameter(p, 'NSangles', 0);
            addParameter(p, 'microPositions', double.empty(0, 2));
            addParameter(p, 'WFSdiagramIndices', []);
            addParameter(p, 'WFSdiagramNames', []);
            addParameter(p, 'microDiagramIndices', []);
            addParameter(p, 'microDiagramNames', []);
            addParameter(p, 'backgroundFileName', []);
            addParameter(p, 'WFSpositions', []);
            addParameter(p, 'WFSangles', []);
            addParameter(p, 'NSsymbol', 'loudspeakerSound');
            addParameter(p, 'WFSsymbol', 'loudspeaker');
            addParameter(p, 'microSymbol', 'microphone');
            addParameter(p, 'microSize', 0.8);
            addParameter(p, 'microColor', {'000000'});
            addParameter(p, 'WFScolor', {'000000'});
            addParameter(p, 'NScolor', {'000000'});
            
            parse(p, varargin{:});
            
            obj.viewBox = p.Results.viewBox;
            obj.NSpositions = p.Results.NSpositions;
            obj.NSangles = p.Results.NSangles;
            obj.microPositions = p.Results.microPositions;
            obj.WFSdiagramIndices = p.Results.WFSdiagramIndices;
            obj.WFSdiagramNames = p.Results.WFSdiagramNames;
            obj.microDiagramIndices = p.Results.microDiagramIndices;
            obj.microDiagramNames = p.Results.microDiagramNames;
            obj.backgroundFileName = p.Results.backgroundFileName;
            obj.NSsymbol = p.Results.NSsymbol;
            obj.WFSsymbol = p.Results.WFSsymbol;
            obj.microSymbol = p.Results.microSymbol;
            obj.microSize = p.Results.microSize;
            obj.microColor = p.Results.microColor;
            obj.WFScolor = p.Results.WFScolor;
            obj.NScolor = p.Results.NScolor;
            
            if ismember('WFSpositions', p.UsingDefaults)
                % Generate positions and orientations of WFS array loudspeakers
                d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
                nb = 8; % Bottom and upper sides of the octogon (2 sides)
                nd = 8; % Diagonal sides of the octogon (4 sides)
                nl = 24; % Lateral side of the octogon (2 sides)
                betabd = 45; % Deviation angle between bottom/upper and diagonal sides
                [ xLoud, yLoud, angleLoud ] = octogon(d, nb, nd, nl, betabd);
                obj.WFSpositions = [xLoud, yLoud];
                obj.WFSangles = angleLoud;
            else
                obj.WFSpositions = p.Results.WFSpositions;
                obj.WFSangles = p.Results.WFSangles;
            end
            
        end
        
        function svgText = getSVG(obj)
            
            % Position and orientation of radiation diagrams
            xDiag = [obj.NSpositions(obj.NSdiagramIndices, 1);...
                obj.microPositions(obj.microDiagramIndices, 1); ...
                obj.WFSpositions(obj.WFSdiagramIndices, 1)];
            yDiag = [obj.NSpositions(obj.NSdiagramIndices, 2);...
                obj.microPositions(obj.microDiagramIndices, 2); ...
                obj.WFSpositions(obj.WFSdiagramIndices, 2)];
            angleDiag = [obj.NSangles(obj.NSdiagramIndices); zeros(numel(obj.microPositions), 1); obj.WFSangles(obj.WFSdiagramIndices)];
            diagName = [obj.NSdiagramNames; obj.microDiagramNames; obj.WFSdiagramNames];
            numDiagrams = numel(xDiag);
            
            % Viewbox and viewport
            viewPort2viewBoxRatio = 100;
            width = obj.viewBox(3)*viewPort2viewBoxRatio;
            height = obj.viewBox(4)*viewPort2viewBoxRatio;
            
            % Create string that will draw them in SVG
            drawLoudspeakerStr = '<use xlink:href="symbolLibrary.svg#%s" x="%g" y="%g" transform="rotate(%g %g %g)" style="fill:%s"/>\n';
            drawMicrophoneStr = '<use xlink:href="symbolLibrary.svg#%s" transform="translate(%g %g) scale(%g)" style="fill:%s"/>\n'; % Example: fill:#AA00FF. Hashtag is necessary
            drawDiagramStr = '<use xlink:href="symbolLibrary.svg#%s" x="%g" y="%g" transform="rotate(%g %g %g)" />\n';
            drawBackgroundStr = '<image xlink:href="%s" x="%g" y="%g" width="%g" height="%g" preserveAspectRatio="none"/>';
            
            if obj.backgroundFlag
                strBackground = sprintf(drawBackgroundStr, obj.backgroundFileName, obj.viewBox(1), obj.viewBox(2), obj.viewBox(3), obj.viewBox(4));
            else
                strBackground = '';
            end
            
            strLoud = cell(obj.numWFS, 1);
            for k = 1:obj.numWFS
                if ~obj.WFScolorUnique
                    strLoud{k} = sprintf(drawLoudspeakerStr, obj.WFSsymbol, obj.WFSpositions(k, 1), obj.WFSpositions(k, 2), obj.WFSangles(k), obj.WFSpositions(k, 1), obj.WFSpositions(k, 2),...
                        obj.WFScolor{k});
                else
                    strLoud{k} = sprintf(drawLoudspeakerStr, obj.WFSsymbol, obj.WFSpositions(k, 1), obj.WFSpositions(k, 2), obj.WFSangles(k), obj.WFSpositions(k, 1), obj.WFSpositions(k, 2),...
                    obj.WFScolor{1});
                end
            end
            strLoud = strjoin(strLoud);
            strLoud = [sprintf('<g>\n'), strLoud, sprintf('</g>\n')];
            
            strLoudNoise = cell(obj.numNS, 1);
            for k = 1:obj.numNS
                if ~obj.NScolorUnique
                    strLoudNoise{k} = sprintf(drawLoudspeakerStr, obj.NSsymbol, obj.NSpositions(k, 1), obj.NSpositions(k, 2), obj.NSangles(k), obj.NSpositions(k, 1), obj.NSpositions(k, 2),...
                            obj.WFScolor{k});   
                else
                    strLoudNoise{k} = sprintf(drawLoudspeakerStr, obj.NSsymbol, obj.NSpositions(k, 1), obj.NSpositions(k, 2), obj.NSangles(k), obj.NSpositions(k, 1), obj.NSpositions(k, 2),...
                            obj.WFScolor{1});
                end
            end
            strLoudNoise = strjoin(strLoudNoise);
            
            strMicro = cell(obj.numMicro, 1);
            for k = 1:obj.numMicro
                if ~obj.microColorUnique
                    strMicro{k} = sprintf(drawMicrophoneStr, obj.microSymbol, obj.microPositions(k, 1), obj.microPositions(k, 2), obj.microSize, obj.microColor{k});
                else
                    strMicro{k} = sprintf(drawMicrophoneStr, obj.microSymbol, obj.microPositions(k, 1), obj.microPositions(k, 2), obj.microSize, obj.microColor{1});
                end
            end
            strMicro = strjoin(strMicro);
            
            strDiag = cell(numDiagrams, 1);
            for k = 1:numDiagrams
                strDiag{k} = sprintf(drawDiagramStr, diagName{k}, xDiag(k), yDiag(k), angleDiag(k), xDiag(k), yDiag(k));
            end
            strDiag = strjoin(strDiag);
            
            % Read file
            file = fopen('C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Documentos\Img\WFSarrayScheme_Template.svg');
            svgText = fread(file, Inf, '*char')';
            fclose(file);
            
            % Set viewport and viewbox
            svgText = strrep(svgText, '[viewPortWidth]', num2str(width));
            svgText = strrep(svgText, '[viewPortHeight]', num2str(height));
            svgText = strrep(svgText, '[viewBox]', num2str(obj.viewBox));
            
            % Insert background image
            svgText = strrep(svgText, '[background]', strBackground);
            
            % Insert loudspeakers in the desired positions and orientations
            svgText = strrep(svgText, '[WFSarray]', strLoud);
            
            % Noise Loudspeakers
            svgText = strrep(svgText, '[NoiseLouds]', strLoudNoise);
            
            % Insert string to draw microphones
            svgText = strrep(svgText, '[MicrophoneArray]', strMicro);
            
            % Insert diagrams
            svgText = strrep(svgText, '[RadDiag]', strDiag);
            
        end
        
        function drawSVG(obj, fileName)
            svgText = obj.getSVG();
            
            % Write file
            destFile = fopen(fileName, 'w', 'n', 'UTF-8');
            fwrite(destFile, svgText, 'char');
            fclose(destFile);
        end
    end
    
end

