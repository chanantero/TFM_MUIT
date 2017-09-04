classdef simulator < handle
    
    properties
        
        % Sources
        sourcePositions
        sourceCoefficients % (numSources x numFrequencies) matrix
        sourceOrientations % Mx4 matrix. Rotation vector: [angle of rotation, Xaxis, Yaxis, Zaxis]. Counter clock-wise
        radPatFuns
        
        % Receivers
        recPositions
        recRadPatFuns
        recOrientations
        
        % Other
        freq
        acPath
        
        % Graphical
        ax
        XLim
        YLim
        XnumPoints
        YnumPoints
        measurePointsImage
        acPathImage
        z
        drawingOption % image or scatter
        
    end
    
    properties(SetAccess = private)
        field % Result of simulation
        imag % Result of simulation
        fieldImage
    end
    
    properties(Constant)
        c = 340;
    end
    
    properties(Dependent)
        numSources
        numFrequencies
        numRec
        numPointsImage
        k % Propagation constant. numFrequencies-element vector
    end
    
    % Getters and setters
    methods
        function numSources = get.numSources(obj)
            numSources = size(obj.sourcePositions, 1);
        end
        
        function numRec = get.numRec(obj)
            numRec = size(obj.recPositions, 1);
        end
        
        function numPointsImage = get.numPointsImage(obj)
            numPointsImage = size(obj.measurePointsImage, 1);
        end
        
        function numFrequencies = get.numFrequencies(obj)
            numFrequencies = numel(obj.freq);
        end
        
        function k = get.k(obj)
            k = 2*pi*obj.freq/obj.c;
        end
    end
    
    methods
        
        function obj = simulator(parentFig)
            if nargin == 1
                obj.ax = axes(parentFig, 'CLim', [-0.1 0.1]);
                colormap(obj.ax, 'gray');
            end
            obj.setDefaultProperties();
        end
        
        function setDefaultProperties(obj)
            sourcePos = [0 0 0];
            sourceCoeff = 1;
            theta = 0; rotAxis = [1 0 0];
            sourceOrient = [theta, rotAxis];
            radPatFun = {@simulator.monopoleRadPat};
            obj.setSources(sourcePos, 'coefficient', sourceCoeff, 'orientation', sourceOrient,...
                'radiationPattern', radPatFun);
            
            recPos = [0 0 0];
            theta = 0; rotAxis = [1 0 0];
            recOrient = [theta, rotAxis];
            radPatFun = {@simulator.monopoleRadPat};
            obj.setReceivers(recPos, 'orientation', recOrient,...
                'radiationPattern', radPatFun);
            
            obj.freq = 340;
            obj.XLim = [-2 2];
            obj.YLim = [-2 2];
            obj.XnumPoints = 20;
            obj.YnumPoints = 20;
            obj.z = 0;
            obj.drawingOption = 'image';
        end
        
        function simulate(obj)
            %             obj.generateMeasurePoints();
            
            % Simulate
            obj.updateTheoricAcousticPaths();
            obj.updateField();
            
            obj.draw();
        end
        
        function generateMeasurePoints(obj)
            % Create measure points
            obj.measurePointsImage = plane(obj.XLim, obj.YLim, obj.XnumPoints, obj.YnumPoints, [], [], []);
        end
        
        function draw(obj)
            if ~isempty(obj.imag) && isvalid(obj.imag);
                type = obj.imag.Type;
            else
                type = 'invalid';
            end
            
            switch obj.drawingOption
                case 'image'
                    % Reshape as an image                    
                    U = reshape(sum(obj.fieldImage, 2), obj.XnumPoints, obj.YnumPoints).';
                    if strcmp(type, 'image')
                        % Change color according to that
                        obj.imag.CData = real(U);
                    else
                        if strcmp(type, 'scatter')
                            delete(obj.imag)
                        end
                        obj.ax.NextPlot = 'Add';
                        obj.imag = image(obj.ax, 'XData', obj.XLim, 'YData', obj.YLim,...
                            'CData', real(U), 'CDataMapping', 'scaled');
                        obj.imag.HitTest = 'off';
                        obj.ax.Children = obj.ax.Children([2:end, 1]);
                        obj.ax.NextPlot = 'Replace';
                    end
                    
                case 'scatter'
                    
                    U = sum(obj.field, 2);
                    if strcmp(type, 'scatter')
                        % Change color according to that
                        obj.imag.CData = real(U);
                    else
                        if strcmp(type, 'image')
                            delete(obj.imag)
                        end
                        obj.ax.NextPlot = 'Add';
                        obj.imag = scatter(obj.ax, obj.recPositions(:, 1), obj.recPositions(:, 2), 50, real(U),...
                            'MarkerEdgeColor', [0 0 0]);
                        obj.imag.HitTest = 'off';
                        obj.ax.Children = obj.ax.Children([2:end, 1]);
                        obj.ax.NextPlot = 'Replace';
                    end
                    
            end
        end
        
        function animate(obj, t, stretchingFactor)
            t_ = t*stretchingFactor;
            
            incrT = diff(t_);
            incrT = [incrT, incrT(end)];
            
            for l = 1:numel(t)
                U = obj.field.*repmat(exp(1i*2*pi*obj.freq'*t(l)), [obj.numMeasurePoints, 1]);
                U = reshape(sum(U, 2), obj.XnumPoints, obj.YnumPoints).';
                obj.imag.CData = real(U);
                pause(incrT(l))
            end
        end
        
        function updateField(obj)
            
            U = zeros(obj.numRec, obj.numFrequencies);
            for f = 1:obj.numFrequencies
                U(:, f) = (obj.acPath(:, :, f) * obj.sourceCoefficients(:, f));
            end
            
            obj.field = U;
            
        end
        
        function updateFieldImage(obj)
            U = zeros(obj.numPointsImage, obj.numFrequencies);
            for f = 1:obj.numFrequencies
                U(:, f) = (obj.acPathImage(:, :, f) * obj.sourceCoefficients(:, f));
            end
            
            obj.fieldImage = U;
        end
        
        function updateTheoricAcousticPaths(obj)
            obj.acPath = simulator.calculateTheoricAcousticPaths(obj.sourcePositions, obj.radPatFuns, obj.sourceOrientations, obj.recPositions, obj.recRadPatFuns, obj.recOrientations, obj.freq, obj.c);
        end
        
        function updateTheoricAcousticPathsImage(obj)
            recRadPat = repmat({@simulator.monopoleRadPat}, [obj.numPointsImage, 1]);
            recOrient = repmat([0 0 0 1], [obj.numPointsImage, 1]);
                    
            obj.acPathImage = simulator.calculateTheoricAcousticPaths(...
                obj.sourcePositions, obj.radPatFuns, obj.sourceOrientations,...
                obj.measurePointsImage, recRadPat, recOrient,...
                obj.freq, obj.c);        
        end
        
        function setSources(obj, varargin)
            
            p = inputParser;
            
            addRequired(p, 'position');
            addParameter(p, 'coefficient', []);
            addParameter(p, 'orientation', []);
            addParameter(p, 'radiationPattern', []);
            
            parse(p, varargin{:})
            %             disp('Parse completed')
            
            obj.sourcePositions = p.Results.position;
            N = obj.numSources;
            coefficient = p.Results.coefficient;
            orientation = p.Results.orientation;
            radiationPattern = p.Results.radiationPattern;
            
            if ismember('coefficient', p.UsingDefaults) || any(size(coefficient) ~= [N, 1])
                obj.sourceCoefficients = ones(N, 1);
            else
                obj.sourceCoefficients = coefficient;
            end
            
            if ismember('orientation', p.UsingDefaults) || any(size(orientation) ~= [N, 4])
                obj.sourceOrientations = repmat([0 1 0 0], obj.numSources, 1);
            else
                obj.sourceOrientations = orientation;
            end
            
            if ismember('radiationPattern', p.UsingDefaults) || any(size(radiationPattern) ~= [N, 1])
                obj.radPatFuns = cell(obj.numSources, 1);
                for s = 1:N
                    obj.radPatFuns{s} = @(x) ones(size(x, 1), 1);
                end
            else
                obj.radPatFuns = radiationPattern;
            end
            
        end
        
        function setReceivers(obj, varargin)
            
            p = inputParser;
            
            addRequired(p, 'position');
            addParameter(p, 'orientation', []);
            addParameter(p, 'radiationPattern', []);
            
            parse(p, varargin{:})
            
            obj.recPositions = p.Results.position;
            N = obj.numRec;
            orientation = p.Results.orientation;
            radiationPattern = p.Results.radiationPattern;
            
            if ismember('orientation', p.UsingDefaults) || any(size(orientation) ~= [N, 4])
                obj.recOrientations = repmat([0 1 0 0], obj.numSources, 1);
            else
                obj.recOrientations = orientation;
            end
            
            if ismember('radiationPattern', p.UsingDefaults) || any(size(radiationPattern) ~= [N, 1])
                obj.recRadPatFuns = cell(obj.numSources, 1);
                for s = 1:N
                    obj.recRadPatFuns{s} = @(x) ones(size(x, 1), 1);
                end
            else
                obj.recRadPatFuns = radiationPattern;
            end
            
        end
        
        function setKirchhoffIntegralScenario(obj, sourcePos, sourceCoeff, surfacePointsPos, surfacePointsNormal, surfacePointsArea)
            [monopoleCoeff, dipoleCoeff] = simulator.KirchhoffIntegralCoeff(obj.k, sourcePos, sourceCoeff, surfacePointsPos, surfacePointsNormal, surfacePointsArea);
            
            numSurfacePoints = size(surfacePointsPos, 1);
            numSourc = size(sourcePos, 1);
            
            originalSourcePos = sourcePos;
            monopoleSourcePos = surfacePointsPos;
            dipoleSourcePos = surfacePointsPos;
            sourcePos = [originalSourcePos; monopoleSourcePos; dipoleSourcePos];
            
            originalSourceCoeff = -sourceCoeff;
            monopoleSourceCoeff = monopoleCoeff;
            dipoleSourceCoeff = dipoleCoeff;
            sourceCoef = [originalSourceCoeff; monopoleSourceCoeff; dipoleSourceCoeff];
            
            originalSourceOrient = repmat([0 1 0 0], numSourc, 1);
            monopoleSourceOrient = simulator.vec2rotVec(surfacePointsNormal);
            dipoleSourceOrient = simulator.vec2rotVec(surfacePointsNormal);
            sourceOrient = [originalSourceOrient; monopoleSourceOrient; dipoleSourceOrient];
            
            originalSourceRadPat = repmat({@(x) simulator.monopoleRadPat(x)}, numSourc, 1);
            monopoleSourceRadPat = repmat({@(x) simulator.monopoleRadPat(x)}, numSurfacePoints, 1);
            dipoleSourceRadPat = repmat({@(x) simulator.dipolePreciseRadPat(x, obj.k)}, numSurfacePoints, 1);
            sourceRadPat = [originalSourceRadPat; monopoleSourceRadPat; dipoleSourceRadPat];
            
            obj.setSources(sourcePos, 'coefficient', sourceCoef, 'orientation', sourceOrient, 'radiationPattern', sourceRadPat);
        end
        
        function setLinearArrayScenario(obj, sourcePos, sourceCoeff, surfacePointsPos)
            numSour = size(sourcePos, 1);
            numMeas = size(surfacePointsPos, 1);
            
            radPat = repmat({@(x) simulator.monopoleRadPat(x)}, numSour, 1);
            aux = simulator;
            aux.setSources(sourcePos, 'coefficient', sourceCoeff, 'radiationPattern', radPat);
            U = aux.calculate(surfacePointsPos);
            delete(aux);
            numSurfacePoints = size(surfacePointsPos, 1);
            numSourc = size(sourcePos, 1);
            
            originalSourcePos = sourcePos;
            monopoleSourcePos = surfacePointsPos;
            sourcePos = [originalSourcePos; monopoleSourcePos];
            
            originalSourceCoeff = sourceCoeff;
            monopoleSourceCoeff = U;
            
            sourceCoef = [originalSourceCoeff; monopoleSourceCoeff];
            
            originalSourceOrient = repmat([0 1 0 0], numSourc, 1);
            monopoleSourceOrient = repmat([0 1 0 0], numMeas, 1);
            sourceOrient = [originalSourceOrient; monopoleSourceOrient];
            
            originalSourceRadPat = repmat({@(x) simulator.monopoleRadPat(x)}, numSourc, 1);
            monopoleSourceRadPat = repmat({@(x) simulator.monopoleRadPat(x)}, numSurfacePoints, 1);
            sourceRadPat = [originalSourceRadPat; monopoleSourceRadPat];
            
            obj.setSources(sourcePos, 'coefficient', sourceCoef, 'orientation', sourceOrient, 'radiationPattern', sourceRadPat);
        end
        
    end
    
    methods(Static)
        
        function [monopoleCoeff, dipoleCoeff] = KirchhoffIntegralCoeff(k, sourcePos, sourceCoeff, surfacePointsPos, surfacePointsNormal, surfacePointsArea)
            % We will assume that the sources are monopoles
            
            % Calulate the value of the field and the gradient at the
            % surfacePoints
            % Field
            numSources = size(sourcePos, 1);
            radPat = repmat({@(x) simulator.monopoleRadPat(x)}, numSources, 1);
            objAux = simulator;
            objAux.setSources(sourcePos, 'coefficient', sourceCoeff, 'radiationPattern', radPat);
            U = objAux.calculate(surfacePointsPos);
            delete(objAux);
            
            % Directional derivative
            grad_U = sphericalWave_grad(k, sourceCoeff, sourcePos, surfacePointsPos);
            directDer_U = sum(grad_U.*surfacePointsNormal, 2);
            
            % Calculate the coeffients of monopoles and dipoles
            % (Kirchhoff Integral)
            monopoleCoeff = 1/(4*pi)*surfacePointsArea.*directDer_U;
            dipoleCoeff = -1/(4*pi)*surfacePointsArea.*U;
            
            
        end
        
        function radPat = monopoleRadPat(relDir, freq)
            % Input arguments:
            % - relDir. (numDir x 3) matrix
            % - freq. (numFreq x 3) matrix
            % Output arguments:
            % - radPat. (numDir x numFreq) matrix 
            
            if nargin == 1
                freq = 0;
            end
            
            numDir = size(relDir, 1);
            numFreq = numel(freq);
            
            radPat = ones(numDir, numFreq);
        end
        
        function radPat = dipoleRadPat(relDir, freq)
            % Input arguments:
            % - relDir. (numDir x 3) matrix
            % - freq. (numFreq x 3) matrix
            % Output arguments:
            % - radPat. (numDir x numFreq) matrix 
            
            if nargin == 1
                freq = 0;
            end
            
            numDir = size(relDir, 1);
            numFreq = numel(freq);
            
            % Broadside at Z axis
            broadside = [0 0 1];
            
            normRelDir = sqrt(sum(relDir.^2, 2));
            
            radPat = dot(repmat(broadside, [numDir, 1]), relDir, 2)./normRelDir;
            radPat = repmat(radPat, [1, numFreq]);
        end
        
        function radPat = dipolePreciseRadPat(relDir, k)
            % Broadside at Z axis
            broadside = [0 0 1];
            
            normRelDir = sqrt(sum(relDir.^2, 2));
            radPat = (-1i*k - 1./normRelDir).*dot(repmat(broadside, size(relDir, 1), 1), -relDir, 2)./normRelDir;
        end
        
        function quat = rotVec2quat(rotVec)
            % Rotation vector to quaternion
            % rotVec is a Mx4 matrix. Each row is a rotation vector
            % quat is a Mx4 matrix. Each row is a quaternion
            rotVecNorm = sqrt(sum(rotVec(:, 2:end).^2, 2));
            rotVec(:, 2:end) = rotVec(:, 2:end)./repmat(rotVecNorm, 1, 3);
            
            theta = -rotVec(:, 1);
            x = rotVec(:, 2);
            y = rotVec(:, 3);
            z = rotVec(:, 4);
            quat = [cos(theta/2), repmat(sin(theta/2), 1, 3).*[x, y, z]];
        end
        
        function rotVec = vec2rotVec(x)
            % x is a 3D vector that indicates where the Z axis has been
            % moved
            N = size(x, 1);
            origin = [0 0 1];
            normOrigin = sqrt(sum(origin.^2));
            normX = sqrt(sum(x.^2, 2));
            
            rotAxis = cross(repmat(origin, N, 1), x);
            singularity = all(rotAxis == 0, 2);
            rotAxis(singularity, :) = repmat([1 0 0], sum(singularity), 1);
            
            theta = acos(normOrigin*normX.*dot(repmat(origin, N, 1), x, 2));
            
            rotVec = [theta, rotAxis];
        end
        
        function vec = rotVec2BroadsideVec(x)
            % Axis Z is assumed to be the original broadside direction
            z = [0 0 1];
            
            quat = simulator.rotVec2quat(x);
            vec = quatrotate_custom(quat, z);
        end
        
        function acPath = calculateTheoricAcousticPaths(sourcePos, sourceRadPat, sourceOrient, recPos, recRadPat, recOrient, frequencies, propagVeloc)
            % Calculate the difference vectors between sources and measure
            % points
            % First dimension along recPositions, second dimension along
            % the X, Y and Z coordinates, third dimension along sources
            numSour = size(sourcePos, 1);
            numFreq = numel(frequencies);            
            numRec = size(recPos, 1);
            
            k = 2 * pi * frequencies / propagVeloc;
            
            % Quaterions for the rotation of the relative directions based
            % on orientations
            rotVec = sourceOrient;
            rotVec(:, 1) = -rotVec(:, 1); % Negative angle to compensate
            quat_source = simulator.rotVec2quat(rotVec);
            
            rotVec = recOrient;
            rotVec(:, 1) = -rotVec(:, 1); % Negative angle to compensate
            quat_rec = simulator.rotVec2quat(rotVec);
            
            acPath = zeros(numRec, numSour, numFreq);
            for s = 1:numSour
                sourceRadPatFun = sourceRadPat{s};

                diffVec = recPos - repmat(sourcePos(s, :), numRec, 1);
                dist = sqrt(sum(diffVec.^2, 2));

                % Quatrotate no va incluído en el matlab del laboratorio. relDir = quatrotate(quat, diffVec); Usar quatrotate_custom.
                relDir = quatrotate_custom( quat_source(s, :), diffVec );
                try
                sourceRadPatCoef = permute(sourceRadPatFun(relDir, frequencies), [1, 3, 2]);
                catch
                    4
                end
                
                % Receiver radiation pattern coefficent
                relDir = quatrotate_custom( quat_rec, -diffVec );
                recRadPatCoef = zeros(numRec, 1, numFreq);
                for r = 1:numRec
                    recRadPatFun = recRadPat{r};
                    recRadPatCoef(r, 1, :) = recRadPatFun(relDir(r, :), frequencies);
                end
                
                % Calculate
                aux = zeros(numRec, 1, numFreq);
                for f = 1:numFreq
                    aux(:, 1, f) = (sourceRadPatCoef(:, 1, f).*exp(-1i*k(f)*dist)./dist .* recRadPatCoef(:, 1, f));
                end
                acPath(:, s, :) = aux;
            end
            
            
        end
        
    end
    
end

