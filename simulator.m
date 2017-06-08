classdef simulator < handle
        
    properties
        ax
        sourcePositions
        sourceCoefficients
        sourceOrientations % Mx4 matrix. Rotation vector: [angle of rotation, Xaxis, Yaxis, Zaxis]
        radPatFuns
        k % Propagation constant
        
        % Graphical
        XLim
        YLim
        XnumPoints
        YnumPoints
        z
        measurePoints
        
        % Results of simulation
        field
        imag
    end
    
    properties(Dependent)
        numSources
    end
    
    % Getters and setters
    methods
        function numSources = get.numSources(obj)
            numSources = size(obj.sourcePositions, 1);
        end
    end
    
    methods
        function obj = simulator(parentFig)
            if nargin == 1
                obj.ax = axes(parentFig, 'CLim', [-1 1]);
                colormap(obj.ax, 'gray');
            end
            obj.setDefaultProperties();
        end
        
        function setDefaultProperties(obj)
            sourcePos = [0 0 0];
            sourceCoeff = 1;
            theta = 0; rotAxis = [1 0 0];
            sourceOrient = [theta, rotAxis];
            radPatFun = {@(~) 1};
            obj.setSources(sourcePos, 'coefficient', sourceCoeff, 'orientation', sourceOrient,...
                'radiationPattern', radPatFun);
            
            obj.k = 1;
            obj.XLim = [-2 2];
            obj.YLim = [-2 2];
            obj.XnumPoints = 20;
            obj.YnumPoints = 20;
            obj.z = 0;
        end
        
        function simulate(obj)
            obj.generateMeasurePoints();
            
            % Simulate
            obj.field = obj.calculate(obj.measurePoints);
                     
            obj.draw();
        end
        
        function generateMeasurePoints(obj)
            % Create measure points
            obj.measurePoints= plane(obj.XLim, obj.YLim, obj.XnumPoints, obj.YnumPoints, [], [], []);
        end
        
        function draw(obj)
            % Reshape as an image
            U = reshape(obj.field, obj.XnumPoints, obj.YnumPoints).';
                       
            cla(obj.ax)
            obj.imag = image(obj.ax, 'XData', obj.XLim, 'YData', obj.YLim,...
                'CData', real(U), 'CDataMapping', 'scaled');
            
%             % Draw sources
%             obj.ax.NextPlot = 'Add';
%             scatter(obj.ax, obj.sourcePositions(:, 1), obj.sourcePositions(:, 2), 50, [1 1 1],...
%                 'MarkerEdgeColor', [0 0 0]);
%             obj.ax.NextPlot = 'Replace';
        end
       
        function animate(obj, f, t)
            incrT = diff(t);
            incrT = [incrT, incrT(end)];
            
            for l = 1:numel(t)
                U = obj.field*exp(1i*2*pi*f*t(l));
                U = reshape(U, obj.XnumPoints, obj.YnumPoints).';
                obj.imag.CData = real(U);
                pause(incrT(l))
            end
        end
        
        function U = calculate(obj, measurePointPositions)
            % Calculate the difference vectors between sources and measure
            % points
            % First dimension along measurePoints, second dimension along
            % the X, Y and Z coordinates, third dimension along sources
            numSour = obj.numSources;
            numMeasPoints = size(measurePointPositions, 1);

            U_aux = zeros(numMeasPoints, numSour);
            for s = 1:numSour
                sourceCoef = obj.sourceCoefficients(s);
                sourcePos = obj.sourcePositions(s, :);
                sourceOrient = obj.sourceOrientations(s, :);
                radPatFun = obj.radPatFuns{s};
                
                diffVec = measurePointPositions - repmat(sourcePos, numMeasPoints, 1);
                dist = sqrt(sum(diffVec.^2, 2));
                
                rotVec = sourceOrient;
                rotVec(1) = -rotVec(1);
                quat = simulator.rotVec2quat(rotVec);
                relDir = quatrotate(quat, diffVec);
                radPatCoef = radPatFun(relDir);
                
                U_aux(:, s) = sourceCoef*radPatCoef.*exp(-1i*obj.k*dist)./dist;
            end
            U = sum(U_aux, 2);
                       
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
            dipoleCoeff = 1/(4*pi)*surfacePointsArea.*U;
            
            
        end
        
        function radPat = monopoleRadPat(relDir)
            % relDir is a Mx3 matrix
            numPoints = size(relDir, 1);
            radPat = ones(numPoints, 1);
        end
        
        function radPat = dipoleRadPat(relDir)
            % Broadside at Z axis
            broadside = [0 0 1];
            
            normRelDir = sqrt(sum(relDir.^2, 2));
            
            radPat = dot(repmat(broadside, size(relDir, 1), 1), relDir, 2)./normRelDir;
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
            rotVecNorm = sqrt(sum(rotVec(2:end).^2, 2));
            rotVec(2:end) = rotVec(2:end)./repmat(rotVecNorm, 1, 3);
            
            theta = rotVec(:, 1);
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
        
    end
    
end

