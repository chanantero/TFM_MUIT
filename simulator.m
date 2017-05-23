classdef simulator < handle
        
    properties
        ax
        sourcePositions
        sourceCoefficients
        k % Propagation constant
        
        % Graphical
        XLim
        YLim
        XnumPoints
        YnumPoints
        z
    end
    
    methods
        function obj = simulator(parentFig)
            obj.ax = axes(parentFig, 'CLim', [0 10]);
            obj.sourcePositions = [0 0 0];
            obj.sourceCoefficients = 1;
            obj.k = 1;
            obj.XLim = [-2 2];
            obj.YLim = [-2 2];
            obj.XnumPoints = 20;
            obj.YnumPoints = 20;
            obj.z = 0;
        end
        
        function simulate(obj)
            % Create measure points
            xVec = linspace(obj.XLim(1), obj.XLim(2), obj.XnumPoints)';
            yVec = linspace(obj.YLim(1), obj.YLim(2), obj.YnumPoints)';
            % rImage = [kron(xVec, ones(YnumPoints, 1)), repmat(yVec, XnumPoints, 1), zeros(XnumPoints*YnumPoints, 1)];
            [Y, X] = ndgrid(yVec, xVec);
            rImage = [X(:), Y(:), obj.z*ones(obj.XnumPoints*obj.YnumPoints, 1)];
            
            % Simulate
            U = obj.calculate(rImage);
            
            % Reshape as an image
            U = reshape(U, obj.YnumPoints, obj.XnumPoints);
            
%             U(real(U) > 10) = 10;
            
            image(obj.ax, 'XData', obj.XLim, 'YData', obj.YLim, 'CData', abs(U), 'CDataMapping', 'scaled')
        end
        
        function U = calculate(obj, measurePointPositions)
            U = sum(sphericalWave(obj.k, obj.sourceCoefficients, obj.sourcePositions, measurePointPositions), 2);
        end
    end
    
end

