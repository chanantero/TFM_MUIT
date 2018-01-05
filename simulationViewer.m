classdef simulationViewer < handle
    
    % Public properties
    properties
        % Results
        expInfStruct
        % Experiment information structure. Fields:
        % - NScoef.
        % - recOnlyNoiseCoef.
        % - WFScoef.
        % - recCoef.
        
        % 2D map
        ax2Dmap
        scatNS
        scatWFS
        scatRec
        
        expScenInd % Experimental scenario index
        cancelInd % Cancellation (WFS configuration) index
    end
    
    % Public properties
    properties(Dependent)
        numExpScenarios
        numCancellations
    end
    
    % Getters and setters
    methods
        function numExpScen = get.numExpScenarios(obj)
            numExpScen = numel(obj.expInfStruct);
        end
        
        function numExpScen = get.numCancellations(obj)
            numExpScen = size(obj.expInfStruct.WFScoef, 2);
        end
        
        function set.expScenInd(obj, value)
            obj.expScenInd = value;
            obj.changeVisualization2D();
        end
        
        function set.cancelInd(obj, value)
            obj.cancelInd = value;
            obj.changeVisualization2D();
        end
        
        function set.ax2Dmap(obj, value)
            obj.ax2Dmap = value;
            obj.findScatObjects();
        end
        
        function set.expInfStruct(obj, value)
            obj.expInfStruct = value;
            obj.setIndices(1, 1);
        end
    end
    
    methods
        function obj = simulationViewer(ax, expInfStruct)
            obj.ax2Dmap = ax;
            obj.expInfStruct = expInfStruct;
        end
        
        function setIndices(obj, expScenInd, cancelInd)
            obj.expScenInd = expScenInd;
            obj.cancelInd = cancelInd;
        end
        
        function findScatObjects(obj)
            obj.scatRec = findobj(obj.ax2Dmap, 'Tag', 'receiver');
            obj.scatWFS = findobj(obj.ax2Dmap, 'Tag', 'loudspeakers');
            obj.scatNS = findobj(obj.ax2Dmap, 'Tag', 'source');
        end
        
        function changeVisualization2D(obj)
            
            recCoef = obj.expInfStruct(obj.expScenInd).recCoef(:, obj.cancelInd);
            WFScoef = obj.expInfStruct(obj.expScenInd).WFScoef(:, obj.cancelInd);
            NScoef = obj.expInfStruct(obj.expScenInd).NScoef;
                 
            powRec = (abs(recCoef).^2)/2;
            powWFS = (abs(WFScoef).^2)/2;
            powNS = (abs(NScoef).^2)/2;
            
            powRecDB = 10*log10(powRec);
            powWFSDB = 10*log10(powWFS);
            powNSDB = 10*log10(powNS);
                        
            % % Set RGB color
            % recColorMap = colormap('gray');
            % indices = scaled2indexedColors(size(recColorMap, 1), [-100 20], powRecDB);
            % scatRec.CData = recColorMap(indices, :);
            %
            % WFSColorMap = colormap('gray');
            % indices = scaled2indexedColors(size(WFSColorMap, 1), [-100 20], powWFSDB);
            % scatWFS.CData = WFSColorMap(indices, :);
            
            % Set color and
            obj.ax2Dmap.CLim = [-100, 20];
            obj.scatRec.CData = powRecDB;
            obj.scatWFS.CData = powWFSDB;
            obj.scatNS.CData = powNSDB;
            
        end
    end
    
     
end

