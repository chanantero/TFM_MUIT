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
        
        % Histogram
        visibleConfHistInd % Indices of cancellation configurations that are visible in the histogram graphs
    end
        
    % Public properties
    properties(Dependent)
        numExpScenarios
        numCancellations
        scenInd % Experimental scenario index
        confInd % Cancellation (WFS configuration) index      
    end
    
    % Private properties
    properties(Access = private)
        % 2D map
        scatNS
        scatWFS
        scatRec
        expScenInd % Experimental scenario index
        cancelInd % Cancellation (WFS configuration) index      
        
        % Histogram
        f
        axHistCoef
        axHistCancel
        mapVisibleConfToCoefHist
        mapVisibleConfToCancelHist
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
            obj.updateVisualization2D();
        end
        
        function set.cancelInd(obj, value)
            obj.cancelInd = value;
            obj.updateVisualization2D();
        end
        
        function set.ax2Dmap(obj, value)
            obj.ax2Dmap = value;
            obj.findScatObjects();
        end
        
        function set.expInfStruct(obj, value)
            obj.expInfStruct = value;
            obj.setIndices(1, 1);
            obj.generateVisualizationHistogram();
        end
        
        function scenInd = get.scenInd(obj)
            scenInd = obj.expScenInd;
        end
        
        function set.scenInd(obj, value)
            if value > 0 && value <= obj.numExpScenarios
                obj.setIndices(value, obj.cancelInd);
            else
                warning('simulationViewer:wrongIndex', 'The scenInd propertie must be an integer bigger than 0 and less or equal than numExpScenarios');
            end
        end
        
        function confInd = get.confInd(obj)
            confInd = obj.cancelInd;
        end
        
        function set.confInd(obj, value)
            if value > 0 && value <= obj.numCancellations
                obj.setIndices(obj.expScenInd, value);
            else
                warning('simulationViewer:wrongIndex', 'The confInd propertie must be an integer bigger than 0 and less or equal than numCancellations');
            end
        end
        
        function set.visibleConfHistInd(obj, value)
            obj.visibleConfHistInd = value;
            obj.updateVisibilityHistogram();
        end
        
    end
    
    methods
        function obj = simulationViewer(ax, expInfStruct)
            obj.f = figure;
            obj.axHistCoef = axes(obj.f);
            obj.axHistCancel = axes(obj.f);
            obj.axHistCoef.NextPlot = 'Add';
            obj.axHistCoef.OuterPosition = [0 0 0.5 1];
            obj.axHistCancel.NextPlot = 'Add';
            obj.axHistCancel.OuterPosition = [0.5 0 0.5 1];
            
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
        
        function updateVisualization2D(obj)
            
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
        
        function generateVisualizationHistogram(obj)
            simulationViewer.visualizeSignalCoefficients(obj.expInfStruct.NScoef, obj.expInfStruct.WFScoef, obj.expInfStruct.recOnlyNoiseCoef, obj.expInfStruct.recCoef, obj.axHistCoef, obj.axHistCancel);
            
            obj.mapVisibleConfToCoefHist.originInd = reshape(repmat(1:obj.numCancellations, 2, 1), [obj.numCancellations*2, 1]);
            obj.mapVisibleConfToCoefHist.destinationInd = (obj.numCancellations*2):-1:1;
            
            obj.mapVisibleConfToCancelHist.originInd = 1:obj.numCancellations;
            obj.mapVisibleConfToCancelHist.destinationInd = obj.numCancellations:-1:1;
            
            obj.visibleConfHistInd = true(obj.numCancellations, 1);
        end
        
        function updateVisibilityHistogram(obj)
            
            visibleCoefHist = false(obj.numCancellations*2+1, 1);
            visibleCoefHist(end) = true;
            visibleCoefHist(obj.mapVisibleConfToCoefHist.destinationInd) = obj.visibleConfHistInd(obj.mapVisibleConfToCoefHist.originInd);
                      
            for k = 1:numel(visibleCoefHist)
                if visibleCoefHist(k)
                    obj.axHistCoef.Children(k).Visible = 'on';
                else
                    obj.axHistCoef.Children(k).Visible = 'off';
                end
            end
            
            visibleCancelHist = false(obj.numCancellations, 1);
            visibleCancelHist(obj.mapVisibleConfToCancelHist.destinationInd) = obj.visibleConfHistInd(obj.mapVisibleConfToCancelHist.originInd);
                           
            for k = 1:numel(visibleCancelHist)
                if visibleCancelHist(k)
                    obj.axHistCancel.Children(k).Visible = 'on';
                else
                    obj.axHistCancel.Children(k).Visible = 'off';
                end
            end
        end
    end
   
    methods(Static)
        function visualizeSignalCoefficients(NScoef, WFScoef, recOnlyNoiseCoef, recCoef, axHistCoef, axHistCancel)
            % recOnlyNoiseCoef. P-element vector. Coefficients when there is only the effect of the noise source.
            % coef. (P x N) matrix. Coefficients after cancellation
            
            N = size(recCoef,2);
            
            % Get power of the sinusoidal signals
            powNS = (abs(NScoef).^2)/2;
            powWFS = (abs(WFScoef).^2)/2;
            powOnlyNoise = (abs(recOnlyNoiseCoef).^2)/2;
            powRec = (abs(recCoef).^2)/2;
            
            % Get Cancellation
            cancel = powRec./repmat(powOnlyNoise, 1, N);
            cancelGlobal = sum(powRec, 1)./sum(powOnlyNoise);
            
            % Convert to dB
            powNSLog = 10*log10(powNS);
            powWFSLog = 10*log10(powWFS);
            powOnlyNoiseLog = 10*log10(powOnlyNoise);
            powRecLog = 10*log10(powRec);
            cancelLog = 10*log10(cancel);
            
            % Get the strings with the statistical information
            strRecOnlyNoise = simulationViewer.getStatsInfoStr( powOnlyNoise );
            strWFS = simulationViewer.getStatsInfoStr( powWFS );
            strRec = simulationViewer.getStatsInfoStr( powRec );
            strCancel = simulationViewer.getStatsInfoStr( cancel );
                % Add the global cancellation information to the end of
                % each cancellation information string
            for n = 1:N
                strCancel{n} = [strCancel{n}, ', $C_{\mathit{global}} = ', num2str(cancelGlobal(n), 2), '$'];
            end
            
            % Clear axes
            cla(axHistCoef)
            cla(axHistCancel)
            
            % Delete legends
            legend(axHistCoef, 'off');
            legend(axHistCancel, 'off');
              
            % Logaritmic Histogram
            edgesCoef = -110:10:40; % Edges for the coefficients histogram
            edgesCancel = -110:10:20; % Edges for the cancellation histogram
            
            powOnlyNoiseLog(powOnlyNoiseLog < edgesCancel(1)) = edgesCancel(1);
            powWFSLog(powWFSLog < edgesCoef(1)) = edgesCoef(1);
            powRecLog(powRecLog < edgesCoef(1)) = edgesCoef(1);
            cancelLog(cancelLog < edgesCancel(1)) = edgesCancel(1);
            
            powOnlyNoiseLog(powOnlyNoiseLog > edgesCancel(end)) = edgesCancel(end);
            powWFSLog(powWFSLog > edgesCoef(end)) = edgesCoef(end);
            powRecLog(powRecLog > edgesCoef(end)) = edgesCoef(end);
            cancelLog(cancelLog > edgesCancel(end)) = edgesCancel(end);
            
            map = colormap('lines');
            
            histogram(axHistCoef, powOnlyNoiseLog, edgesCoef, 'Normalization', 'Probability', 'DisplayStyle', 'stairs', 'EdgeColor', 'black');
            for n = 1:N
                histogram(axHistCoef, powWFSLog(:, n), edgesCoef, 'Normalization', 'Probability', 'DisplayStyle', 'stairs', 'LineStyle', '--', 'EdgeColor', map(n, :));
                histogram(axHistCoef, powRecLog(:, n), edgesCoef, 'Normalization', 'Probability', 'DisplayStyle', 'stairs', 'EdgeColor', map(n, :));
                histogram(axHistCancel, cancelLog(:, n), edgesCancel, 'Normalization', 'Probability', 'DisplayStyle', 'stairs', 'EdgeColor', map(n, :));
            end
            
            axHistCoef.XTick = edgesCoef;
            axHistCoef.XTickLabels{1} = '-Inf';
            axHistCoef.XTickLabels{end} = 'Inf';
            axHistCoef.XLim = [edgesCoef(1), edgesCoef(end)];
            axHistCoef.XLabel.String = 'Normalized Power (dB)';
            legend(axHistCoef, [strRecOnlyNoise; strWFS; strRec], 'Interpreter', 'Latex');
            
            axHistCancel.XTick = edgesCancel;
            axHistCancel.XTickLabels{1} = '-Inf';
            axHistCancel.XTickLabels{end} = 'Inf';
            axHistCancel.XLim = [edgesCancel(1), edgesCancel(end)];
            axHistCancel.XLabel.String = 'Cancellation (dB)';
            legend(axHistCancel, strCancel, 'Interpreter', 'Latex');

        end
        
        function strs = getStatsInfoStr( coefMat )
            
            N = size(coefMat, 2);
            
            coefMean = mean(coefMat, 1);
            coefMax = max(coefMat, [], 1);
            coefMin = min(coefMat, [], 1);
            coefStd = std(coefMat, 0, 1);
            
            strs = cell(N, 1);
            for n = 1:N
                strs{n} = ['$\bar{p} = ', num2str(coefMean(n), 2), '$, $\max{p} = ', num2str(coefMax(n), 2), '$, $\min{p} = ', num2str(coefMin(n), 2), '$, $\sigma_p = ', num2str(coefStd(n), 2), '$'];
            end
            
        end
    end
    
end

