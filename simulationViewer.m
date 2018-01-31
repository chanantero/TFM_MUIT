classdef simulationViewer < handle
    
    % Public properties
    properties
        % Results
        expInfStruct
        % Experiment information structure. Fields:
        % - NScoef
        % - WFScoef
        % - recNScoef
        % - recWFSCoef
        % - recCoef
        % - recNSCoefExp
        % - recWFSCoefExp
        % - recCoefExp
        
        % 2D map
        ax2Dmap
        representationType = categorical({'dB'}, simulationViewer.represTypeCat, simulationViewer.represTypeCat, 'Protected', true);
        magnitude = categorical({'Field'}, simulationViewer.magnitudeCategories, simulationViewer.magnitudeCategories, 'Protected', true);
        axHist
        hist
        
        % Histogram
        visibleConfHistInd % Indices of cancellation configurations that are visible in the histogram graphs
    end
    
    % Public properties
    properties(Dependent)
        numExpScenarios
        scenInd % Experimental scenario index
    end
    
    properties(Dependent, Access=private)
        WFScoef
        recNScoef
        recCoef
        numWFS
    end
    
    % Private properties
    properties(Access = private)
        % 2D map
        scatNS
        scatWFS
        scatRec
        expScenInd % Experimental scenario index
        
        % Histogram
        f
        axHistCoef
        axHistCancel
        mapVisibleConfToCoefHist
        mapVisibleConfToCancelHist
    end
    
    properties(Access = private, Constant)
        % Auxiliar purposes during construction of object
        magnitudeCategories = {'WFS', 'NS', 'Field', 'Cancellation', 'WFS2NSratio', 'WFS_exp', 'NS_exp', 'Cancellation_exp', 'Field_exp'}
        represTypeCat = {'dB', 'power', 'abs', 'phase', 'real', 'imag'};
    end
    
    % Getters and setters
    methods
        function numExpScen = get.numExpScenarios(obj)
            numExpScen = numel(obj.expInfStruct);
        end
        
        function set.ax2Dmap(obj, value)
            obj.ax2Dmap = value;
            obj.findScatObjects();
        end
        
        function set.expInfStruct(obj, value)
            obj.expInfStruct = value;
            obj.setIndices(1);
%             obj.generateVisualizationHistogram();
        end
        
        function scenInd = get.scenInd(obj)
            scenInd = obj.expScenInd;
        end
        
        function set.scenInd(obj, value)
            if value > 0 && value <= obj.numExpScenarios
                obj.setIndices(value);
            else
                warning('simulationViewer:wrongIndex', 'The scenInd propertie must be an integer bigger than 0 and less or equal than numExpScenarios');
            end
        end
        
        function set.visibleConfHistInd(obj, value)
            obj.visibleConfHistInd = value;
            obj.updateVisibilityHistogram();
        end
        
        function numWFS = get.numWFS(obj)
            CData = obj.scatWFS.CData;
            
            if isvector(CData)
                numWFS = numel(obj.scatWFS.CData);
            else
                numWFS = size(CData, 1);
            end
        end
        
        function set.representationType(obj, value)
            obj.representationType(1) = value; % Indexing is necessary when working with categorical arrays
            obj.updateVisualization2D();
        end
        
        function set.magnitude(obj, value)
            obj.magnitude(1) = value; % Indexing is necessary when working with categorical arrays
            obj.updateVisualization2D();
        end
        
        function WFScoef = get.WFScoef(obj)
            WFScoef = [obj.expInfStruct(:).WFScoef];
        end
        
        function recNScoef = get.recNScoef(obj)
            recNScoef = [obj.expInfStruct(:).recNScoef];
        end
        
        function recCoef = get.recCoef(obj)
            recCoef = [obj.expInfStruct(:).recCoef];
        end
    end
    
    % Public methods
    methods
        function obj = simulationViewer(ax2D, expInfStruct, histogramAxes)
            
%             obj.f = figure;
%             obj.axHistCoef = axes(obj.f);
%             obj.axHistCancel = axes(obj.f);
%             obj.axHistCoef.NextPlot = 'Add';
%             obj.axHistCoef.OuterPosition = [0 0 0.5 1];
%             obj.axHistCancel.NextPlot = 'Add';
%             obj.axHistCancel.OuterPosition = [0.5 0 0.5 1];
            
            obj.ax2Dmap = ax2D;
            obj.ax2Dmap.CLim = [-100, 20];
            
            if exist('histogramAxes', 'var') == 1
                obj.axHist = histogramAxes;
                obj.hist = histogram(obj.axHist, [], [0 1], 'Normalization', 'Probability');
            end
            
            obj.expInfStruct = expInfStruct;
            
            obj.magnitude = categorical({'WFS'}, obj.magnitudeCategories, obj.magnitudeCategories, 'Protected', true);
            obj.representationType = categorical({'dB'}, obj.represTypeCat, obj.represTypeCat, 'Protected', true);
            
        end
        
        function adjustColorbar(obj)
            minCLim = min(obj.scatRec.CData);
            maxCLim = max(obj.scatRec.CData);
            obj.ax2Dmap.CLim = [minCLim, maxCLim];
        end
    end
    
    % Private methods
    methods(Access = private)
        function setIndices(obj, expScenInd)
            obj.expScenInd = expScenInd;
            
            % Update positions
            obj.scatWFS.XData = obj.expInfStruct(obj.expScenInd).WFSposition(:, 1);
            obj.scatWFS.YData = obj.expInfStruct(obj.expScenInd).WFSposition(:, 2);
            obj.scatNS.XData = obj.expInfStruct(obj.expScenInd).NSposition(:, 1);
            obj.scatNS.YData = obj.expInfStruct(obj.expScenInd).NSposition(:, 2);
            obj.scatRec.XData = obj.expInfStruct(obj.expScenInd).recPosition(:, 1);
            obj.scatRec.YData = obj.expInfStruct(obj.expScenInd).recPosition(:, 2);
            
            obj.updateVisualization2D();
        end
        
        function findScatObjects(obj)
            obj.scatRec = findobj(obj.ax2Dmap, 'Tag', 'receiver');
            obj.scatWFS = findobj(obj.ax2Dmap, 'Tag', 'loudspeakers');
            obj.scatNS = findobj(obj.ax2Dmap, 'Tag', 'source');
        end
        
        function updateVisualization2D(obj)
            
            sScen = obj.expInfStruct(obj.expScenInd);
            
            WFScoeff = sScen.WFScoef;
            NScoef = sScen.NScoef;
            
            powWFS = (abs(WFScoeff).^2)/2;
            powNS = (abs(NScoef).^2)/2;
            
            powWFSDB = 10*log10(powWFS);
            powNSDB = 10*log10(powNS);
            
            % Set color
            % % Example of RGB:
            % recColorMap = colormap('gray');
            % indices = scaled2indexedColors(size(recColorMap, 1), [-100 20], powRecDB);
            % scatRec.CData = recColorMap(indices, :);
            obj.scatWFS.CData = powWFSDB;
            obj.scatNS.CData = powNSDB;
            
            switch obj.magnitude
                case 'WFS'
                    repCoef = sScen.recWFScoef;
                case 'NS'
                    repCoef = sScen.recNScoef;
                case 'Field'
                    repCoef = sScen.recCoef;
                case 'Cancellation'
                    repCoef = sScen.recCoef./sScen.recNScoef;
                case 'WFS2NSratio'
                    repCoef = sScen.recWFScoef./sScen.recNScoef;
                case 'WFS_exp'
                    repCoef = sScen.recWFScoefExp;
                case 'NS_exp'
                    repCoef = sScen.recNScoefExp;
                case 'Field_exp'
                    repCoef = sScen.recCoefExp;
                case 'Cancellation_exp'
                    repCoef = sScen.recCoefExp./sScen.recNScoefExp;
            end
            
            switch obj.representationType
                case 'dB'
                    powWFS = (abs(WFScoeff).^2)/2;
                    powNS = (abs(NScoef).^2)/2;
                    
                    xWFS = 10*log10(powWFS);
                    xNS = 10*log10(powNS);
                    
                    pow = (abs(repCoef).^2)/2;
                    x = 10*log10(pow);
                    
                    histTitle = 'dB';
                    edges = -110:10:40;
                case 'power'
                    xWFS = (abs(WFScoeff).^2)/2;
                    xNS = (abs(NScoef).^2)/2;
                    
                    x = (abs(repCoef).^2)/2;
                    
                    histTitle = 'Power';
                    edges = linspace(0, max(x), 10);
                case 'abs'
                    xWFS = abs(WFScoeff);
                    xNS = abs(NScoef);
                    
                    x = abs(repCoef);
                    
                    histTitle = 'Magnitude';
                    edges = linspace(0, max(x), 10);
                case 'phase'
                    xWFS = rad2deg(angle(WFScoeff));
                    xNS = rad2deg(angle(NScoef));
                    
                    x = rad2deg(angle(repCoef));
                    
                    histTitle = 'Phase';
                    edges = -180:10:180;
                case 'real'
                    xWFS = real(WFScoeff);
                    xNS = real(NScoef);
                    
                    x = real(repCoef);
                    
                    histTitle = 'Real';
                    edges = linspace(min(x), max(x), 10);
                case 'imag'
                    xWFS = imag(WFScoeff);
                    xNS = imag(NScoef);
                    
                    x = imag(repCoef);
                    
                    histTitle = 'Imaginary';
                    edges = linspace(min(x), max(x), 10);
            end
            
            obj.scatWFS.CData = xWFS;
            obj.scatNS.CData = xNS;
            obj.scatRec.CData = x;
            obj.adjustColorbar();
            
            if isgraphics(obj.hist)
                obj.hist.Data = x;
                obj.hist.BinEdges = edges;
                obj.axHist.Title.String = histTitle;
            end
            
        end
        
        function generateVisualizationHistogram(obj)
            simulationViewer.visualizeSignalCoefficients(obj.expInfStruct(1).NScoef, obj.WFScoef, obj.recNScoef(:, 1), obj.recCoef, obj.axHistCoef, obj.axHistCancel);
            
            obj.mapVisibleConfToCoefHist.originInd = reshape(repmat(1:obj.numExpScenarios, 2, 1), [obj.numExpScenarios*2, 1]);
            obj.mapVisibleConfToCoefHist.destinationInd = (obj.numExpScenarios*2):-1:1;
            
            obj.mapVisibleConfToCancelHist.originInd = 1:obj.numExpScenarios;
            obj.mapVisibleConfToCancelHist.destinationInd = obj.numExpScenarios:-1:1;
            
            obj.visibleConfHistInd = true(obj.numExpScenarios, 1);
        end
        
        function updateVisibilityHistogram(obj)
            
            visibleCoefHist = false(obj.numExpScenarios*2+1, 1);
            visibleCoefHist(end) = true;
            visibleCoefHist(obj.mapVisibleConfToCoefHist.destinationInd) = obj.visibleConfHistInd(obj.mapVisibleConfToCoefHist.originInd);
            
            for k = 1:numel(visibleCoefHist)
                if visibleCoefHist(k)
                    obj.axHistCoef.Children(k).Visible = 'on';
                else
                    obj.axHistCoef.Children(k).Visible = 'off';
                end
            end
            
            visibleCancelHist = false(obj.numExpScenarios, 1);
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
    
    % Static methods
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

