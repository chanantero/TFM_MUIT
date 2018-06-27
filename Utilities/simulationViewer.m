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
        
        % GUI
        fig
        ax2Dmap
        representationType = categorical({'dB'}, simulationViewer.represTypeCat, simulationViewer.represTypeCat, 'Protected', true);
        magnitude = categorical({'Field'}, simulationViewer.magnitudeCategories, simulationViewer.magnitudeCategories, 'Protected', true);
        axHist
        hist
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
        function obj = simulationViewer(ax2D, expInfStruct)
            
            function scenIndCallback(ind)
                obj.scenInd = ind;
            end
            
            function magnitudeCallback(magnitude)
                obj.magnitude = magnitude;
            end
            
            function reprTypeCallback(representationType)
                obj.representationType = representationType;
            end
                        
            [fig, ax2Dmap, axHist, hist] = simulationViewer.createGUI(ax2D, @(ind) scenIndCallback(ind), @(mag) magnitudeCallback(mag), @(reprType) reprTypeCallback(reprType));
            
            obj.fig = fig;
            obj.ax2Dmap = ax2Dmap;
            
            obj.axHist = axHist;
            obj.hist = hist;
  
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

            % Set color
            % % Example of RGB:
            % recColorMap = colormap('gray');
            % indices = scaled2indexedColors(size(recColorMap, 1), [-100 20], powRecDB);
            % scatRec.CData = recColorMap(indices, :);

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
                    
                    pow = (abs(repCoef).^2);
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
        
    end
    
    % Static methods
    methods(Static)
         
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
        
        function [fig, ax2Dmap, axHist, hist, buttonPannel] = createGUI(ax2DmapOrig, scenIndCallback, magnitudeCallback, represTypeCallback)
            fig = figure;
            
            ax2Dmap = copyobj(ax2DmapOrig, fig);
            ax2Dmap.Units = 'Normalized';
            ax2Dmap.OuterPosition = [0.5, 0, 0.5, 0.8];
            ax2Dmap.CLim = [-100, 20];
            colormap(ax2Dmap, 'jet')
            colorbar(ax2Dmap)
            
            % Create histogram axes
            axHist = axes(fig, 'Units', 'normalized', 'OuterPosition', [0 0 0.5 1]);
            hist = histogram(axHist, [], [0 1], 'Normalization', 'Probability');
            
            % Create buttons
            buttonPannel = uipanel(fig, 'Units', 'Normalized', 'Position', [0.52 0.82 0.46 0.16]);
            
            scenIndButton = uicontrol(buttonPannel, 'Style', 'Edit',...
                'Units', 'normalized', 'Position', [0 0 0.3 0.5],...
                'String', '1', 'Callback', @(hObject, eventData, handles) scenIndCallback(scenIndProcessing()));
            
            magnitudeButton = uicontrol(buttonPannel, 'Style', 'popupmenu',...
                'Units', 'normalized', 'Position', [0.4 0 0.3 0.5],...
                'String', simulationViewer.magnitudeCategories, 'Callback', @(hObject, eventData, handles) magnitudeCallback(magnitudeProcessing()));
            
            represTypeButton = uicontrol(buttonPannel, 'Style', 'popupmenu',...
                'Units', 'normalized', 'Position', [0.7 0 0.3 0.5],...
                'String', simulationViewer.represTypeCat, 'Callback', @(hObject, eventData, handles) represTypeCallback(represTypeProcessing()));
            
            function magnitude = magnitudeProcessing()
                magnitude = simulationViewer.magnitudeCategories(magnitudeButton.Value);
            end
            
            function representationType = represTypeProcessing()
                representationType = simulationViewer.represTypeCat(represTypeButton.Value);
            end
            
            function scenInd = scenIndProcessing()
                scenInd = str2double(scenIndButton.String);
                if isnan(scenInd) || floor(scenInd) ~= scenInd
                    warning('simulationViewer:scenIndButtonCallback', 'You must insert an integer number');
                    scenInd = [];
                end
            end
        end
        
    end
    
end

