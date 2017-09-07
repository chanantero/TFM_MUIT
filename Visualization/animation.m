classdef animation < handle
    % Given an array of doubles and the correspondent parameters, we can
    % define a method of visualization.
    % First, we define which of the parameters will be the independent
    % variables in the representation.
    % Then, we select which indices of the parameters we will use, in the
    % independent variables as well as in the rest.
    % When this is done we select and permute to have the independent
    % variables as the lower dimensions, and the rest as the higher
    % dimensions
    % Finally, just particularizing the non-independent parameters, we can
    % represent the data.
    
    properties
        % Essential variables for precalculated mode
        parameters
        dataArrays
        indepVecInd
        indepVecInd_last
        singularSelInd
        singularSelInd_last
        histSingularSelInd
        % Essential variables for real time mode
        functs
        funcNotIndepDims % Special function for when the are not independen variables
        indepVec
        indepVec_last
        singularSel
        singularSel_last
        histSingularSel % For multiple lines at the same time
        % Essential variables for both modes
%         indepDims % It is in the observable variables
        indepDims_last
        % Contextual information
        parameterLabels
        dataLabels
        % Derived variables that are convenient
        numDim
        sizeArray
        numDataArrays
        dimOrder
        % User Control Interface
        paramUIControl
        indepToggle
        modeFlag % False: precalculated mode. True: real time mode.
        modeFlag_last
        modeToggle
        % Representation
            % Precalculated mode {
        currParticSelInd % Current indices with respect to the paramters of the particularized array. Only the non independent matter
        histCurrParticSelInd % In case there are added plots
            % }
        discretized2D % Logical
        discretized2D_last
        discretized2DToggle % Uicontrol
        clearButton % Uicontrol
        addPlotButton % Uicontrol
        fig
        axs % Vector with the handles of the represented axes
        graphics % Vector of cells with the handles of the represented objects. Each cell is for each axes in axs
    end
    properties (SetObservable)
        numIndepDim % Determines if plot or surf
        indepDims
    end
    
    methods
        function obj = animation(parameters, dataArrays, parameterLabels, dataLabels, functions, funcNotIndepDims)
            if(nargin == 2)
                parameterLabels = cell(numel(parameters), 1);
                dataLabels = cell(numel(dataArrays), 1);
            end
            obj.parameters = parameters;
            obj.parameterLabels = parameterLabels;
            obj.dataArrays = dataArrays;
            obj.dataLabels = dataLabels;
            obj.functs = functions;
            obj.funcNotIndepDims = funcNotIndepDims;
            obj.numDim = numel(parameters);
            obj.sizeArray = cellfun(@(p) numel(p), parameters);
            obj.numDataArrays = numel(dataArrays);
            
%             obj.particParameters = parameters;
%             obj.particDataArrays = dataArrays;
%             obj.sizeParticArray = obj.sizeArray;
%             obj.particSelInd = cellfun(@(p) 1:numel(p), parameters, 'UniformOutput', false);
            obj.dimOrder = 1:obj.numDim;
%             
%             obj.candIndepDim = [];
%             obj.candCurrParticSelInd = ones(obj.numDim, 1);
            
            obj.numIndepDim = 0;
            obj.currParticSelInd = ones(obj.numDim, 1);
            obj.histCurrParticSelInd = [];
            obj.discretized2D = false;
            obj.modeFlag = false;
            
            obj.singularSelInd = ones(obj.numDim, 1);
            obj.indepVecInd = cellfun(@(p) 1:numel(p), parameters, 'UniformOutput', false);
            obj.indepVec = parameters;
            obj.singularSel = cellfun(@(p) p(1), parameters, 'UniformOutput', false);
            
            addlistener(obj, 'indepDims', 'PostSet',...
                @obj.callbackIndepDims);
            
            addlistener(obj, 'numIndepDim', 'PostSet',...
                @obj.callbackNumIndepDim);
            
            obj.createGUI();
        end
        
        function createGUI(obj)
            obj.fig = figure('Units', 'pixels', 'Position', [100 50 1200 600]);
            
            paramPanel = uipanel(obj.fig,'Title','Parameters', 'Units', 'normalized',...
                'Position',[.05 .05 .4 .9]);
            
            obj.paramUIControl = animation.createUIparamControl(paramPanel, obj.parameters, obj.parameterLabels, @(hObj, evntData) obj.readUserParamSelection());
            
            % Correction of the callback of the value editable texts,
            % because the function createUIparamControl doesn't allow so
            % many options
            for k = 1:obj.numDim
                if ~iscategorical(obj.parameters{k})
                    obj.paramUIControl{k}.currValueText.Callback = @(hObj, evndData) obj.singularSelEditedCallback(k);
                end
            end
                        
            p = uipanel(obj.fig, 'Units', 'normalized', 'Position',[.5 .05 .5 .9]);
            obj.axs = gobjects(obj.numDataArrays, 1);
            for k = 1:obj.numDataArrays
                pos = [0 (obj.numDataArrays - k)/obj.numDataArrays 1 1/obj.numDataArrays];
                obj.axs(k) = axes(p, 'Units', 'normalized', 'OuterPosition', pos);
                obj.axs(k).UIContextMenu = uicontextmenu;
                uimenu(obj.axs(k).UIContextMenu, 'Label', 'Export Axes', 'Callback', @(source, eventData) obj.exportSpecifiedAxes(k));
            end
            rotate3d(obj.fig); % Enable rotation
            
            obj.graphics = cell(obj.numDataArrays, 1);
            
            obj.discretized2DToggle  = uicontrol(obj.fig, 'Visible', 'off', 'Style', 'togglebutton', 'String', 'Dicretize 2D', 'Units', 'normalized', ...
                'Position', [0.5 0.95 0.1 0.05], 'Callback', @(hObject, eventData) obj.setDiscretized2D(hObject.Value));
            
            obj.clearButton  = uicontrol(obj.fig, 'Visible', 'off', 'Style', 'pushbutton', 'String', 'Clear', 'Units', 'normalized', ...
                'Position', [0.7 0.95 0.09 0.05], 'Callback', @(hObject, eventData) obj.clearCallback());
            
            obj.addPlotButton = uicontrol(obj.fig, 'Visible', 'off', 'Style', 'pushbutton', 'String', 'Add plot', 'Units', 'normalized', ...
                'Position', [0.8 0.95 0.09 0.05], 'Callback', @(hObject, eventData) obj.addPlot());
            
            obj.modeToggle = uicontrol(obj.fig, 'Visible', 'on', 'Style', 'togglebutton', 'String', 'Real Time', 'Units', 'normalized', ...
                'Position', [0.05 0.95 0.1 0.05], 'Callback', @(hObject, eventData) obj.setMode(hObject.Value));
            
        end
        
        function readUserParamSelection(obj)
            
            [values, independDims, numPoints, minims, maxims] = animation.getValuesFromParamUIControl(obj.paramUIControl);
            obj.indepDims = animation.setIndependentDimensions(obj.indepDims, independDims);
            iscat = cellfun(@(p) iscategorical(p), obj.parameters);
            if obj.modeFlag
                % Real time mode. Values are taken exactly.
                for k = 1:obj.numDim
                    if iscat(k)
                        obj.singularSel{k} = categorical(values(k), categories(obj.parameters{k}), 'Ordinal', true);
                    else
                        obj.singularSel{k} = values{k};
                    end
                end
                obj.setIndepVecs(independDims, numPoints, minims, maxims);
            else
                % Precalculated mode. Indices are approximated
                for k = 1:obj.numDim
                    if iscat(k)
                        obj.singularSelInd(k) = find(obj.parameters{k} == values{k}); 
                    else
                        obj.singularSelInd(k) = findIndBestApprox(obj.parameters{k}, values{k});
                    end
                end
                obj.setIndepVecInds(independDims, minims, maxims);
            end
            
                      
            % In case the new values need corrections, we update the UI
            obj.updateUserControl();
            
            obj.updateRepresentation();
        end
        
        function singularSelEditedCallback(obj, k)
            obj.paramUIControl{k}.slider.Value = str2double(obj.paramUIControl{k}.currValueText.String);
            obj.readUserParamSelection();
        end
        
        function updateUserControl(obj) % From variables to UI
            if obj.modeFlag
                obj.updateUserControl_RealTime();
            else
                obj.updateUserControl_PreCalc();
            end
        end
        
        function updateUserControl_PreCalc(obj)
            % From candidates chosen by user to graphical objects
            animation.setParamUIControl_basedOnIndices(obj.paramUIControl, obj.parameters, obj.indepVecInd, obj.singularSelInd, obj.indepDims);
        end
  
        function updateUserControl_RealTime(obj)                        
            % From candidates chosen by user to graphical objects
            animation.setParamUIControl_basedOnValues(obj.paramUIControl, obj.parameters, cell(obj.numDim, 1), cell(obj.numDim, 1), obj.indepVec, obj.singularSel, obj.indepDims);
        end               
        
        function [fig, ax] = exportSpecifiedAxes(obj, number, fig, saveLegend)
            % [fig, ax] = exportSpecifiedAxes(obj, number, fig, saveLegend)
            
            
            
            if nargin < 3
                fig = figure;
                inputFigOrAxes = true;
            else
                if strcmp(fig.Type, 'axes')
                    axesDest = fig;
                    fig = axesDest.Parent;
                    inputFigOrAxes = false;
                else
                    inputFigOrAxes = true;
                end
            end
            if nargin < 4
                saveLegend = false;
            end
            
            axes = obj.axs(number);
            
            if saveLegend
                legend = getappdata(axes,'LegendPeerHandle');
            else
                legend = [];
            end
            
            if inputFigOrAxes
                [obs] = copyobj([axes legend], fig);
                ax = obs(1);
                ax.Units = 'normalized';
                ax.OuterPosition = [0 0 1 1];
            else
                copyobj(axes.Children, axesDest);
                ax = axesDest;
            end
        end
        
        function setMode(obj, flag)    
            obj.modeFlag = flag;
            if flag
                obj.fromPreCalcToRealTime();
            else
                obj.fromRealTimeToPreCalc();
            end
            obj.updateUserControl();
            obj.updateRepresentation();
        end
        
        function fromPreCalcToRealTime(obj)
            for k = 1:obj.numDim
                obj.indepVec{k} = obj.parameters{k}(obj.indepVecInd{k});
                obj.singularSel{k} = obj.parameters{k}(obj.singularSelInd(k));
            end            
        end
        
        function fromRealTimeToPreCalc(obj)
            for k = 1:obj.numDim
                obj.indepVecInd{k} = unique(findIndBestApprox(double(obj.parameters{k}), double(obj.indepVec{k})));
                obj.singularSelInd(k) = findIndBestApprox(double(obj.parameters{k}), double(obj.singularSel{k}));
            end
        end
        
        function addPlot(obj)
            if(obj.numIndepDim == 1)
                for k = 1:obj.numDataArrays;
                    lastLine = obj.graphics{k}(end);
                    obj.graphics{k} = [obj.graphics{k}; gobjects(1, 1)];
                    obj.axs(k).NextPlot = 'add';
                    obj.graphics{k}(end) = plot(obj.axs(k), lastLine.XData, lastLine.YData);
                    obj.axs(k).NextPlot = 'replace';
                    obj.histSingularSel = [obj.histSingularSel, obj.singularSel];
                    obj.histSingularSelInd = [obj.histSingularSelInd, obj.singularSelInd];
                end
            end
        end

        function clearCallback(obj)
            for k = 1:obj.numDataArrays
                for l = 1:numel(obj.graphics{k})-1
                    delete(obj.graphics{k}(l));
                end
                obj.graphics{k} = obj.graphics{k}(end);
            end
        end
        
        function callbackIndepDims(obj, src, evnt)
            if(strcmp(src.Name, 'indepDims'))
                value = evnt.AffectedObject.indepDims(:);
                obj.numIndepDim = numel(value);
                obj.dimOrder = [value; find(~ismember((1:obj.numDim)', value))];
            end
        end
        
        function callbackNumIndepDim(obj, src, evnt)
            if(strcmp(src.Name, 'numIndepDim'))
                value = evnt.AffectedObject.numIndepDim;
                if(value == 2)
                    obj.discretized2DToggle.Visible = 'on';
                    obj.clearButton.Visible = 'off';
                    obj.addPlotButton.Visible = 'off';
                elseif(value == 1)
                    obj.discretized2DToggle.Visible = 'off';
                    obj.clearButton.Visible = 'on';
                    obj.addPlotButton.Visible = 'on';
                else
                    obj.discretized2DToggle.Visible = 'off';
                    obj.clearButton.Visible = 'off';
                    obj.addPlotButton.Visible = 'off';
                end
            end
        end
        
        function updateRepresentation(obj)
            
            if obj.numIndepDim == 0 && obj.modeFlag
                % Show the spectrum only
                spectrumData = obj.funcNotIndepDims(obj.singularSel{:});
                for k = 1:obj.numDataArrays
                    % Clear axes
                    cla(obj.axs(k), 'reset');
                    obj.graphics{k} = plot(obj.axs(k), spectrumData.investigated_theta, spectrumData.spectrum);
                    obj.axs(k).NextPlot = 'add';
                    for l = 1:numel(spectrumData.calcDoA)
                        plot(obj.axs(k), [spectrumData.calcDoA(l), spectrumData.calcDoA(l)], obj.axs(k).YLim, 'r--')
                    end
                    obj.axs(k).NextPlot = 'replace';
                end
                
                % Update "_last" variables
                obj.indepDims_last = obj.indepDims;
                obj.modeFlag_last = obj.modeFlag;
                obj.indepVec_last = obj.indepVec;
                obj.singularSel_last = obj.singularSel;
                obj.indepVecInd_last = obj.indepVecInd;
                obj.singularSelInd_last = obj.singularSelInd;
                obj.discretized2D_last = obj.discretized2D;
                return;
            end
            
            % Calculate the array to represent
            isIndep = false(obj.numDim, 1);
            isIndep(obj.indepDims) = true;
            if obj.modeFlag                
                reprParamValues = cell(obj.numDim, 1);
                reprParamValues(~isIndep) = obj.singularSel(~isIndep);
                reprParamValues(isIndep) = obj.indepVec(isIndep);
                valueArrays = animation.calculateValues(reprParamValues, obj.functs);
                       
            else
                reprParamValues = cell(obj.numDim, 1);
                reprParamInd = cell(obj.numDim, 1);
                for k = 1:obj.numDim
                    if isIndep(k)
                        reprParamInd{k} = obj.indepVecInd{k};
                        reprParamValues{k} = obj.parameters{k}(obj.indepVecInd{k});
                    else
                        reprParamInd{k} = obj.singularSelInd(k);
                        reprParamValues{k} = obj.parameters{k}(obj.singularSelInd(k));
                    end
                end
                
                valueArrays = cell(obj.numDataArrays, 1);
                for k = 1:obj.numDataArrays
                    valueArrays{k} = obj.dataArrays{k}(reprParamInd{:});
                end            
            end
            
            % Permute dimensions
            params = reprParamValues(obj.dimOrder);
            valueArrays = cellfun(@(array) permute(array, obj.dimOrder), valueArrays, 'UniformOutput', false);
            
            
            % Check if we need to reset the axes. This is necessary when:
            % - The independent dimensions have changed.
            % - The representation mode has changed (real time,
            % precalculated).
            % - Discretized2D has changed
            % - The represented independent vectors change (indepVec,
            % indepVecInd).
            % Reset the axes implies to delete all drawn objects and put
            % the correct axes labels.
            indepDimChange = ~isequaln(obj.indepDims_last, obj.indepDims);
            reprModeChange = obj.modeFlag_last ~= obj.modeFlag;
            if ~indepDimChange
                if obj.modeFlag
                    reprIndepVecChange = ~isequaln(obj.indepVec_last(obj.indepDims_last), obj.indepVec(obj.indepDims));
                else
                    reprIndepVecChange = ~isequaln(obj.indepVecInd_last(obj.indepDims_last), obj.indepVecInd(obj.indepDims));
                end
            else
                reprIndepVecChange = false;
            end
            discretized2Dchange = obj.discretized2D ~= obj.discretized2D_last;
            
            resetAxes = indepDimChange || reprModeChange || reprIndepVecChange ...
                || discretized2Dchange;
                
            if resetAxes
                for k = 1:obj.numDataArrays
                    % Clear axes
                    cla(obj.axs(k), 'reset');
                    
                    % Reset graphs vector
                    if obj.numIndepDim == 1     
                        obj.graphics{k} = plot(obj.axs(k), params{1}, zeros(numel(params{1}), 1));
                    elseif obj.numIndepDim == 2
                        if obj.discretized2D
                            obj.graphics{k} = gobjects(numel(params{2}), 1);
                            obj.axs(k).NextPlot = 'add';
                            for l = 1:numel(params{2})
                                obj.graphics{k}(l) = plot(obj.axs(k), params{1}, zeros(numel(params{1}), 1));
                            end
                            obj.axs(k).NextPlot = 'replace';
                            l = legend(obj.axs(k), toCellstr(params{2}));
                            l.Title.String = obj.parameterLabels{obj.indepDims(2)};
                        else
                            obj.graphics{k} = surf(obj.axs(k), double(params{1}), double(params{2}), zeros(numel(params{2}), numel(params{1})));
                        end
                    end
                    
                % Add axes labels
                if obj.numIndepDim == 1
                        obj.axs(k).XLabel.String = obj.parameterLabels{obj.indepDims(1)};
                        obj.axs(k).YLabel.String = obj.dataLabels{k};                    
                elseif obj.numIndepDim == 2
                    if obj.discretized2D
                        obj.axs(k).XLabel.String = obj.parameterLabels{obj.indepDims(1)};
                        obj.axs(k).YLabel.String = obj.dataLabels{k};
                    else
                        obj.axs(k).XLabel.String = obj.parameterLabels{obj.indepDims(1)};
                        obj.axs(k).YLabel.String = obj.parameterLabels{obj.indepDims(2)};
                        obj.axs(k).ZLabel.String = obj.dataLabels{k};
                    end
                end
                end
                
                
            end

            % The only thing that is left is to change de dependent dimension
            % data of the current graph.                               
            if obj.numIndepDim == 1
                for k = 1:obj.numDataArrays
                    obj.graphics{k}(end).YData = valueArrays{k};
                end
            elseif obj.numIndepDim == 2
                if obj.discretized2D
                    for k = 1:obj.numDataArrays
                        for l = 1:numel(params{2})
                            obj.graphics{k}(l).YData = valueArrays{k}(:,l);
                        end
                    end
                else
                    for k = 1:obj.numDataArrays
                        obj.graphics{k}.ZData = valueArrays{k}';
                    end
                end
            end
            
            % Update "_last" variables
            obj.indepDims_last = obj.indepDims;
            obj.modeFlag_last = obj.modeFlag;
            obj.indepVec_last = obj.indepVec;
            obj.singularSel_last = obj.singularSel;
            obj.indepVecInd_last = obj.indepVecInd;
            obj.singularSelInd_last = obj.singularSelInd;
            obj.discretized2D_last = obj.discretized2D;
            
        end
        
        function showLevelCurve(obj, levels, number, tolerance)
            % Create params
            isIndep = ismember(1:obj.numDim, obj.indepDims);
            params = cell(obj.numDim, 1);
            params(isIndep) = obj.indepVec(isIndep);
            params(~isIndep) = obj.singularSel(~isIndep);
            
            % Create dataArray
            permInd = zeros(obj.numDim, 1);
            permInd(obj.indepDims) = [1 2];
            permInd(~isIndep) = (3:obj.numDim);
            dataArray = permute(obj.graphics{number}.ZData.', permInd);
            
            f = figure; ax = axes(f, 'NextPlot', 'add');
            for k = 1:numel(levels)
                [x, y] = animation.calcLevelCurvePoints(params, dataArray, obj.functs{number}, obj.indepDims, levels(k), tolerance );
                scatter(ax, x, y)
            end
            legend(ax, numToCellstr(levels))
           
        end          
                         
        function setDiscretized2D(obj, value)
            if(obj.discretized2D ~= value)
                obj.discretized2D = value;
                obj.updateRepresentation();
            end
        end
        
        function setIndepVecs(obj, dimensions, numPoints, minims, maxims)
            for k = 1:numel(dimensions)
            vec = obj.indepVec{dimensions(k)};
            if ~iscategorical(vec)
                obj.indepVec{dimensions(k)} = linspace( minims(k), maxims(k), numPoints(k));
            end
            end
        end
        
        function setIndepVecInds(obj, dimensions, minims, maxims)
            for k = 1:numel(dimensions)
                vec = obj.parameters{dimensions(k)};
                if ~iscategorical(vec)
                    obj.indepVecInd{dimensions(k)} = find(vec >= minims(k) & vec <= maxims(k));
                end
            end
        end
            
    end
    
    methods(Static)
        function values = calculateValues(parameters, functions)
            
            numDim = numel(parameters);
            
            paramMats = cell(numDim, 1);
            [paramMats{:}] = ndgrid_multClass(parameters{:});
%             % Convert categoricals to double
%             parametersDouble = parameters;
%             iscat = find(cellfun(@(p) iscategorical(p), parameters));
%             for k = 1:numel(iscat)
%                 parametersDouble{iscat(k)} = double(parameters{iscat(k)}); % ngrid can't work with categorical arrays
%             end
%             % Calculate grids
%             [paramMats{:}] = ndgrid(parametersDouble{:});
%             % Restore the categoricals again
%             for k = 1:numel(iscat)
%                 categ = categories(parameters{iscat(k)});
%                 paramMats{iscat(k)} = categorical(paramMats{iscat(k)}, 1:numel(categ), categ); % Turn again to categorical arrays
%             end

            % New way
            values = funcResult2grid(paramMats, functions);
            % funcResult2grid returns a cell array. But, as we know
            % (assume) that the functions return scalars, we convert them
            % to double with cell2mat
            numResults = numel(functions);
            for k = 1:numResults
                values{k} = cell2mat(values{k}{1});
            end
            
            % Old way
%             % Convert every array to cell array
%             for k = 1:numDim
%                 paramMats{k} = num2cell(paramMats{k});
%             end
% 
%             aux = cell(numDim, 1);
%             for k = 1:numDim
%                 aux{k} = ones(sizeArray(k), 1);
%             end
%             % It is complicated, but it is like that
%             points = mat2cell(cat(numDim+1, paramMats{:}), aux{:}, numDim);
%             
%             values = cell(numResults, 1);
%             for k = 1:numResults
%                 values{k} = zeros(sizeArray);
%                 for l = 1:prod(sizeArray)
%                     % point = cellfun(@(p) p(l), paramMats, 'UniformOutput', false);
%                     point = points{l};
%                     values{k}(l) = functions{k}(point{:});
%                 end
%             end
            
        end
        
        function externalPosition = changePosCoord(containerPosition, internalPosition)
            externalPosition = [containerPosition(1) + internalPosition(1)*containerPosition(3), containerPosition(2) + internalPosition(2)*containerPosition(4), internalPosition(3)*containerPosition(3), internalPosition(4)*containerPosition(4)];
        end
        
        function paramUIControlStruct = createUIparamControl(parent, parameters, paramLabels, callbackFunc)
            
            % For each parameter, create its control bar
            numDims = numel(parameters);
            paramUIControlStruct = cell(numDims, 1);
            relPosSlider = [0.1 0.05 0.8 0.25];
            relPosLabel = [0.4 0.7 0.2 0.2];
            relPosCurrValue = [0.4 0.4 0.2 0.2];
            relPosMinim = [0.1 0.4 0.2 0.2];
            relPosMaxim = [0.7 0.4 0.2 0.2];
            relPosButtonGroup = [0.1 0.1 0.8 0.4];
            relPosIndepToggle = [0.8 0.6 0.1 0.1];
            
            for k = 1:numDims
                extPos = [0 (numDims - k)/numDims 1 1/numDims];
                
                param = parameters{k};
                if(isa(param, 'double'))
                    
                    maxim = max(param);
                    minim = min(param);
                    
                    s = struct();
                    s.slider = uicontrol(parent, 'Style', 'slider', 'Min', minim, 'Max', maxim, 'Value', minim,...
                        'SliderStep', [1/10; 1/5], 'Units', 'normalized', 'Position', animation.changePosCoord(extPos, relPosSlider), ...
                        'Callback', callbackFunc);
                    
                    s.paramLabel = uicontrol(parent, 'Style', 'text', 'String', paramLabels{k}, 'Units', 'normalized',...
                        'Position', animation.changePosCoord(extPos, relPosLabel));
                    
                    s.currValueText = uicontrol(parent, 'Style', 'edit', 'String', num2str(minim), 'Units', 'normalized',...
                        'Position', animation.changePosCoord(extPos, relPosCurrValue), ...
                        'Callback', callbackFunc);
                    
                    s.numPoints = uicontrol(parent, 'Style', 'edit', 'String', num2str(numel(parameters{k})), 'Units', 'normalized',...
                        'Position', animation.changePosCoord(extPos, relPosCurrValue), ...
                        'Callback', callbackFunc, 'Visible', 'off');
                    
                    s.minim = uicontrol(parent, 'Style', 'edit', 'String', num2str(min(parameters{k})), 'Units', 'normalized',...
                        'Position', animation.changePosCoord(extPos, relPosMinim), ...
                        'Callback', callbackFunc, 'Visible', 'off');
                    
                    s.maxim = uicontrol(parent, 'Style', 'edit', 'String', num2str(max(parameters{k})), 'Units', 'normalized',...
                        'Position', animation.changePosCoord(extPos, relPosMaxim), ...
                        'Callback', callbackFunc, 'Visible', 'off');
                                        
                elseif(isa(param, 'categorical'))
                    
                    s = struct();
                    s.bg = uibuttongroup(parent, 'Visible','off',...
                        'Units', 'normalized', 'Position', animation.changePosCoord(extPos, relPosButtonGroup),...
                        'SelectionChangedFcn', callbackFunc);
                    
                    categoryNames = categories(param);
                    numCat = numel(categoryNames);
                    
                    rb = gobjects(numCat, 1);
                    for l = 1:numCat
                        rb(l) = uicontrol(s.bg, 'Style', 'radiobutton',...
                            'String', categoryNames{l}, 'Units', 'normalized',...
                            'Position',[(l-1)/numCat 0 1/numCat 1],...
                            'HandleVisibility','off', 'Tag', num2str(l));
                    end
                    s.rb = rb;
                    
                    s.bg.Visible = 'on';
                                        
                end
                
                s.indepToggle = uicontrol(parent, 'Style', 'togglebutton', 'Value', 0,...
                    'Units', 'normalized', 'Position', animation.changePosCoord(extPos, relPosIndepToggle),...
                    'Callback', callbackFunc);
                
                paramUIControlStruct{k} = s;
            end
        end
        
        function [values, indepDims, numPoints, minims, maxims] = getValuesFromParamUIControl(paramUIControlVec)
            numDims = numel(paramUIControlVec);
            values = cell(numDims, 1); % Valores para double, indices para categorical
            indepDims = false(numDims, 1);
            numPoints = zeros(numDims, 1);
            minims = zeros(numDims, 1);
            maxims = zeros(numDims, 1);
            for k = 1:numDims
                if isfield(paramUIControlVec{k}, 'slider') % Variable tipo double
                    values{k} = paramUIControlVec{k}.slider.Value;
                    numPoints(k) = str2double(paramUIControlVec{k}.numPoints.String);
                    minims(k) = str2double(paramUIControlVec{k}.minim.String);
                    maxims(k) = str2double(paramUIControlVec{k}.maxim.String);
                elseif isfield(paramUIControlVec{k}, 'bg') % Variable tipo categorical
%                     values(k) = str2double(paramUIControlVec{k}.bg.SelectedObject.Tag);
                    values{k} = paramUIControlVec{k}.bg.SelectedObject.String;
                else
                    error('animation:getValuesFromParamUIControl', 'Unknown error')
                end
                indepDims(k) = paramUIControlVec{k}.indepToggle.Value;
            end
            indepDims = find(indepDims);
            numPoints = numPoints(indepDims);
            minims = minims(indepDims);
            maxims = maxims(indepDims);
        end
        
        function newIndepDims = setIndependentDimensions(currIndepDims, candIndepDims)
            numCand = numel(candIndepDims);
            numCurr = numel(currIndepDims);
            currIndepDims = currIndepDims(:);
            candIndepDims = candIndepDims(:);
            if isequaln(sort(currIndepDims), sort(candIndepDims))
                newIndepDims = currIndepDims;
            else
                if numel(candIndepDims) > 2
                    newIndepDims = currIndepDims;
                elseif numel(candIndepDims) == 2
                    if numCurr == 0
                        newIndepDims = candIndepDims;
                        warning('Cuidado, esto no debería ocurrir')
                    elseif numCurr == 1
                        if any(currIndepDims == candIndepDims)
                            newIndepDims = [currIndepDims; candIndepDims(candIndepDims ~= currIndepDims)];
                        else
                            newIndepDims = candIndepDims;
                            warning('Cuidado, esto no debería ocurrir')
                        end
                    elseif numCurr == 2
                        newIndepDims = candIndepDims;
                        warning('Cuidado, esto no debería ocurrir')
                    end
                elseif numCand == 1
                    newIndepDims = candIndepDims;
                elseif numCand == 0
                    newIndepDims = candIndepDims;
                else
                    error('Ni idea de lo que ocurre')
                end
            end
        end
        
        function setParamUIControl_basedOnIndices(paramUIControlStructure, params, indepVecInds, singularSelectedIndices, indepDimens)
            % For each parameter, update its control bar
            numDims = numel(params);
            for k = 1:numDims
                param = params{k};
                % Distinguish between double and categorical parameters
                if(isa(param, 'double'))
                    % Update slider bar
                    paramUIControlStructure{k}.slider.Value = param(singularSelectedIndices(k));
                    paramUIControlStructure{k}.slider.Min = min(param);
                    paramUIControlStructure{k}.slider.Max = max(param);
                    % Update editable text
                    paramUIControlStructure{k}.currValueText.String = num2str(paramUIControlStructure{k}.slider.Value);
                    paramUIControlStructure{k}.minim.String = min(param(indepVecInds{k}));
                    paramUIControlStructure{k}.maxim.String = max(param(indepVecInds{k}));
                elseif(isa(param, 'categorical'))
                    % Update radio button
                    paramUIControlStructure{k}.bg.SelectedObject = paramUIControlStructure{k}.rb(singularSelectedIndices(k));      
                end
                % Update independence related graphic objects
                paramUIControlStructure{k}.indepToggle.Value = any(k == indepDimens);
                if paramUIControlStructure{k}.indepToggle.Value
                    paramUIControlStructure{k}.slider.Visible = 'off';
                    paramUIControlStructure{k}.currValueText.Visible = 'off';
                    paramUIControlStructure{k}.minim.Visible = 'on';
                    paramUIControlStructure{k}.maxim.Visible = 'on';   
                else
                    paramUIControlStructure{k}.slider.Visible = 'on';
                    paramUIControlStructure{k}.currValueText.Visible = 'on';
                    paramUIControlStructure{k}.minim.Visible = 'off';
                    paramUIControlStructure{k}.maxim.Visible = 'off';
                end
                
                % Hide real time mode related graphic objects
                paramUIControlStructure{k}.numPoints.Visible = 'off';
            end
        end
        
        function setParamUIControl_basedOnValues(paramUIControlStructure, params, minimSel, maximSel, indepVecs, singularSelectedValues, indepDimens)
            % For each parameter, update its control bar
            numDims = numel(paramUIControlStructure);
            for k = 1:numDims
                indepVec = indepVecs{k};
                % Distinguish between double and categorical parameters
                if(isa(params{k}, 'double'))
                    % Update slider bar
                    paramUIControlStructure{k}.slider.Value = singularSelectedValues{k};
                    if ~isempty(minimSel{k})
                        paramUIControlStructure{k}.slider.Min = minimSel{k};
                    end
                    if ~isempty(maximSel{k})
                        paramUIControlStructure{k}.slider.Max = maximSel{k};
                    end
                    % Update editable text
                    paramUIControlStructure{k}.currValueText.String = num2str(paramUIControlStructure{k}.slider.Value);
                    paramUIControlStructure{k}.numPoints.String = numel(indepVec);
                    paramUIControlStructure{k}.minim.String = min(indepVec);
                    paramUIControlStructure{k}.maxim.String = max(indepVec);
                elseif(isa(indepVec, 'categorical'))
                    % Update radio button
                    selected = ismember(params{k}, singularSelectedValues{k});
                    paramUIControlStructure{k}.bg.SelectedObject = paramUIControlStructure{k}.rb(selected);
                end
                % Update independence related graphic objects
                paramUIControlStructure{k}.indepToggle.Value = any(k == indepDimens);
                if paramUIControlStructure{k}.indepToggle.Value
                    paramUIControlStructure{k}.slider.Visible = 'off';
                    paramUIControlStructure{k}.currValueText.Visible = 'off';
                    paramUIControlStructure{k}.numPoints.Visible = 'on';
                    paramUIControlStructure{k}.minim.Visible = 'on';
                    paramUIControlStructure{k}.maxim.Visible = 'on';
                else
                    paramUIControlStructure{k}.slider.Visible = 'on';
                    paramUIControlStructure{k}.currValueText.Visible = 'on';
                    paramUIControlStructure{k}.numPoints.Visible = 'off';
                    paramUIControlStructure{k}.minim.Visible = 'off';
                    paramUIControlStructure{k}.maxim.Visible = 'off';
                end
                
            end
        end
        
         function [x, y] = calcLevelCurvePoints(params, dataArray, func, indepDims, thresValue, tolerance )
        % intedDims must have 2 components
        % x and y corerspont to the first and second indepDim respectively
        
            % Threshols in dimension 1
            thresDim = indepDims(1);

            [ thresholds1, ~, ~, ~] = searchThreshold( params, dataArray, func, thresDim, thresValue, tolerance );
            
            % If threshold is empty, substitute it for max(params(thresDim))
            % if it is all above threshold or min(params(thresDim)) if it
            % is all under threshold.
%             thresholds1(allAboveThres1) = num2cell(max(params{thresDim})*ones(sum(allAboveThres1), 1));
%             thresholds1(allUnderThres1) = num2cell(min(params{thresDim})*ones(sum(allUnderThres1), 1));
            
            % Threshols in dimension 2
            thresDim = indepDims(2);

            [ thresholds2, ~, ~, ~] = searchThreshold( params, dataArray, func, thresDim, thresValue, tolerance );
            
            % If threshold is empty, substitute it for max(params(thresDim))
            % if it is all above threshold or min(params(thresDim)) if it
            % is all under threshold.
%             thresholds2(allAboveThres2) = num2cell(max(params{thresDim})*ones(sum(allAboveThres2), 1));
%             thresholds2(allUnderThres2) = num2cell(min(params{thresDim})*ones(sum(allUnderThres2), 1));
            
            
            % Calculate vectors
%             paramMats = cell(numel(params), 1);
%             paramsDouble = cellfun(@(p) double(p), params, 'UniformOutput', false);
%             
%             paramsDouble1 = paramsDouble;
%             paramsDouble1{indepDims(1)} = 0;
%             [paramMats{:}] = ndgrid(paramsDouble1{:});
%             [ threshVec, paramVec ] = cellContents2Vector( thresholds1, paramMats );
%             x1 = paramVec{indepDims(2)}; y1 = threshVec;
% 
%             paramsDouble2 = paramsDouble;
%             paramsDouble2{indepDims(2)} = 0;
%             [paramMats{:}] = ndgrid(paramsDouble2{:});
%             [ threshVec, paramVec ] = cellContents2Vector( thresholds2, paramMats );
%             x2 = paramVec{indepDims(1)}; y2 = threshVec;
% 
% The previous 3 lines are correct, but they are too general. As we know,
% there are only 2 dimensions with more than one element, then, let's do
% this. Basically, we take advantage of the fact that thresholds1 and
% thresholds 2 are vectors.
            [ threshVec, paramVec ] = cellContents2Vector( thresholds1, params(indepDims(2)) );
            x1 = paramVec{1}; y1 = threshVec;
            
            [ threshVec, paramVec ] = cellContents2Vector( thresholds2, params(indepDims(1)) );
            x2 = paramVec{1}; y2 = threshVec;
            
            % Unify them
            x = [y1; x2];
            y = [x1; y2];
            [x, ind] = sort(x);
            y = y(ind);
            
        end

        function [fig, ax] = exportAxes(axes, fig)
            if(nargin == 1)
                fig = figure;
            end
            ax = copyobj(axes, fig);
            ax.Units = 'normalized';
            ax.OuterPosition = [0 0 1 1];
        end
    end
    
end

