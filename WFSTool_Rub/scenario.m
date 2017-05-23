classdef scenario < handle
            
    properties(SetAccess = private)
        panel
        
        sourcesPosition
        activeSource
        
        loudspeakersPosition
        loudspeakersOrientation
        loudspeakersState % Logical. True is enabled, false is disabled
        
        forcedEnabledLoudspeakers
        forcedDisabledLoudspeakers % Priority over enabled loudspeakers
        
        delays
        attenuations
    end
    
    properties(Dependent)
        numSources
        numLoudspeakers
    end
    
    properties (Constant)
        c = 343;
        sourceColor = [0 0 1];
    end
    
    % Getters and setters
    methods
        function numSources = get.numSources(obj)
            numSources = size(obj.sourcesPosition, 1);
        end
        
        function numLoudspeakers = get.numLoudspeakers(obj)
            numLoudspeakers = size(obj.loudspeakersPosition, 1);
        end
    end
    
    methods
        
        function obj = scenario(parent)
            
            % Default
            sourcesPosition = [0 1 0];
            loudspeakersPosition = [-0.1 0 0; 0.1 0 0];
            loudspeakersOrientation = [1 0 0; -1 0 0];
            roomPosition = [-2, -2, 4, 4];
            
            obj.panel = obj.createGraphics(parent, @(ax) obj.mouseClickCallback(ax));
            obj.setScenario(sourcesPosition, loudspeakersPosition, loudspeakersOrientation, roomPosition);
            obj.activeSource = 1;
        end         
        
        function setNumSources(obj, numSources)
            
            obj.sourcesPosition = zeros(numSources, 3);
            if numSources > 0
                obj.activeSource = 1;
            else
                obj.activeSource = [];
            end
            obj.loudspeakersState = true(obj.numLoudspeakers, obj.numSources);

            source = findobj(obj.panel, 'Tag', 'source');

            source.XData = zeros(1, numSources);
            source.YData = zeros(1, numSources);
            source.ZData = zeros(1, numSources);
            source.CData = repmat(obj.sourceColor, numSources, 1);
            
        end
        
        function setSourcePosition(obj, sourcePosition, index)
            % Graphics
            source = findobj(obj.panel, 'Tag', 'source');
                 
            source.XData(index) = sourcePosition(1);
            source.YData(index) = sourcePosition(2);
            source.ZData(index) = sourcePosition(3);
            
            % Other variables
            obj.sourcesPosition(index, :) = sourcePosition;
            
            obj.updateDelaysAndAttenuations();

        end
        
        function setScenario(obj, sourcesPosition, loudspeakersPosition, loudspeakersOrientation, roomPosition)
            % Graphics
            ax = findobj(obj.panel, 'Type', 'Axes');
            source = findobj(obj.panel, 'Tag', 'source');
            loudspeakers = findobj(ax, 'Tag', 'loudspeakers');
            
            ax.XLim = [roomPosition(1), roomPosition(1)+roomPosition(3)];
            ax.YLim = [roomPosition(2), roomPosition(2)+roomPosition(4)];
            
            source.XData = sourcesPosition(:, 1);
            source.YData = sourcesPosition(:, 2);
            source.ZData = sourcesPosition(:, 3);
            source.CData = repmat(obj.sourceColor, size(sourcesPosition, 1), 1);
            
            loudspeakers.XData = loudspeakersPosition(:, 1);
            loudspeakers.YData = loudspeakersPosition(:, 2);
            loudspeakers.ZData = loudspeakersPosition(:, 3);
            loudspeakers.CData = repmat([0 0 1], size(loudspeakersPosition, 1), 1);
            
            % Other variable
            obj.sourcesPosition = sourcesPosition;
            if obj.numSources > 0
                obj.activeSource = 1;
            else
                obj.activeSource = [];
            end
            obj.loudspeakersPosition = loudspeakersPosition;
            obj.loudspeakersOrientation = loudspeakersOrientation;
            obj.loudspeakersState = true(obj.numLoudspeakers, obj.numSources);
            obj.forcedDisabledLoudspeakers = false(obj.numLoudspeakers, 1);
            obj.forcedEnabledLoudspeakers = false(obj.numLoudspeakers, 1);
            
            obj.updateDelaysAndAttenuations();
        end
            
        function setForcedEnabledLoudspeakers(obj, enabled)
            obj.forcedEnabledLoudspeakers = enabled(:);
        end
        
        function setForcedDisabledLoudspeakers(obj, enabled)
            obj.forcedDisabledLoudspeakers = enabled(:);
        end
    end
    
    methods(Access = private)
        
        function panel = createGraphics(obj, parent, callback)
            panel = uipanel(parent, 'BackgroundColor','white', 'Units', 'normalized', 'Position', [0.55 0.05 0.4 0.9]);
            
            ax = axes(panel, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.8],...
                'ButtonDownFcn', @(hObject, ~) callback(hObject), 'Box', 'on',...
                'XTick', [], 'XTickMode', 'manual', 'YTick', [], 'YTickMode', 'manual',...
                'DataAspectRatio', [1 1 1], 'DataAspectRatioMode', 'manual');
            
            ax.NextPlot = 'Add';
            
            source = scatter(ax, 0, 0, 50, obj.sourceColor);
            source.Tag = 'source';
            source.ButtonDownFcn = @(hObj, ~) obj.sourcesCallback(hObj);
            
            loudspeakers = scatter(ax, 0, 1, 'k', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', 'flat');
            loudspeakers.Tag = 'loudspeakers';
            loudspeakers.ButtonDownFcn = @(hObj, ~) obj.loudspeakersCallback(hObj);
            
            ax.XLim = [-1 1];
            ax.YLim = [-1 1];
            
        end
        
        function mouseClickCallback(obj, axes)
            x = axes.CurrentPoint(1, 1);
            y = axes.CurrentPoint(1, 2);
            
            source = findobj(axes, 'Tag', 'source');
            source.XData(obj.activeSource) = x;
            source.YData(obj.activeSource) = y;
            
            obj.sourcesPosition(obj.activeSource, 1) = x;
            obj.sourcesPosition(obj.activeSource, 2) = y;
                                 
            obj.updateDelaysAndAttenuations();
        end
        
        function loudspeakersCallback(obj, scat)
            point = scat.Parent.CurrentPoint;
            
            % Which is the closest loudspeaker to the current point?
            [~, ind] = min(sum((repmat(point(1,:), size(obj.loudspeakersPosition, 1), 1) - obj.loudspeakersPosition).^2, 2));
            obj.forcedDisabledLoudspeakers(ind) = ~obj.forcedDisabledLoudspeakers(ind);
            obj.updateDelaysAndAttenuations();

            obj.updateLoudspeakersColor();

        end
        
        function sourcesCallback(obj, scat)
            point = scat.Parent.CurrentPoint;
            
            % Which is the closest loudspeaker to the current point?
            [~, ind] = min(sum((repmat(point(1,:), size(obj.sourcesPosition, 1), 1) - obj.sourcesPosition).^2, 2));
            obj.activeSource = ind;
        end
        
        function updateLoudspeakersColor(obj)            
            activeState = 2*ones(obj.numLoudspeakers, 1);
            activeState(obj.loudspeakersState(:, obj.activeSource)) = 1;
            activeState(obj.forcedEnabledLoudspeakers) = 3;
            activeState(obj.forcedDisabledLoudspeakers) = 4;
            
            colors = [0 0 1; % Active
                      0 0 0.25; % Inactive
                      0 1 0; % Forced Enabled
                      1 0 0]; % Forced Disabled
                              
            scat = findobj(obj.panel, 'Tag', 'loudspeakers');
            scat.CData = colors(activeState, :);
        end
        
        function updateDelaysAndAttenuations(obj)
                      
            % Calculate distances between sources and loudspeakers
            dist = obj.calcDistances();
            obj.delays = dist/obj.c;          
            
            % Calculate attenuations based on those distances and
            % orientations
            obj.attenuations = obj.calcAttenuations(dist);
            
            % Update graphics            
            obj.updateLoudspeakersColor();
            
        end
        
        function dist = calcDistances(obj, sourceIndices)
            if nargin == 1
                sourceIndices = 1:obj.numSources;
            end
            
            dist = zeros(obj.numLoudspeakers, numel(sourceIndices));
            for k = 1:numel(sourceIndices)
                relPos = obj.loudspeakersPosition - repmat(obj.sourcesPosition(sourceIndices(k), :), obj.numLoudspeakers, 1);
                dist(:, k) = sqrt(sum(relPos.^2, 2));
            end
            
        end
        
        function atten = calcAttenuations(obj, distances)
            
            % Coseno del ángulo entre la línea que une la fuente virtual 
            % con cada altavoz y la dirección principal en la que el
            % altavoz emite ("broadside")
            cosAlfa = obj.calcCosAlfa();
            
            % Calculate active loudspeakers
%             obj.loudspeakersState = cosAlfa > 0; % Buscamos los altavoces que deberían estar activos
            obj.loudspeakersState(:) = true;

            atten = obj.calcPhysicalAttenuation(distances, cosAlfa);
                      
            atten(~obj.getEnabledLoudspeakers()) = 0;
            
            
        end
        
        function cosAlfa = calcCosAlfa(obj, sourceIndices)
            if nargin == 1
                sourceIndices = 1:obj.numSources;
            end
            
            % Coseno del ángulo entre la línea que une la fuente virtual 
            % con cada altavoz y la dirección principal en la que el
            % altavoz emite ("broadside")
            
            normOrient = modVec(obj.loudspeakersOrientation); % Ideally 1 if the orientations are normalized

            cosAlfa = zeros(obj.numLoudspeakers, numel(sourceIndices));
            for k = 1:numel(sourceIndices)
                relPos = obj.loudspeakersPosition - repmat(obj.sourcesPosition(sourceIndices(k), :), obj.numLoudspeakers, 1);
                dist = sqrt(sum(relPos.^2, 2));
                cosAlfa(:, k) = dot(relPos, obj.loudspeakersOrientation, 2)./(dist.*normOrient);
            end
        end
        
        function attenuation = calcPhysicalAttenuation(obj, distances, cosAlfa)
            if obj.numLoudspeakers == 96
                r0 = 1.44*(0.5+cosd(45));
                A = sqrt(r0./(r0+distances));
                attenuation = A.*cosAlfa./sqrt(distances);
            else
%                 attenuation = ones(size(distances));
                attenuation = min(1, (0.1)./distances);
            end
        end
        
        function ind = getEnabledLoudspeakers(obj, sourceIndices)
            if nargin == 1
                sourceIndices = 1:obj.numSources;
            end
            
            forcedEnabled = repmat(obj.forcedEnabledLoudspeakers, 1, numel(sourceIndices));
            forcedDisabled = repmat(obj.forcedDisabledLoudspeakers, 1, numel(sourceIndices));
            ind = (obj.loudspeakersState | forcedEnabled) & ~forcedDisabled;
        end
        
        
    end
    
end

