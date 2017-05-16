classdef scenario < handle
            
    properties(SetAccess = private)
        panel
        sourcesPosition
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
            
        end         
        
        function setNumSources(obj, numSources)
            
            obj.sourcesPosition = zeros(numSources, 3);
            
            source = findobj(obj.panel, 'Tag', 'source');

            source.XData = zeros(1, numSources);
            source.YData = zeros(1, numSources);
            source.ZData = zeros(1, numSources);
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
            
            loudspeakers.XData = loudspeakersPosition(:, 1);
            loudspeakers.YData = loudspeakersPosition(:, 2);
            loudspeakers.ZData = loudspeakersPosition(:, 3);
            loudspeakers.CData = repmat([0 0 1], size(loudspeakersPosition, 1), 1);
            
            % Other variable
            obj.sourcesPosition = sourcesPosition;
            obj.loudspeakersPosition = loudspeakersPosition;
            obj.loudspeakersOrientation = loudspeakersOrientation;
            obj.loudspeakersState = true(size(loudspeakersPosition, 1), 1);
            obj.forcedDisabledLoudspeakers = false(size(loudspeakersPosition, 1), 1);
            obj.forcedEnabledLoudspeakers = false(size(loudspeakersPosition, 1), 1);
            
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
            
            source = scatter(ax, 0, 0, 'b');
            source.Tag = 'source';
            
            loudspeakers = scatter(ax, 0, 1, 'k', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', 'flat');
            loudspeakers.Tag = 'loudspeakers';
            loudspeakers.ButtonDownFcn = @(hObj, ~) obj.loudspeakersCallback(hObj);
            
            ax.XLim = [-1 1];
            ax.YLim = [-1 1];
            
        end
        
        function mouseClickCallback(obj, axes, index)
            x = axes.CurrentPoint(1, 1);
            y = axes.CurrentPoint(1, 2);
            
            source = findobj(axes, 'Tag', 'source');
            source.XData = x;
            source.YData = y;
            
            obj.sourcesPosition(index, 1) = x;
            obj.sourcesPosition(index, 2) = y;
                                 
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
        
        function updateLoudspeakersColor(obj)
            N = numel(obj.loudspeakersState);
            
            activeState = 2*ones(N, 1);
            activeState(obj.loudspeakersState) = 1;
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

            % Calculate active loudspeakers
            % Coseno del ángulo entre la línea que une la fuente virtual 
            % con cada altavoz y la dirección principal en la que el
            % altavoz emite ("broadside")
            cosAlfa = obj.calcCosAlfa();
            obj.loudspeakersState = cosAlfa > 0; % Buscamos los altavoces que debería estar activos
            
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
            atten = scenario.calcPhysicalAttenuation(distances);
            
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
        
        function ind = getEnabledLoudspeakers(obj)
            ind = (obj.loudspeakersState | obj.forcedEnabledLoudspeakers) & ~obj.forcedDisabledLoudspeakers;
        end
        
        
    end
    
    methods(Static)
        function attenuation = calcPhysicalAttenuation(distances)
            attenuation = ones(size(distances));
        end
    end
    
end

