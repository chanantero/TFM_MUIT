classdef scenario < handle
            
    properties(SetAccess = private)
        panel
        sourcePosition
        loudspeakersPosition
        loudspeakersOrientation
        loudspeakersState % Logical. True is enabled, false is disabled
        forcedEnabledLoudspeakers
        forcedDisabledLoudspeakers % Priority over enabled loudspeakers
        delays
        attenuations
    end
    
    properties (Constant)
        c = 343;
    end
    
    methods
        
        function obj = scenario(parent)
            
            % Default
            sourcePosition = [0 1 0];
            loudspeakersPosition = [-0.1 0 0; 0.1 0 0];
            loudspeakersOrientation = [1 0 0; -1 0 0];
            roomPosition = [-2, -2, 4, 4];
            
            obj.panel = obj.createGraphics(parent, @(ax) obj.mouseClickCallback(ax));
            obj.setScenario(sourcePosition, loudspeakersPosition, loudspeakersOrientation, roomPosition);
            
        end         
        
        function delays = getDelays(obj)
            delays = obj.delays;
        end
        
        function attenuations = getAttenuations(obj)
            attenuations = obj.attenuations;
        end
        
        function setSourcePosition(obj, sourcePosition)
            % Graphics
            source = findobj(obj.panel, 'Tag', 'source');
                 
            source.XData = sourcePosition(1);
            source.YData = sourcePosition(2);
            source.ZData = sourcePosition(3);
            
            % Other variables
            obj.sourcePosition = sourcePosition;
            
            obj.updateDelaysAndAttenuations();

        end
        
        function setScenario(obj, sourcePosition, loudspeakersPosition, loudspeakersOrientation, roomPosition)
            % Graphics
            ax = findobj(obj.panel, 'Type', 'Axes');
            source = findobj(obj.panel, 'Tag', 'source');
            loudspeakers = findobj(ax, 'Tag', 'loudspeakers');
            
            ax.XLim = [roomPosition(1), roomPosition(1)+roomPosition(3)];
            ax.YLim = [roomPosition(2), roomPosition(2)+roomPosition(4)];
            
            source.XData = sourcePosition(1);
            source.YData = sourcePosition(2);
            source.ZData = sourcePosition(3);
            
            loudspeakers.XData = loudspeakersPosition(:, 1);
            loudspeakers.YData = loudspeakersPosition(:, 2);
            loudspeakers.ZData = loudspeakersPosition(:, 3);
            loudspeakers.CData = repmat([0 0 1], size(loudspeakersPosition, 1), 1);
            
            % Other variable
            obj.sourcePosition = sourcePosition;
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
        
        function mouseClickCallback(obj, axes)
            x = axes.CurrentPoint(1, 1);
            y = axes.CurrentPoint(1, 2);
            
            source = findobj(axes, 'Tag', 'source');
            source.XData = x;
            source.YData = y;
            
            obj.sourcePosition(1) = x;
            obj.sourcePosition(2) = y;
                                 
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
            activeState(obj.loudspeakersState | obj.forcedEnabledLoudspeakers) = 3;
            activeState(obj.forcedDisabledLoudspeakers) = 1;
            
            colors = [1 1 1; 0 0 0.5; 0 0 1];
            
            scat = findobj(obj.panel, 'Tag', 'loudspeakers');
            scat.CData = colors(activeState, :);
        end
        
        function updateDelaysAndAttenuations(obj)
            % Each row is a 3 element vector with the 3D position
            N = size(obj.loudspeakersPosition, 1);
            relPos = obj.loudspeakersPosition - repmat(obj.sourcePosition, N, 1);
            dist = sqrt(sum(relPos.^2, 2));
                        
            % Calculate distances
            obj.delays = dist/obj.c;
            
            % Calculate active loudspeakers
            normOrient = modVec(obj.loudspeakersOrientation); % Ideally 1 if the orientations are normalized
            cosAlfa = dot(relPos, obj.loudspeakersOrientation, 2)./(dist.*normOrient); % Coseno del ángulo entre la línea que une la fuente virtual con cada altavoz y
            % la dirección principal en la que el altavoz emite ("broadside")
            obj.loudspeakersState = cosAlfa > 0; % Buscamos los altavoces que debería estar activos
            obj.updateLoudspeakersColor();
            
            % Calculate attenuations
            obj.attenuations = obj.calcAttenuation(dist);
            
        end
        
        function ind = getEnabledLoudspeakers(obj)
            ind = (obj.loudspeakersState | obj.forcedEnabledLoudspeakers) & ~obj.forcedDisabledLoudspeakers;
        end
        
        function attenuation = calcAttenuation(obj, dist)
            attenuation = ones(size(dist));
            attenuation(~obj.getEnabledLoudspeakers()) = 0;
        end
    end
    
end

