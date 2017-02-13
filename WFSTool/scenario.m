classdef scenario < handle
            
    properties(SetAccess = private)
        panel
        sourcePosition
        loudspeakersPosition
        activeLoudspeakers
        delays
        attenuations
    end
    
    methods
        
        function obj = scenario(parent)
            
            % Default
            sourcePosition = [0 1 0];
            loudspeakersPosition = [-1 0 0; 1 0 0]; 
            roomPosition = [-2, -2, 4, 4];
            
            obj.panel = obj.createGraphics(parent, @(ax) obj.mouseClickCallback(ax));
            obj.setScenario(sourcePosition, loudspeakersPosition, roomPosition);
            
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
        end
        
        function setScenario(obj, sourcePosition, loudspeakersPosition, roomPosition)
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
            obj.activeLoudspeakers = true(size(loudspeakersPosition, 1), 1);
            
            obj.updateDelaysAndAttenuations();
        end
             
    end
    
    methods(Access = private)
        
        function panel = createGraphics(obj, parent, callback)
            panel = uipanel(parent, 'BackgroundColor','white', 'Units', 'normalized', 'Position', [0.55 0.05 0.4 0.9]);
            
            ax = axes(panel, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.8],...
                'ButtonDownFcn', @(hObject, ~) callback(hObject), 'Box', 'on',...
                'XTick', [], 'XTickMode', 'manual', 'YTick', [], 'YTickMode', 'manual');
            
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
            axes = scat.Parent;
            point = axes.CurrentPoint;
            
            % Which is the closest loudspeaker to the current point?
            [~, ind] = min(sum((repmat(point(1,:), size(obj.loudspeakersPosition, 1), 1) - obj.loudspeakersPosition).^2, 2));
            obj.activeLoudspeakers(ind) = ~obj.activeLoudspeakers(ind);
            if obj.activeLoudspeakers(ind)
                scat.CData(ind, :) = [0 0 1];
            else
                scat.CData(ind, :) = [1 1 1];
            end
        end
        
        function updateDelaysAndAttenuations(obj)
            % Each row is a 3 element vector with the 3D position
            N = size(obj.loudspeakersPosition, 1);
            relPos = obj.loudspeakersPosition - repmat(obj.sourcePosition, N, 1);
            dist = sqrt(sum(relPos.^2, 2));
                        
            c = 340; % m/s
            obj.delays = dist/c;
            obj.attenuations = obj.calcAttenuation(dist);
        end
        
        function attenuation = calcAttenuation(obj, dist)
            attenuation = ones(size(dist));
            attenuation(~obj.activeLoudspeakers) = 0;
        end
    end
    
end

