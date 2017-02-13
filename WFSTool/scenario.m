classdef scenario < handle
            
    properties(SetAccess = private)
        panel
        sourcePosition
        loudspeakersPosition
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
            
            % Other variable
            obj.sourcePosition = sourcePosition;
            obj.loudspeakersPosition = loudspeakersPosition;
            
            obj.updateDelaysAndAttenuations();
        end
             
    end
    
    methods(Access = private)
        
        function panel = createGraphics(~, parent, callback)
            panel = uipanel(parent, 'BackgroundColor','white', 'Units', 'normalized', 'Position', [0.55 0.05 0.4 0.9]);
            
            ax = axes(panel, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.8],...
                'ButtonDownFcn', @(hObject, ~) callback(hObject), 'Box', 'on',...
                'XTick', [], 'XTickMode', 'manual', 'YTick', [], 'YTickMode', 'manual');
            
            ax.NextPlot = 'Add';
            
            source = scatter(ax, 0, 0, 'b');
            source.Tag = 'source';
            
            loudspeakers = scatter(ax, 0, 1, 'k');
            loudspeakers.Tag = 'loudspeakers';
            
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
        
        function updateDelaysAndAttenuations(obj)
            % Each row is a 3 element vector with the 3D position
            N = size(obj.loudspeakersPosition, 1);
            relPos = obj.loudspeakersPosition - repmat(obj.sourcePosition, N, 1);
            dist = sqrt(sum(relPos.^2, 2));
                        
            c = 340; % m/s
            obj.delays = dist/c;
            obj.attenuations = obj.calcAttenuation(dist);
        end
        
    end
    
    methods(Static)
        function attenuation = calcAttenuation(dist)
            attenuation = zeros(size(dist));
        end
    end
    
end

