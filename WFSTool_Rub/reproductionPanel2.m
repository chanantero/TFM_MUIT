classdef reproductionPanel2 < handle
    
    properties
        panel
        list
        orderCallback
    end
    
    properties(SetAccess = private)
        selectedTrack
    end
    
    
    methods
        
        function obj = reproductionPanel2(parent, orderCallback)
            obj.panel = obj.createPanel(parent, @(buttonTag) obj.buttonCallback(buttonTag),...
                @(eventData) obj.changeSelection(eventData),...
                @() obj.setTrackCallback(),...
                @(eventData) obj.cellEditCallback(eventData));
            obj.list = findobj(obj.panel, 'Tag', 'list');
            obj.orderCallback = orderCallback;
        end
                         
        function trackName = getTrackName(obj, ind)
            trackNames = obj.getTrackNames();
            trackName = trackNames{ind};
        end
        
        function trackNames = getTrackNames(obj)
            N = size(obj.list.Data, 1); % Number of tracks/rows
            trackNames = cell(N, 1);
            for k = 1:N
                trackNames{k} = [obj.list.UserData.paths{k}, obj.list.Data{k, 1}];
            end
        end
        
        function setTrack(obj, FileName, PathName, index)
            obj.list.Data(index, 1) = {FileName};
            obj.list.UserData.paths(index) = {PathName};
        end
              
    end
    
    methods(Access = private)
        function panel = createPanel(~, parent, buttonCallback, listSelectionCallback, setTrackCallback, cellEditCallback)
            panel = uipanel(parent, 'BackgroundColor','white', 'Units', 'normalized', 'Position', [0.05, 0.5, 0.4, 0.4]);
            menu = uicontextmenu;
            list_ = uitable(panel, 'Units', 'Normalized', 'Position', [0.05, 0.2, 0.9, 0.7], ...
                'ColumnName',{'Signal', 'Real/Virtual'}, 'ColumnFormat', {'char', 'logical'}, 'Tag', 'list',...
                'ColumnEditable', [true, true], 'UIContextMenu', menu,...
                'CellSelectionCallback', @(hObject, eventData) listSelectionCallback(eventData),...
                'CellEditCallback', @(hObject, eventData) cellEditCallback(eventData));
            list_.UserData = struct('paths', []);
            
            uimenu(menu, 'Label', 'Select Audio Track', 'Callback', setTrackCallback());
            
            y0 = 0.05; h = 0.1; w = 0.1;
            playPos = [0.35, y0, w, h];
            stopPos = [0.2, y0, w, h];
            pausePos = [0.5, y0, w, h];
                        
            playButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', playPos, 'Callback', @(hObject, eventData) buttonCallback('play'));
            stopButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', stopPos, 'Callback', @(hObject, eventData) buttonCallback('stop'));
            pauseButton = uicontrol(panel, 'Style', 'pushbutton',...
                'Units', 'normalized', 'Position', pausePos, 'Callback', @(hObject, eventData) buttonCallback('pause'));
                        
            
            buttons = [playButton, stopButton, pauseButton];
            iconNames = {'images/playIcon.bmp', 'images/stopIcon.bmp', 'images/pauseIcon.jpg'};
            
            for k = 1:numel(buttons)
                button = buttons(k);
                imName = iconNames{k};
                
                % Get size in pixels
                button.Units = 'Pixels';
                width = button.Position(3);
                height = button.Position(4);
                % Set image as icon
                icon = imread(imName); % Read image
                icon = imresize_Rub(icon, [height, width]); % Resize
                % Convert to binary image
                if size(icon, 3) == 3
                    icon = rgb2gray(icon);
                end
                icon = ind2rgb(icon, [0 0 0; 1 1 1]);
                
                button.CData = icon;
                
                button.Units = 'Normalized';
            end
        end
        
        function buttonCallback(obj, buttonTag)
            % Create order structure
            switch buttonTag
                case 'play'
                    order.action = 'play';
                case 'stop'
                    order.action = 'stop';
                case 'pause'
                    order.action = 'pause';
            end
            
            % Send order structure to the callback
            obj.orderCallback(order);
        end
                
        function changeSelection(obj, eventdata)
            if ~isempty(eventdata.Indices)
                obj.selectedTrack = eventdata.Indices(1);
            else
                obj.selectedTrack = [];
            end
        end
        
        function setTrackCallback(obj)
            musicDirectory = 'C:\Users\Rubén\Music\Varias\';
            [FileName, PathName, ~] = uigetfile({'*.mp3', 'MP3 File'; '*.wav', 'WAV File'}, 'Select audio track', 'MultiSelect', 'on', musicDirectory);
                        
            if FileName ~= 0
                obj.setTrack(FileName, PathName, obj.selectedTrack);
            end
        end
        
        function cellEditCallback(~, eventData)
            
            if eventData.Indices(2) == 1
            % If the name of the track has been changed, check:
%             % - Is it a valid formulation for a complex number?
%             a = regexp(eventData.EditData, '(?<Amplitude>(\d+\.\d+|\d+))\*exp\(j\*(?<Phase>(\d+\.\d+|\d+))\)', 'once');
%             if isempty(a)
%                 error('reproductionPanel:cellEditCallback', 'Wrong format for a complex number');
%             end

            % - Is it in the valid format to represent a complex number?
            a = regexp(eventData.EditData, 'A:(?<Amplitude>(\d+\.\d+|\d+)) Ph:(?<Phase>(\d+\.\d+|\d+)) f:(?<Frequency>(\d+\.\d+|\d+))', 'once');
            if isempty(a)
                error('reproductionPanel:cellEditCallback', 'Wrong format for a complex number');
            end
            
            end
                        
        end
        
        
    end
   
    
end

