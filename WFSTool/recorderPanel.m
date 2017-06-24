classdef recorderPanel < handle
    
    properties(SetAccess = private)
        activeChannels
    end
    
    properties(SetAccess = private)
        panel
        list
    end
    
    events
        updatedValues
    end
    
    methods
        function obj = recorderPanel(parent, position)
            obj.panel = obj.createPanel(parent, position,...
                @(eventData) obj.cellEditCallback(eventData),...
                @(x) obj.setNumActiveChannelsCallback(x));
            obj.list = findobj(obj.panel, 'Tag', 'list');
            
            obj.setNumActiveChannelsCallback(0);
        end
        
        function disableGUI(obj)
            obj.list.Enable = 'off';
            numSourcesPopUp = findobj(obj.panel, 'Tag', 'numSources');
            numSourcesPopUp.Enable = 'off';
        end
        
        function enableGUI(obj)
            obj.list.Enable = 'on';
            numSourcesPopUp = findobj(obj.panel, 'Tag', 'numSources');
            numSourcesPopUp.Enable = 'on';
        end
        
    end
    
    methods(Access = private)
        
        function panel = createPanel(~, parent, position, cellEditCallback, setNumChannelsCallback)
            panel = uipanel(parent, 'BackgroundColor','white', 'Units', 'normalized', 'Position', position);
            menu = uicontextmenu;
            list_ = uitable(panel, 'Units', 'Normalized', 'Position', [0.05, 0.2, 0.9, 0.7], ...
                'ColumnName', {'Active Recorder Channel'}, 'ColumnFormat', {'numeric'}, 'Tag', 'list',...
                'ColumnEditable', true, 'UIContextMenu', menu,...
                'CellEditCallback', @(hObject, eventData) cellEditCallback(eventData));
            list_.UserData = struct('paths', []);
            
            % Number of sources uicontrol
            uicontrol(panel, 'Style', 'text', 'String', 'Nº active recorder channels', 'Units', 'normalized',...
                'Position', [0 0.9 0.2 0.1]);
            uicontrol(panel, 'Style', 'edit', 'String', '0', 'Units', 'normalized', 'Position', [0.2 0.9 0.1 0.1],...
                'Tag', 'numSources', 'Callback', @(hObject, ~, ~) setNumChannelsCallback(str2double(hObject.String)));
        end
                      
        function setNumActiveChannelsCallback(obj, numChannels)
            obj.activeChannels = (1:numChannels)';
            
            obj.list.Data = obj.activeChannels;
            evntData = updatedValuesEvntData('activeRecorderChannels');
            notify(obj, 'updatedValues', evntData);
        end
        
        function cellEditCallback(obj, eventData)
            
            row = eventData.Indices(1);
            column = eventData.Indices(2);
            
            newChannel = eventData.EditData;
            
            if column == 1 % There is no other possibility
                chan = str2double(newChannel);
                
                if isnan(chan)
                    obj.list.Data(row, column) = eventData.PreviousData;
                else
                    candidate = obj.activeChannels;
                    candidate(row) = [];
                    if ismember(chan, candidate)
                        obj.list.Data(row, column) = eventData.PreviousData;
                    else
                        obj.activeChannels(row) = chan;
                        evntData = updatedValuesEvntData('activeRecorderChannels');
                        notify(obj, 'updatedValues', evntData);
                    end
                end
                
                
            end
            
        end
        
    end
    
    
end