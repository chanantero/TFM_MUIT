classdef propertiesPanel < handle
    
    properties
        panel
        menuLabels
        frameDurationMenus
        deviceMenus
        driverMenus
    end
    
    methods
        function obj = propertiesPanel(reprodObj, fig, position, setFrameDurationFunc, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels)            
            obj.panel = uipanel(fig, 'Units', 'normalized', 'Position', position);
            
            if nargin == 7
                obj.setFunctions(setFrameDurationFunc, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels);
            end
            
            addlistener(reprodObj, 'playingState', 'PostSet', @(~, eventData) obj.playingStateSetCallback(eventData.AffectedObject.playingState));
        end
        
        function setFunctions(obj, setFrameDurationFunc, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels)
            % Delete previous menus
            delete(obj.frameDurationMenus);
            delete(obj.deviceMenus);
            delete(obj.driverMenus);
                
            % Create new ones
            N = numel(labels);
            obj.frameDurationMenus = gobjects(1, 1);
            obj.deviceMenus = gobjects(N, 1);
            obj.driverMenus = gobjects(N, 1);
            obj.menuLabels = gobjects(N, 1);
            
            % Positions
            yPos = (0:N-1)'/N;
            yDim = 1/N;
            xPos_label = 0.05;
            xDim_label = 0.15;
            xPos_frameDurationMenu = 0.3;
            xDim_frameDurationMenu = 0.15;
            xPos_deviceMenu = 0.55;
            xDim_deviceMenu = 0.15;
            xPos_driverMenu = 0.80;
            xDim_driverMenu = 0.15;
            
%             position_label = [xPos_label*ones(N, 1), yPos, xDim_label*ones(N, 1), yDim*ones(N, 1)];
%             position_frameDurationMenu = [xPos_frameDurationMenu*ones(N, 1), yPos, xDim_frameDurationMenu*ones(N, 1), yDim*ones(N, 1)];
%             position_deviceMenu = [xPos_deviceMenu*ones(N, 1), yPos, xDim_deviceMenu*ones(N, 1), yDim*ones(N, 1)];
%             position_driverMenu= [xPos_driverMenu*ones(N, 1), yPos, xDim_driverMenu*ones(N, 1), yDim*ones(N, 1)];
            
            obj.frameDurationMenus = uicontrol(obj.panel, 'Style', 'popupmenu', 'String', {'0.5', '1', '2'}, 'Tag', 'frameDuration',...
                    'Units', 'normalized', 'Position', [xPos_frameDurationMenu, yPos(1), xDim_frameDurationMenu, yDim], 'Callback', @(hObject, ~) setFrameDurationFunc(str2double(hObject.String{hObject.Value})));
             
            for k = 1:N         
                obj.deviceMenus(k) = uicontrol(obj.panel, 'Style', 'popupmenu', 'String', getAudioDevices(audioDeviceWriter()), 'Tag', 'device', 'Units', 'normalized', 'Position', [xPos_deviceMenu, yPos(k), xDim_deviceMenu, yDim],...
                    'Callback', @(hObject, ~) setDeviceFunc(hObject.String{hObject.Value}, k));
                
                obj.driverMenus(k) = uicontrol(obj.panel, 'Style', 'popupmenu', 'String', {'DirectSound', 'ASIO'}, 'Tag', 'driver', 'Units', 'normalized', 'Position', [xPos_driverMenu, yPos(k), xDim_driverMenu, yDim],...
                    'Callback', @(hObject, ~) driverMenuCallback(hObject.String{hObject.Value}, k));
                
                obj.menuLabels(k) = uicontrol(obj.panel, 'Style', 'text', 'String', labels{k}, 'Tag', 'label', 'Units', 'normalized', 'Position', [xPos_label, yPos(k), xDim_label, yDim]);
            end
            
            % Execute callbacks for the first time
            frameDuration = obj.frameDurationMenus.String(obj.frameDurationMenus.Value);
            setFrameDurationFunc(frameDuration);
            
            for k = 1:N
                driverOptions = obj.driverMenus(k).String;
                driver = driverOptions{obj.deviceMenus(k).Value};
                driverMenuCallback(driver, k);
                
                deviceOptions = obj.deviceMenus(k).String;
                device = deviceOptions{obj.deviceMenus(k).Value};
                setDeviceFunc(device, k);                
            end
            
            
            function driverMenuCallback(driver, ind)
                % Update deviceMenu
                setDriverFunc(driver, ind);
                
                % Get devices for the new driver
                obj.deviceMenus(ind).String = getAvailableDevicesFunc(ind);
                obj.deviceMenus(ind).Value = 1;
            end
            
        end
        
       
    end
    
    methods(Access = private)
        function playingStateSetCallback(obj, newPlayingState)
            if newPlayingState == playingStateClass('stopped')
                for k = 1:numel(obj.driverMenus)
                    obj.frameDurationMenus(k).Visible = 'on';
                    obj.deviceMenus(k).Visible = 'on';
                    obj.driverMenus(k).Visible = 'on';
                    obj.menuLabels(k).Visible = 'on';
                end
            else
                for k = 1:numel(obj.driverMenus)
                    obj.frameDurationMenus(k).Visible = 'off';
                    obj.deviceMenus(k).Visible = 'off';
                    obj.driverMenus(k).Visible = 'off';
                    obj.menuLabels(k).Visible = 'off';
                end
            end
        end
    end
    
end


