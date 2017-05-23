classdef propertiesPanel < handle
    
    properties
        panel
        menuLabels
        frameDurationMenu
        volumeControls
        deviceMenus
        driverMenus
    end
    
    methods
        function obj = propertiesPanel(fig, position, setFrameDurationFunc, setVolumeFunc, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels)            
            obj.panel = uipanel(fig, 'Units', 'normalized', 'Position', position);
            
            if nargin == 8
                obj.setFunctions(setFrameDurationFunc, setVolumeFunc, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels);
            end
            
        end
        
        function setFunctions(obj, setFrameDurationFunc, setVolumeFunc, setDeviceFunc, setDriverFunc, getAvailableDevicesFunc, labels)
            % Delete previous menus
            delete(obj.frameDurationMenu);
            delete(obj.volumeControls);
            delete(obj.deviceMenus);
            delete(obj.driverMenus);
            delete(obj.menuLabels);
                
            % Create new ones
            N = numel(labels);
            obj.frameDurationMenu = gobjects(1, 1);
            obj.volumeControls = gobjects(N, 1);
            obj.deviceMenus = gobjects(N, 1);
            obj.driverMenus = gobjects(N, 1);
            obj.menuLabels = gobjects(N, 1);
                        
            % Positions
            yPos = 0.8*(0:N-1)'/N;
            yDim = 0.8/max(3, N);
            xPos_label = 0.05;
            xDim_label = 0.15;
            xPos_Volume = 0.3;
            xDim_Volume = 0.15;
            xPos_deviceMenu = 0.55;
            xDim_deviceMenu = 0.15;
            xPos_driverMenu = 0.80;
            xDim_driverMenu = 0.15;
            
%             position_label = [xPos_label*ones(N, 1), yPos, xDim_label*ones(N, 1), yDim*ones(N, 1)];
%             position_frameDurationMenu = [xPos_frameDurationMenu*ones(N, 1), yPos, xDim_frameDurationMenu*ones(N, 1), yDim*ones(N, 1)];
%             position_deviceMenu = [xPos_deviceMenu*ones(N, 1), yPos, xDim_deviceMenu*ones(N, 1), yDim*ones(N, 1)];
%             position_driverMenu= [xPos_driverMenu*ones(N, 1), yPos, xDim_driverMenu*ones(N, 1), yDim*ones(N, 1)];

            obj.frameDurationMenu = uicontrol(obj.panel, 'Style', 'popupmenu', 'String', {'0.5', '1', '2'}, 'Tag', 'frameDuration',...
                    'Units', 'normalized', 'Position', [0.3, 0.9, 0.15, 0.1], 'Callback', @(hObject, ~) setFrameDurationFunc(str2double(hObject.String{hObject.Value})));
             
            for k = 1:N     
                obj.volumeControls(k) = uicontrol(obj.panel, 'Style', 'Edit', 'String', '1', 'Tag', 'Volume', 'Units', 'normalized', 'Position', [xPos_Volume, yPos(k), xDim_Volume, yDim],...
                    'Callback', @(hObject, ~) setVolumeFunc(str2double(hObject.String), k));
                
                obj.deviceMenus(k) = uicontrol(obj.panel, 'Style', 'popupmenu', 'String', getAudioDevices(audioDeviceWriter()), 'Tag', 'device', 'Units', 'normalized', 'Position', [xPos_deviceMenu, yPos(k), xDim_deviceMenu, yDim],...
                    'Callback', @(hObject, ~) setDeviceFunc(hObject.String{hObject.Value}, k));
                
                obj.driverMenus(k) = uicontrol(obj.panel, 'Style', 'popupmenu', 'String', {'DirectSound', 'ASIO'}, 'Tag', 'driver', 'Units', 'normalized', 'Position', [xPos_driverMenu, yPos(k), xDim_driverMenu, yDim],...
                    'Callback', @(hObject, ~) driverMenuCallback(hObject.String{hObject.Value}, k));
                
                obj.menuLabels(k) = uicontrol(obj.panel, 'Style', 'text', 'String', labels{k}, 'Tag', 'label', 'Units', 'normalized', 'Position', [xPos_label, yPos(k), xDim_label, yDim]);
            end
            
            % Execute callbacks for the first time
            frameDuration = str2double(obj.frameDurationMenu.String(obj.frameDurationMenu.Value));
            setFrameDurationFunc(frameDuration);
            
            for k = 1:N
                driverOptions = obj.driverMenus(k).String;
                driver = driverOptions{obj.driverMenus(k).Value};
                driverMenuCallback(driver, k);
                
                deviceOptions = obj.deviceMenus(k).String;
                device = deviceOptions{obj.deviceMenus(k).Value};
                setDeviceFunc(device, k);      
                
                volume = str2double(obj.volumeControls(k).String);
                setVolumeFunc(volume, k);
            end
            
            
            function driverMenuCallback(driver, ind)
                % Update deviceMenu
                setDriverFunc(driver, ind);
                
                % Get devices for the new driver
                obj.deviceMenus(ind).String = getAvailableDevicesFunc(ind);
                obj.deviceMenus(ind).Value = 1;
            end
            
        end
        
        function enableGUI(obj)
            for k = 1:numel(obj.driverMenus)
                obj.frameDurationMenu.Visible = 'on';
%                 obj.volumeControls(k).Visible = 'on';
                obj.deviceMenus(k).Visible = 'on';
                obj.driverMenus(k).Visible = 'on';
%                 obj.menuLabels(k).Visible = 'on';
            end
        end
        
        function disableGUI(obj)
            for k = 1:numel(obj.driverMenus)
                obj.frameDurationMenu.Visible = 'off';
%                 obj.volumeControls(k).Visible = 'off';
                obj.deviceMenus(k).Visible = 'off';
                obj.driverMenus(k).Visible = 'off';
%                 obj.menuLabels(k).Visible = 'off';
            end
        end
       
    end
    
    
end


