function [playDeviceID,chanList] = select_play(device_select)
% Original code: YU Bofei. Part of her Undergraduate Project "Software tool
% for Spatial Audio reproduction systems"
% Modified by Rubén Chuliá Mena (rchulia@outlook.com) on 22/10/2016.

a = 0;
device_list = playrec('getDevices');

if ischar(device_select)
    str = device_select;
elseif isnumeric(device_select)
    switch device_select
        case 1
%             str = '¶ú»ú (Cirrus Logic CS4208 (AB 1';
            str = 'Asignador de sonido Microsoft - Output';
        otherwise
            str = 'MOTU PCI ASIO';
    end
else
    error('select_play:wrongInput', 'Input must be of type char or numeric') 
end

 for i = 1:size(device_list,2);
    if strcmp(device_list(i).name,str)
        playDeviceID = device_list(i).deviceID;
        chanList = 1:device_list(i).outputChans;
        a = 1;
        break;
    end
end

if a == 0
    fprintf('Cannot find selected devices!\n');
    h = dialog('name','Error! ','position',[640,480,200,70]);  
    uicontrol('parent',h,'style','text','string','Cannot find selected devices! ','position',[15 10 180 50],'fontsize',10);  
    uicontrol('parent',h,'style','pushbutton','string','OK','position',[75 5 50 30],'callback','delete(gcbf)');
    playDeviceID = -1;
    chanList = [];
end
end

