function [playDeviceID,chanList] = select_play(device_select)
a = 0;
device_list = playrec('getDevices');

if device_select == 1
    str = '¶ú»ú (Cirrus Logic CS4208 (AB 1';
else
    str = 'MOTU PCI ASIO';
end

 for i = 1:size(device_list,2);
    if strcmp(device_list(i).name,str)
        playDeviceID = device_list(i).deviceID
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

