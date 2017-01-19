function [ playDeviceID, chanList ] = get_device( device_select )
% Created by Rubén Chuliá Mena (18/01/2017)
deviceWriter = audioDeviceWriter;
device_list = deviceWriter.getAudioDevices;

if ischar(device_select)
    str = device_select;
elseif isnumeric(device_select)
    switch device_select
        case 1
            str = 'Default';
        otherwise
            str = 'MOTU PCI ASIO';
    end
else
    error('select_play:wrongInput', 'Input must be of type char or numeric') 
end

flag = ismember(str, device_list);
if flag == 0
    fprintf('Cannot find selected devices!\n');
    h = dialog('name','Error! ','position',[640,480,200,70]);  
    uicontrol('parent',h,'style','text','string','Cannot find selected devices! ','position',[15 10 180 50],'fontsize',10);  
    uicontrol('parent',h,'style','pushbutton','string','OK','position',[75 5 50 30],'callback','delete(gcbf)');
    playDeviceID = '';
    chanList = [];
else
    playDeviceID = str;
    deviceWriter.Device = playDeviceID;
    inf = info(deviceWriter);
    chanList = 1:inf.MaximumOutputChannels;
end



end

