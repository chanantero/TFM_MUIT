function but_playall_Callback(hObject,eventdata,handles)
global playTimer endSample startSample isTable music;
global playDeviceID chanList pageSize Fs endPoint namelist;

% Get info about which device is selected for every time user click on
% button 'play all'. Device is not allowed to change until the music is
% finished
device_select = get(handles.menu_device,'Value');
[playDeviceID,chanList] = select_play(device_select);

% Get info about page size for every time user click on button 'play all'.
% Page size is not allowed to change during the playing process
val = get(handles.menu_buffersize,'Value');
str = get(handles.menu_buffersize,'String');
str = cell2mat(str(val));
pageSize = str2num(str);
pageSize = check_pageSize(pageSize);

%if isTable == -1
%    errorpage('Empty Table! ');
%    return; 
%end
%if isTable == 0
%    errorpage('No active sound source!\n');
%    return;
%end

if playDeviceID ~= -1
    device_init(Fs,playDeviceID,chanList); % Initialize playrec configuration
    [endPoint,music] = read_music_file(namelist); % Read audio files
    t = 2 * size(music,2) * pageSize / Fs;
    t = round(t,3) - 0.14 * size(music,2);
    set(playTimer,'period',t);
    startSample = 1 + pageSize;
    endSample = startSample + pageSize - 1;
    start(playTimer);
end
fprintf('Click PlayAll button. \n');
