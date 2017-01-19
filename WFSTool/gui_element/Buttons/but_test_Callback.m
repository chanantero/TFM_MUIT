function but_test_Callback(hObject, eventdata, handles)
% Original code: YU Bofei. Part of her Undergraduate Project "Software tool
% for Spatial Audio reproduction systems"
% Modified by Rubén Chuliá Mena (rchulia@outlook.com) on 18/01/2016.

% Instead of using playrec, I'm goint to use the Matlab DSP Toolbox

% hObject    handle to but_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% test_wav = waveForm();
global Fs;
test_wav = 0.8 * sin((1:Fs)*2*pi*1000/Fs)';
device_select = get(handles.menu_device,'Value');
[playDeviceID,chanList] = get_device(device_select);
if ~isempty(playDeviceID)
    % Init playrec with specific sample frequency, device and channels
    deviceWriter = audioDeviceWriter('Device', playDeviceID, 'SampleRate', Fs);
    fprintf('Test start! \n');
    for i = 1:size(chanList,2)
        fprintf('Test channel: %d\n',chanList(i));
        signal = zeros(numel(test_wav), max(chanList));
        signal(:, i) = test_wav;
        play(deviceWriter, signal);
    end
    pause(2)
    release(deviceWriter)
end

% % hObject    handle to but_test (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % test_wav = waveForm();
% global Fs;
% test_wav = 0.8 * sin((1:Fs)*2*pi*1000/Fs)';
% device_select = get(handles.menu_device,'Value');
% [playDeviceID,chanList] = select_play(device_select);
% if playDeviceID ~= -1
%     % Init playrec with specific sample frequency, device and channels
%     device_init(Fs,playDeviceID,chanList);
%     fprintf('Test start! \n');
%     for i = 1:size(chanList,2)
%         fprintf('Test channel: %d\n',chanList(i));
%         playrec('play',test_wav,chanList(i));
%     end
% end