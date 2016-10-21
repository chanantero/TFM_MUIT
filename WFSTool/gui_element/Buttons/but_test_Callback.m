function but_test_Callback(hObject, eventdata, handles)
% hObject    handle to but_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% test_wav = waveForm();
global Fs;
test_wav = 0.8 * sin((1:Fs)*2*pi*1000/Fs)';
device_select = get(handles.menu_device,'Value');
[playDeviceID,chanList] = select_play(device_select);
if playDeviceID ~= -1
    % Init playrec with specific sample frequency, device and channels
    device_init(Fs,playDeviceID,chanList);
    fprintf('Test start! \n');
    for i = 1:size(chanList,2)
        fprintf('Test channel: %d\n',chanList(i));
        playrec('play',test_wav,chanList(i));
    end
end