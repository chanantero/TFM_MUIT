function input_Callback(hObject, eventdata, handles)
% hObject    handle to input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sg_count;
tem = sg_count;
wave_input;
if sg_count ~= tem
    sg_name = ['signal',num2str(sg_count)];
    sg_info = {sg_name,0,0,true};
    filelist = get(handles.table,'data');
    if isempty(filelist)
        filelist = sg_info;
    else
        filelist(end + 1,:) = sg_info;
    end
    set(handles.table,'data',filelist);
end