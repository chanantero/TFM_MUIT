function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename] = uigetfile('.wav','Open wav file');
if isempty(strfind(filename,'wav')) == 1
    fprintf('Invalid Source File! \n');
    h = dialog('name','Error! ','position',[640,480,200,70]);  
    uicontrol('parent',h,'style','text','string','Invalid Source File! ','position',[15 10 180 50],'fontsize',10);  
    uicontrol('parent',h,'style','pushbutton','string','OK','position',[75 5 50 30],'callback','delete(gcbf)');
    return;  
end
filelist = get(handles.table,'data');
if isempty(filelist)
    filelist = {filename 0 0 true};
else
    filelist(end + 1,:) = {filename 0 0 true};
end
set(handles.table,'data',filelist);