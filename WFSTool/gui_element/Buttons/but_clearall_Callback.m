function but_clearall_Callback(hObject, eventdata, handles)
% this function used to clear the content of the table
global sg_count signal_list;
sg_count = 0;
data = {};
signal_list = {};
set(handles.table,'Data',data);
fprintf('Click ClearAll button. \n');