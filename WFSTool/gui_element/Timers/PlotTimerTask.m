function [] = PlotTimerTask(obj,event)
% This function is used to read the current position of the mouse and put
% the coordinate into GUI table
global x y selected parray_act an tn;
handles = obj.UserData;
list = get(handles.table,'data');

if selected > 0
    list(selected,2) = num2cell(x);
    list(selected,3) = num2cell(y);
end
set(handles.table,'data',list);
