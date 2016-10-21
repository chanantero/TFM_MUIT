function [] = StartPlotTimerTask(obj,event)
% This function is used to read the current position of the mouse and put
% the coordinate into GUI table
global x y selected;
handles = obj.UserData;
list = get(handles.table,'data');
row = size(list,1);
selected = 0;

if row ~= 0
    a = cell2mat(list(:,2));
    b = cell2mat(list(:,3));
    % Select the nearby source
    for i = 1:row
        if abs(a(i)-x) < 0.5 && abs(b(i)-y) < 0.5
            selected = i;
        end
    end
end

if selected > 0
    list(selected,2) = num2cell(x);
    list(selected,3) = num2cell(y);
end
set(handles.table,'data',list);
fprintf('Source selected: %d\n',selected);