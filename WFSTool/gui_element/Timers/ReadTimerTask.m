function [ ] = ReadTimerTask(obj,event)
global playDeviceID isTable isActive;
global an tn parray_act namelist angles dir;

% Pass GUI structure handles.
handles = obj.UserData;

% Renew uitable data. If the table is empty, call playrec to pause the
% music and stop reading the list
[namelist,vsx,vsy,isActive] = renew_vs(handles);

dir = str2num(get(handles.ed_angle,'string'));
if dir < 0 || dir > 360 || isempty(dir) || isnan(dir)
    dir = 0;
end

if isTable > 0
    % Calculate related coefficients based on sound source position
    [parray_act,an,tn] = WFS_parameters(playDeviceID,vsx,vsy);
    % Calculate the included angle of the source and the center of arrays
    angles = find_patter_angles(vsx,vsy,dir);
else
    parray_act = {};
end

plot_axes(parray_act,handles,vsx,vsy);
end


