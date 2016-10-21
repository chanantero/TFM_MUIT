function but_faceset_Callback(hObject, eventdata, handles)
global dir z_pat amp_z;

updown = str2num(get(handles.ed_z,'string'));

amp_z = 1;
for o = 1:size(z_pat,2)
    if z_pat(1,o) == updown
        amp_z = z_pat(2,o);
        break;
    end
end

if updown < 0
    updown = 360 - abs(updown);
end

ima1 = image_rotate(dir,'1.jpg');
ima2 = image_rotate(updown,'2.jpg');

axes(handles.axes5);
imshow(ima1);

axes(handles.axes6);
imshow(ima2);

