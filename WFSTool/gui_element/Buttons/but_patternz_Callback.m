function but_patternz_Callback(hObject, eventdata, handles)
global z_pat;
var_pat = get(handles.ed_patternz,'string');
cor_x = (-90:1:90);
tem = evalin('base',var_pat);
if size(tem,2) ~= 181
    errorpage('Invalid workspace variable!');
else
    plot(handles.axes7,cor_x,tem,'-b');
    xlim(handles.axes7,[-90,90]);
    ylim(handles.axes7,[-max(tem)-1,max(tem)+1]);
end
z_pat = [cor_x;tem];

