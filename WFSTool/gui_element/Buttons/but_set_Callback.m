function but_set_Callback(hObject, eventdata, handles)

global pattern;
var_pat = get(handles.ed_pattern,'string');
cor_x = (-179:1:180);
tem = evalin('base',var_pat);
if size(tem,2) ~= 360
    errorpage('Invalid workspace variable!');
else
    plot(handles.ax_pattern,cor_x,tem,'-b');
    xlim(handles.ax_pattern,[-179,180]);
    ylim(handles.ax_pattern,[-max(tem)-1,max(tem)+1]);
end
pattern = [cor_x;tem];