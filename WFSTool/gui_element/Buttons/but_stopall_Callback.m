function but_stopall_Callback(hObject, eventdata, handles)

global playTimer;
stop(playTimer);
playrec('pause',1);
playrec('delPage');
fprintf('Click StopAll button. \n');