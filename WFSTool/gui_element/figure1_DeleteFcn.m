function figure1_DeleteFcn(hObject, eventdata, handles)
% Function called when the main interface shut down
global readTimer playTimer;
stop(readTimer);
stop(playTimer);
fprintf('Window closed. \n');
if(playrec('isInitialised'))
    playrec('delPage');
    playrec('pause',1);
end