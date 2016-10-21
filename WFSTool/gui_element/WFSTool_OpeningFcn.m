function WFSTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This is the initiation function of the software. 
% All the configuration set at the begining can be changed here. 
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Clear the table
filelist=[];
set(handles.table,'data',filelist);

% Initialize the axes
alt = generate_array();
plot(handles.axes1,alt(1,:),alt(2,:),'ob');
xlim(handles.axes1,[-2,6]);
ylim(handles.axes1,[-2,8]);
grid(handles.axes1);

global Fs;
Fs = 44100;

global pageSize;
pageSize = 4096;

% Initialize the timer
global readTimer;
readTimer = timer('Name','ReadTimer',...
'TimerFcn',@ReadTimerTask,...
'Period',0.5,...
'ExecutionMode','fixedspacing');
readTimer.UserData = handles;
start(readTimer);

global playTimer;
playTimer = timer('Name','PlayTimer',...
'TimerFcn',@PlayTimer,...
'Period',0.18,...
'ExecutionMode','fixedspacing');

global plotTimer;
plotTimer = timer('Name','ReadTimer',...
'StartFcn',@StartPlotTimerTask,...
'TimerFcn',@PlotTimerTask,...
'Period',0.5,...
'ExecutionMode','fixedspacing');
plotTimer.UserData = handles; 

% Set windowbuttondownfcn for click on the whole figure
set(gcf, 'WindowButtonDownFcn', 'tmouse down');

global isFinished;
isFinished = 0;

global sg_count signal_list;
sg_count = 0;
signal_list = {};

xlim(handles.ax_pattern,[-180,180]);
ylim(handles.ax_pattern,[-2,2]);

axes(handles.axes5);
ima = image_rotate(0,'1.jpg');
imshow(ima);

axes(handles.axes6);
ima2 = image_rotate(0,'2.jpg');
imshow(ima2);